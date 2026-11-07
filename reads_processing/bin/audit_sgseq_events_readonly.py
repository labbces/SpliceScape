#!/usr/bin/env python3
"""Read-only audit of SGSeq event identities using original per-sample TSV exports.

Inputs: 3-column SpliceScape manifest: sample<TAB>MAJIQ_directory<TAB>SGSeq_directory.
Reads only SGSeq TSVs (never SQLite). Generates diagnostic files under --output.

SGSeq eventID/variantID are analysis-local; this audit constructs *candidate* cross-sample
structural identities using the exported event graph boundaries and alternative roles.
It does not assume PSI of different alternatives is a single event-wide number.

No input is modified. The output directory must not already exist.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys

FILES = {
    'sgseq_SE.tsv': 'SE', 'sgseq_RI.tsv': 'RI',
    'sgseq_A5SS.tsv': 'A5SS', 'sgseq_A3SS.tsv': 'A3SS',
    'sgseq_MXE.tsv': 'MXE', 'sgseq_AFE.tsv': 'AFE',
    'sgseq_ALE.tsv': 'ALE', 'sgseq_constitutive.tsv': 'NALT',
}
BAD = {'', 'NA', 'NaN', 'nan', 'None', 'NULL', '<NA>'}

def clean(value):
    x = str(value if value is not None else '').strip()
    return '' if x in BAD else x

def digest(parts):
    return hashlib.sha256(json.dumps(parts, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:24]

def numeric(value):
    value = clean(value)
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None

def valid_coord(row):
    chrom = clean(row.get('seqnames')) or clean(row.get('seqID'))
    strand = clean(row.get('strand'))
    if not strand:
        node = clean(row.get('from'))
        ss = node.split(':')
        if len(ss) == 4:
            strand = ss[3]
            chrom = chrom or ss[1]
    s = clean(row.get('junc_start')) or clean(row.get('ri_junc_start')) or clean(row.get('start'))
    e = clean(row.get('junc_end')) or clean(row.get('ri_junc_end')) or clean(row.get('end'))
    if not (chrom and strand in ('+', '-') and s and e):
        return ''
    try:
        a, b = int(float(s)), int(float(e))
    except ValueError:
        return ''
    if a < 0 or b <= a:
        return ''
    return f'{chrom}:{a}-{b}:{strand}'

def metadata_from_row(r):
    native = clean(r.get('eventID')) or clean(r.get('event_id'))
    gene = clean(r.get('geneID')) or clean(r.get('geneName'))
    role = clean(r.get('event_type')) or clean(r.get('variantType'))
    src_var = clean(r.get('variantID')) or clean(r.get('variantName'))
    fr, to = clean(r.get('from')), clean(r.get('to'))
    junc = valid_coord(r)
    return native, gene, role, src_var, fr, to, junc

def write_tsv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as h:
        w = csv.DictWriter(h, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', required=True, type=Path)
    ap.add_argument('--genome-version', required=True)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--max-examples', type=int, default=15)
    args = ap.parse_args()
    if args.output.exists():
        ap.error(f'Refusing existing output directory: {args.output}')
    with args.manifest.open(newline='', encoding='utf-8') as h:
        manifest = [r for r in csv.reader(h, delimiter='\t') if r]
    if not manifest or any(len(r) != 3 for r in manifest):
        ap.error('Manifest must contain sample<TAB>MAJIQ_dir<TAB>SGSeq_dir without header')
    if len({r[0] for r in manifest}) != len(manifest):
        ap.error('Duplicate sample IDs')
    if not args.genome_version.strip():
        ap.error('--genome-version cannot be empty')
    if args.max_examples < 0:
        ap.error('--max-examples must be nonnegative')

    # Each source event is grouped using (sample, class, original SGSeq eventID, gene).
    # Missing IDs are deliberately *not* joined into fake multi-row events.
    events = {}
    summary_by_file = defaultdict(Counter)
    header_report = []
    missing_files = []
    for sample, _, sg_dir in manifest:
        for filename, etype in FILES.items():
            path = Path(sg_dir) / filename
            if not path.is_file():
                missing_files.append(dict(sample=sample, filename=filename, path=str(path)))
                continue
            with path.open(newline='', encoding='utf-8-sig') as h:
                reader = csv.DictReader(h, delimiter='\t')
                cols = reader.fieldnames or []
                header_report.append(dict(sample=sample, filename=filename, columns=','.join(cols),
                                          has_eventID=int('eventID' in cols), has_event_id=int('event_id' in cols),
                                          has_variantType=int('variantType' in cols), has_from_to=int('from' in cols and 'to' in cols)))
                c = summary_by_file[filename]
                c['files'] += 1
                for line, row in enumerate(reader, start=2):
                    c['rows'] += 1
                    native, gene, role, src_var, fr, to, junction = metadata_from_row(row)
                    psi = numeric(row.get('PSI'))
                    if psi is not None:
                        c['rows_with_psi'] += 1
                    if native:
                        c['rows_with_native_id'] += 1
                    else:
                        c['rows_without_native_id'] += 1
                    if clean(row.get('eventID')):
                        c['rows_with_original_eventID'] += 1
                    if clean(row.get('event_id')):
                        c['rows_with_derived_event_id'] += 1
                    if junction:
                        c['rows_with_resolved_junction'] += 1
                    if fr and to:
                        c['rows_with_from_to'] += 1
                    if not native:
                        continue
                    key = (sample, etype, gene, native)
                    ev = events.setdefault(key, dict(sample=sample, event_type=etype, gene=gene,
                                  native_event_id=native, rows=0, source_files=set(),
                                  junctions=set(), graph_alternatives=set(), structural_alternatives=set(), roles=set(),
                                  role_psi=defaultdict(set), variant_psi=defaultdict(set),
                                  missing_junction_rows=0, missing_graph_rows=0,
                                  example_rows=[]))
                    ev['rows'] += 1
                    ev['source_files'].add(filename)
                    if junction:
                        ev['junctions'].add(junction)
                        ev['structural_alternatives'].add((role, junction))
                    else:
                        ev['missing_junction_rows'] += 1
                    if fr and to:
                        # Graph alternative signature ignores local featureID and native IDs.
                        # Multiplicity and variant roles are retained separately.
                        ev['graph_alternatives'].add((role, fr, to))
                    else:
                        ev['missing_graph_rows'] += 1
                    ev['roles'].add(role)
                    if psi is not None:
                        if junction:
                            ev['role_psi'][(role, junction)].add(psi)
                        if src_var:
                            ev['variant_psi'][(role, src_var)].add(psi)
                    if len(ev['example_rows']) < 2:
                        ev['example_rows'].append(f'{filename}:{line}')

    results = []
    grouped = defaultdict(list)
    node_grouped = defaultdict(set)
    structural_grouped = defaultdict(list)
    for ev in events.values():
        # Candidate event fingerprint includes role and graph endpoints, not local event IDs.
        # The RI retention alternative can have no true junction: it is not discarded.
        parts = sorted([list(t) for t in ev['graph_alternatives']])
        graph_ok = ev['missing_graph_rows'] == 0 and bool(parts)
        graph_sig = digest([args.genome_version, ev['event_type'], parts]) if graph_ok else ''
        junc_ok = ev['missing_junction_rows'] == 0 and bool(ev['junctions'])
        structural_sig = digest([args.genome_version, ev['event_type'], sorted([list(t) for t in ev['structural_alternatives']])]) if junc_ok else ''
        junction_sig = digest([args.genome_version, ev['event_type'], sorted(ev['junctions'])]) if junc_ok else ''
        variant_conflict = sum(len(v) > 1 for v in ev['variant_psi'].values())
        junction_role_conflict = sum(len(v) > 1 for v in ev['role_psi'].values())
        out = dict(sample=ev['sample'], event_type=ev['event_type'],gene=ev['gene'],
                   native_event_id=ev['native_event_id'], source_files=','.join(sorted(ev['source_files'])),
                   source_rows=ev['rows'], n_roles=len(ev['roles']),
                   roles=';'.join(sorted(ev['roles'])), n_graph_alternatives=len(parts),
                   graph_signature=graph_sig, n_junctions=len(ev['junctions']),
                   junction_signature=junction_sig, structural_signature=structural_sig, missing_graph_rows=ev['missing_graph_rows'],
                   missing_junction_rows=ev['missing_junction_rows'],
                   variant_psis_disagree=variant_conflict,
                   junction_role_psis_disagree=junction_role_conflict,
                   graph_alternatives=json.dumps(parts,separators=(',',':')),
                   junction_coordinates=';'.join(sorted(ev['junctions'])))
        results.append(out)
        if graph_sig:
            grouped[(ev['event_type'],graph_sig)].append(out)
        if structural_sig:
            structural_grouped[(ev['event_type'],structural_sig)].append(out)
        if junction_sig:
            node_grouped[(ev['event_type'],junction_sig)].add(ev['sample'])

    for ev in results:
        struct_matched = structural_grouped.get((ev['event_type'],ev['structural_signature']),[]) if ev['structural_signature'] else []
        ev['samples_with_same_structure'] = len({v['sample'] for v in struct_matched})
        matched = grouped.get((ev['event_type'],ev['graph_signature']),[]) if ev['graph_signature'] else []
        ev['samples_with_same_graph'] = len({v['sample'] for v in matched})
        ev['junction_signatures_with_same_graph'] = len({v['junction_signature'] for v in matched})
        ev['native_ids_with_same_graph'] = len({v['native_event_id'] for v in matched})

    recurrence = defaultdict(Counter)
    by_type = defaultdict(Counter)
    graph_conflicts=[]
    for (etype,sig), group in structural_grouped.items():
        ss={r['sample'] for r in group}
        recurrence[etype]['distinct_structural_signatures'] += 1
        recurrence[etype]['shared_structural_signatures'] += len(ss)>1
        recurrence[etype]['structural_signatures_in_all_samples'] += len(ss)==len(manifest)
    for (etype,sig), group in grouped.items():
        samples = {r['sample'] for r in group}
        c = recurrence[etype]
        c['distinct_graph_signatures'] += 1
        if len(samples)>1:
            c['shared_graph_signatures'] += 1
        if len(samples)==len(manifest):
            c['graph_signatures_in_all_samples'] += 1
        junctions={r['junction_signature'] for r in group}
        if len(junctions)>1:
            c['graph_signatures_with_discordant_junction_export']+=1
            if len(graph_conflicts)<args.max_examples:
                graph_conflicts.append(dict(event_type=etype,graph_signature=sig,
                                            n_samples=len(samples),junction_signatures=';'.join(sorted(junctions)),
                                            examples=';'.join(f"{r['sample']}:{r['native_event_id']}" for r in group[:5])))
    for r in results:
        c=by_type[r['event_type']]
        c['event_occurrences']+=1
        c['event_occurrences_with_graph']+=bool(r['graph_signature'])
        c['event_occurrences_with_complete_structural_signature']+=bool(r['structural_signature'])
        c['event_occurrences_with_complete_junctions']+=bool(r['junction_signature'])
        c['occurrences_with_multiple_variant_PSI_values']+=bool(r['variant_psis_disagree'])
        c['occurrences_with_multiple_role_junction_PSI_values']+=bool(r['junction_role_psis_disagree'])
    for etype, cs in node_grouped.items():
        pass
    for (etype,_), samples in node_grouped.items():
        recurrence[etype]['distinct_complete_junction_signatures']+=1
        recurrence[etype]['shared_complete_junction_signatures']+=len(samples)>1

    args.output.mkdir(parents=True)
    write_tsv(args.output/'source_headers.tsv',header_report,['sample','filename','columns','has_eventID','has_event_id','has_variantType','has_from_to'])
    fields=['sample','event_type','gene','native_event_id','source_files','source_rows','n_roles','roles','n_graph_alternatives','graph_signature','n_junctions','junction_signature','structural_signature','missing_graph_rows','missing_junction_rows','variant_psis_disagree','junction_role_psis_disagree','samples_with_same_graph','samples_with_same_structure','junction_signatures_with_same_graph','native_ids_with_same_graph','graph_alternatives','junction_coordinates']
    write_tsv(args.output/'sgseq_event_occurrences.tsv',sorted(results,key=lambda r:(r['event_type'],r['sample'],r['gene'],r['native_event_id'])),fields)
    write_tsv(args.output/'graph_export_conflicts.tsv',graph_conflicts,['event_type','graph_signature','n_samples','junction_signatures','examples'])
    summary=dict(reference=args.genome_version, samples=[r[0] for r in manifest],
                 n_native_event_occurrences=len(results), missing_files=missing_files,
                 source_file_statistics={k:dict(v) for k,v in sorted(summary_by_file.items())},
                 event_type_statistics={k:dict(by_type[k] | recurrence[k]) for k in sorted(by_type)},
                 caveats=[
                   'Native eventID is used only to group rows within a sample, never for cross-sample matching.',
                   'Graph signature uses exact event type and role/from/to boundary triples, excluding native identifiers.',
                   'Structural signature uses event type and the COMPLETE role-labeled junction set; this is the preferred candidate match for SE/MXE/AFE/ALE exports without from/to.',
                   'Matching graph signatures are candidate equivalent definitions; shared boundaries alone do not prove equivalent denominator or isoform topology.',
                   'SE inclusion export duplicates one SGSeq variant PSI across multiple constituent junction rows; do not sum these rows as independent measurements.',
                   'Missing/invalid junction coordinates do not invalidate graph-defined events, notably RI retained alternatives.',
                   'A3SS/A5SS source exports include variant-level records with eventID but have no synthetic event_id column.',
                   'Partial event graphs and sample-specific variant discovery can prevent exact signatures from recurring.',
                 ])
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    print('Outputs:',args.output/'summary.json',args.output/'source_headers.tsv',args.output/'sgseq_event_occurrences.tsv',args.output/'graph_export_conflicts.tsv')
    return 0

if __name__=='__main__':
    sys.exit(main())