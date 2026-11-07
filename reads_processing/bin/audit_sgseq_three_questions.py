#!/usr/bin/env python3
"""Read-only audit of SGSeq MXE export, A3/A5 structural conflicts, and variant PSI.

Reads a headerless sample<TAB>MAJIQ_dir<TAB>SGSeq_dir manifest and SGSeq TSVs.
Does not read SQLite, run SGSeq, or modify input files. Requires a NEW output directory.
Cross-sample event IDs are never treated as stable: within-sample IDs group source rows.
The graph-signature audit checks exact role/from/to triples, not biological equivalence.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

FILES = ['sgseq_MXE.tsv', 'sgseq_A3SS.tsv', 'sgseq_A5SS.tsv',
         'sgseq_SE.tsv', 'sgseq_RI.tsv', 'sgseq_AFE.tsv', 'sgseq_ALE.tsv']
MISSING = {'', 'NA', 'NaN', 'nan', 'NULL', 'None', '<NA>'}


def clean(x):
    x = '' if x is None else str(x).strip()
    return '' if x in MISSING else x


def number(x):
    x = clean(x)
    if not x:
        return None
    try:
        n = float(x)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def digest(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:24]


def coord(row):
    ch = clean(row.get('seqnames')) or clean(row.get('seqID'))
    strand = clean(row.get('strand'))
    if (not ch or not strand) and clean(row.get('from')):
        t = clean(row['from']).split(':')
        if len(t) == 4:
            ch, strand = ch or t[1], strand or t[3]
    a = clean(row.get('junc_start')) or clean(row.get('ri_junc_start')) or clean(row.get('start'))
    b = clean(row.get('junc_end')) or clean(row.get('ri_junc_end')) or clean(row.get('end'))
    try:
        if not (ch and strand in ('+', '-') and a and b):
            return ''
        a, b = int(float(a)), int(float(b))
        if 0 <= a < b:
            return f'{ch}:{a}-{b}:{strand}'
    except ValueError:
        pass
    return ''


def write_table(path, rows, cols):
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def load_manifest(path):
    with path.open(newline='', encoding='utf-8') as fh:
        rows = [r for r in csv.reader(fh, delimiter='\t') if r]
    if not rows or any(len(r) != 3 or not all(r) for r in rows):
        raise ValueError('Manifest must have sample<TAB>MAJIQ_dir<TAB>SGSeq_dir, no header')
    if len({r[0] for r in rows}) != len(rows):
        raise ValueError('Repeated sample in manifest')
    return rows


def audit(manifest, reference, output, max_examples):
    if output.exists():
        raise ValueError(f'Refusing existing output directory: {output}')
    if max_examples < 0:
        raise ValueError('max-examples must be nonnegative')
    samples = [x[0] for x in manifest]
    files_stats = defaultdict(Counter)
    events = {}
    variants = {}
    mxe_bad = []
    missing_files = []
    # Distinct eventID and group sizes diagnose the MXE 192k-row / six-event anomaly.
    for sample, _, sgdir in manifest:
        for filename in FILES:
            path = Path(sgdir) / filename
            if not path.is_file():
                missing_files.append(f'{sample}:{path}')
                continue
            typ = filename.removeprefix('sgseq_').removesuffix('.tsv')
            c = files_stats[filename]
            c['files'] += 1
            with path.open(newline='', encoding='utf-8-sig') as fh:
                r = csv.DictReader(fh, delimiter='\t')
                c['has_eventID_column'] += int('eventID' in (r.fieldnames or []))
                c['has_variantID_column'] += int('variantID' in (r.fieldnames or []))
                c['has_event_id_column'] += int('event_id' in (r.fieldnames or []))
                for line, row in enumerate(r, 2):
                    c['rows'] += 1
                    original_id = clean(row.get('eventID'))
                    derived_id = clean(row.get('event_id'))
                    eid = original_id or derived_id
                    gene = clean(row.get('geneID')) or clean(row.get('geneName'))
                    role = clean(row.get('event_type')) or clean(row.get('variantType'))
                    from_node, to_node = clean(row.get('from')), clean(row.get('to'))
                    junction = coord(row)
                    psi = number(row.get('PSI'))
                    for col in ('featureID', 'featureID5p', 'featureID3p', 'junc_start', 'junc_end', 'seqnames', 'start', 'end', 'strand'):
                        if clean(row.get(col)):
                            c[f'has_{col}'] += 1
                    if eid:
                        c['rows_with_native_eventID'] += 1
                    if psi is not None:
                        c['rows_with_PSI'] += 1
                    if junction:
                        c['rows_with_resolved_junction'] += 1
                    if from_node and to_node:
                        c['rows_with_from_to'] += 1
                    if filename == 'sgseq_MXE.tsv' and (not junction or not eid) and len(mxe_bad) < max_examples:
                        mxe_bad.append({'sample':sample,'file':filename,'line':line,
                            'eventID':eid,'gene':gene,'event_type':role,
                            'start':clean(row.get('start')),'end':clean(row.get('end')),
                            'seqnames':clean(row.get('seqnames')),'strand':clean(row.get('strand')),
                            'featureID':clean(row.get('featureID')),
                            'featureID5p':clean(row.get('featureID5p')),
                            'featureID3p':clean(row.get('featureID3p')),
                            'PSI':clean(row.get('PSI'))})
                    if not eid:
                        continue
                    key = (sample, typ, gene, eid)
                    e = events.setdefault(key, {'rows':0, 'graph':set(), 'junction_roles':set(),
                        'missing_junction_rows':0,'roles':set(), 'variant_keys':set(), 'first_line':line})
                    e['rows'] += 1
                    e['roles'].add(role)
                    if from_node and to_node:
                        e['graph'].add((role,from_node,to_node))
                    if junction:
                        e['junction_roles'].add((role,junction))
                    else:
                        e['missing_junction_rows'] += 1
                    # Native SGSeq variantID may refer to a variant rather than a unique
                    # quantification row; include role in identity and track exact values.
                    vid = clean(row.get('variantID')) or clean(row.get('variantName'))
                    if vid:
                        var_key = ('native_variantID',vid,role)
                    elif from_node and to_node:
                        var_key = ('graph_alternative',role,from_node,to_node)
                    else:
                        var_key = ('role_only',role)
                    vk = key + (var_key,)
                    v = variants.setdefault(vk, {'source_rows':0, 'psi_values':set(),
                        'junctions':set(), 'line':line, 'missing_psi_rows':0,
                        'read5':set(),'read3':set()})
                    v['source_rows'] += 1
                    if psi is None:
                        v['missing_psi_rows'] += 1
                    else:
                        v['psi_values'].add(psi)
                    if junction:
                        v['junctions'].add(junction)
                    for field, label in [('countsVariant5p','read5'),('countsVariant3p','read3')]:
                        value = number(row.get(field))
                        if value is not None:
                            v[label].add(value)
                    e['variant_keys'].add(vk)

    event_stats = defaultdict(Counter)
    mxe_groups = []
    graph_groups = defaultdict(list)
    for (sample,typ,gene,eid), ev in events.items():
        c = event_stats[typ]
        c['native_event_occurrences'] += 1
        c['source_rows'] += ev['rows']
        c['events_missing_some_junctions'] += ev['missing_junction_rows'] > 0
        c['events_with_no_resolved_junctions'] += not ev['junction_roles']
        c['events_with_multiple_roles'] += len(ev['roles']) > 1
        if typ == 'MXE':
            mxe_groups.append({'sample':sample,'gene':gene,'eventID':eid,
                'rows':ev['rows'],'distinct_roles':len(ev['roles']),
                'distinct_junctions':len({x[1] for x in ev['junction_roles']}),
                'unresolved_rows':ev['missing_junction_rows'],
                'variant_groups':len(ev['variant_keys'])})
        if typ in ('A3SS','A5SS') and ev['graph']:
            gs = digest([reference,typ,sorted(ev['graph'])])
            js = digest([reference,typ,sorted(ev['junction_roles'])])
            graph_groups[(typ,gs)].append({'sample':sample,'gene':gene,'eventID':eid,
                'junction_signature':js,'junction_roles':ev['junction_roles'],
                'graph':ev['graph'], 'source_rows':ev['rows'],
                'missing_junction_rows':ev['missing_junction_rows']})

    conflicts = []
    group_conflicts = Counter()
    for (typ,gs), group in graph_groups.items():
        different = {r['junction_signature'] for r in group}
        if len(different) <= 1:
            continue
        group_conflicts[typ] += 1
        union = set().union(*(r['junction_roles'] for r in group))
        for r in group:
            conflicts.append({'event_type':typ,'graph_signature':gs,'sample':r['sample'],
                'gene':r['gene'],'native_eventID':r['eventID'],'source_rows':r['source_rows'],
                'missing_junction_rows':r['missing_junction_rows'],
                'junction_signature':r['junction_signature'],
                'n_group_occurrences':len(group),
                'n_distinct_junction_exports':len(different),
                'junctions_in_this_export':';'.join('::'.join(v) for v in sorted(r['junction_roles'])),
                'junctions_missing_vs_group_union':';'.join('::'.join(v) for v in sorted(union-r['junction_roles'])),
                'graph_alternatives':json.dumps(sorted(r['graph']))})

    variant_stats = defaultdict(Counter)
    variant_conflicts=[]
    psi_examples=[]
    for (sample,typ,gene,eid,vkey),v in variants.items():
        st=variant_stats[typ]
        st['variant_groups'] += 1
        st['source_rows'] += v['source_rows']
        st['extra_rows_due_to_repeated_variant_groups'] += v['source_rows'] - 1
        st['groups_with_repeated_rows'] += v['source_rows'] > 1
        st['groups_with_discordant_PSI'] += len(v['psi_values']) > 1
        st['groups_with_missing_PSI'] += v['missing_psi_rows'] > 0
        st['groups_with_discordant_variant_5p_counts'] += len(v['read5']) > 1
        st['groups_with_discordant_variant_3p_counts'] += len(v['read3']) > 1
        st['groups_with_multiple_junctions'] += len(v['junctions']) > 1
        if v['source_rows'] > 1 and len(psi_examples) < max_examples:
            psi_examples.append({'sample':sample,'event_type':typ,'gene':gene,
                'native_eventID':eid,'variant_key':repr(vkey),'source_rows':v['source_rows'],
                'distinct_junctions':len(v['junctions']),'psi_values':';'.join(map(str,sorted(v['psi_values']))),
                'distinct_read5_values':len(v['read5']),'distinct_read3_values':len(v['read3'])})
        if len(v['psi_values']) > 1:
            variant_conflicts.append({'sample':sample,'event_type':typ,'gene':gene,
                'native_eventID':eid,'variant_key':repr(vkey),'source_rows':v['source_rows'],
                'distinct_junctions':len(v['junctions']),'psi_values':';'.join(map(str,sorted(v['psi_values'])))})
    # The exporter can copy PSI from a single variant to several junction rows.
    # Do not interpret total per-row PSI or per-junction PSI as independent variant calls.
    # Test sums only for events with exactly two distinct role-only groups with unambiguous PSI.
    complement = defaultdict(Counter)
    complement_examples=[]
    by_event = defaultdict(list)
    for k,v in variants.items():
        by_event[k[:4]].append((k[4],v))
    for (sample,typ,gene,eid),entries in by_event.items():
        if typ not in ('SE','A3SS','A5SS','RI'):
            continue
        role_psis=defaultdict(set)
        for key,v in entries:
            role = key[1] if key[0] in ('native_variantID','graph_alternative','role_only') else ''
            role_psis[role].update(v['psi_values'])
        if len(role_psis)!=2 or any(len(s)!=1 for s in role_psis.values()):
            continue
        a,b = [next(iter(s)) for s in role_psis.values()]
        complement[typ]['two_role_events_with_unambiguous_PSI'] += 1
        if math.isclose(a+b,1,rel_tol=1e-5,abs_tol=1e-3):
            complement[typ]['approximately_sum_to_one'] += 1
        else:
            complement[typ]['not_sum_to_one'] += 1
            if len(complement_examples)<max_examples:
                complement_examples.append({'sample':sample,'event_type':typ,'gene':gene,'eventID':eid,
                    'roles':';'.join(sorted(role_psis)),'psi_sum':a+b})

    output.mkdir(parents=True)
    write_table(output/'mxe_unresolved_examples.tsv',mxe_bad, ['sample','file','line','eventID','gene','event_type','start','end','seqnames','strand','featureID','featureID5p','featureID3p','PSI'])
    write_table(output/'mxe_event_groups.tsv', sorted(mxe_groups,key=lambda x:-x['rows']), ['sample','gene','eventID','rows','distinct_roles','distinct_junctions','unresolved_rows','variant_groups'])
    write_table(output/'a3_a5_graph_conflicts_detailed.tsv',conflicts, ['event_type','graph_signature','sample','gene','native_eventID','source_rows','missing_junction_rows','junction_signature','n_group_occurrences','n_distinct_junction_exports','junctions_in_this_export','junctions_missing_vs_group_union','graph_alternatives'])
    write_table(output/'variant_psi_conflicts.tsv',variant_conflicts,['sample','event_type','gene','native_eventID','variant_key','source_rows','distinct_junctions','psi_values'])
    write_table(output/'variant_repeated_examples.tsv',psi_examples,['sample','event_type','gene','native_eventID','variant_key','source_rows','distinct_junctions','psi_values','distinct_read5_values','distinct_read3_values'])
    write_table(output/'two_role_psi_sum_exceptions.tsv',complement_examples,['sample','event_type','gene','eventID','roles','psi_sum'])
    summary = {'reference':reference,'samples':samples,'missing_files':missing_files,
       'source_files':{k:dict(v) for k,v in sorted(files_stats.items())},
       'event_groups_by_type':{k:dict(v) for k,v in sorted(event_stats.items())},
       'mxe':{'event_groups':len(mxe_groups),'max_rows_per_native_event_group':max((x['rows'] for x in mxe_groups),default=0),
              'groups_with_zero_resolved_junctions':sum(x['distinct_junctions']==0 for x in mxe_groups),
              'groups_with_any_unresolved_rows':sum(x['unresolved_rows']>0 for x in mxe_groups)},
       'a3_a5_graph_conflict_counts':dict(group_conflicts),
       'a3_a5_graph_conflict_occurrences':len(conflicts),
       'variant_psi_by_type':{k:dict(v) for k,v in sorted(variant_stats.items())},
       'two_role_psi_sums':{k:dict(v) for k,v in sorted(complement.items())},
       'interpretation_limits':[
          'SGSeq eventID and variantID are used only within one sample, not across samples.',
          'Missing junctions can be a legitimate feature of some alternative definitions (e.g. RI), not automatically invalid biology.',
          'Multiple rows per event variant can be duplicated export of one quantification; row counts do not equal independent PSI observations.',
          'Role-level PSI sum is a diagnostic only, not a universal SGSeq invariant.',
          'MXE group sizes and coordinate validity must be checked against source exporter and preferably sgseq_variants.rds.',
          'Graph-boundary identity is not proof that two definitions have identical variant sets or PSI denominators.',
          'No input file or production database was modified.'
       ]}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--genome-version',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-examples',type=int,default=30)
    args=p.parse_args()
    try:
        audit(load_manifest(args.manifest),args.genome_version,args.output,args.max_examples)
    except (ValueError,OSError) as exc:
        p.error(str(exc))

if __name__=='__main__':
    main()