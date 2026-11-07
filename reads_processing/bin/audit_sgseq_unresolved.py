#!/usr/bin/env python3
"""Read-only second-stage audit of SGSeq source TSV exports.

Input: headerless sample<TAB>MAJIQ_dir<TAB>SGSeq_dir manifest.
Only reads sgseq_MXE.tsv, sgseq_A3SS.tsv, sgseq_A5SS.tsv,
sgseq_AFE.tsv, sgseq_ALE.tsv. Does not open SQLite or RDS.
Writes TSV/JSON reports into a *new*, previously nonexistent output directory.

Important: AFE/ALE rows have no variantID in the current exported TSVs;
identical PSI/assay values are evidence of repeated export measurements but
are NOT sufficient to establish biological alternative identity.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

MISSING = {'', 'NA', 'NaN', 'nan', 'NULL', 'None', '<NA>'}
FILES = {'MXE':'sgseq_MXE.tsv','A3SS':'sgseq_A3SS.tsv','A5SS':'sgseq_A5SS.tsv',
         'AFE':'sgseq_AFE.tsv','ALE':'sgseq_ALE.tsv'}

def clean(v):
    s = '' if v is None else str(v).strip()
    return '' if s in MISSING else s

def finite(v):
    s=clean(v)
    if not s: return None
    try:
        v=float(s)
        return v if math.isfinite(v) else None
    except ValueError: return None

def coordinate(row):
    chrom=clean(row.get('seqnames')) or clean(row.get('seqID'))
    strand=clean(row.get('strand'))
    if not chrom or strand not in {'+','-'}: return ''
    start=finite(row.get('junc_start') or row.get('start'))
    end=finite(row.get('junc_end') or row.get('end'))
    if start is None or end is None or not (start.is_integer() and end.is_integer() and 0 <= start < end): return ''
    return f'{chrom}:{int(start)}-{int(end)}:{strand}'

def hash_signature(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()[:24]

def tabulate(path,records,fields):
    with path.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader(); w.writerows(records)

def manifest_rows(path):
    with path.open(encoding='utf-8',newline='') as h:
        rows=[row for row in csv.reader(h,delimiter='\t') if row]
    if not rows or any(len(r)!=3 or not all(r) for r in rows):
        raise ValueError('Manifest must have sample<TAB>MAJIQ_dir<TAB>SGSeq_dir without header')
    if len({r[0] for r in rows}) != len(rows): raise ValueError('Duplicate samples in manifest')
    return rows

def diagnostic(rows,reference,output,example_limit):
    if output.exists(): raise ValueError(f'Output already exists: {output}')
    if example_limit < 0: raise ValueError('example-limit must be nonnegative')
    totals=defaultdict(Counter)
    mxe_groups=defaultdict(lambda:{'rows':0,'resolved':0,'genes':Counter(),'roles':Counter(),
        'featureID':Counter(),'featureID5p':Counter(),'featureID3p':Counter(),
        'raw_eventID':Counter(),'raw_variantType':Counter(),'first_rows':[]})
    mxe_examples=[]
    alternatives=defaultdict(lambda:{'rows':[],'positions':set()})
    graph_events=defaultdict(lambda:{'coords':defaultdict(set),'alternatives':set(),'examples':[], 'row_count':0})
    missing_files=[]
    for sample,_,directory in rows:
        for typ,basename in FILES.items():
            filename=Path(directory)/basename
            if not filename.is_file():
                missing_files.append(f'{sample}:{filename}')
                continue
            totals[typ]['files']+=1
            with filename.open(encoding='utf-8-sig',newline='') as h:
                reader=csv.DictReader(h,delimiter='\t')
                fields=reader.fieldnames or []
                if not fields: raise ValueError(f'Empty TSV: {filename}')
                for line,r in enumerate(reader,2):
                    c=totals[typ]; c['rows']+=1
                    eid=clean(r.get('eventID'))
                    gene=clean(r.get('geneID')) or clean(r.get('geneName'))
                    role=clean(r.get('event_type')) or clean(r.get('variantType'))
                    location=coordinate(r)
                    psi=finite(r.get('PSI'))
                    c['with_eventID']+=bool(eid)
                    c['with_coordinates']+=bool(location)
                    c['with_psi']+=psi is not None
                    if typ=='MXE':
                        # Unlike the previous audit, record the actual raw eventID,
                        # including missing or nonnumeric IDs, and per-sample group sizes.
                        groupkey=(sample,gene,eid)
                        g=mxe_groups[groupkey];g['rows']+=1;g['resolved']+=bool(location)
                        for field in ('featureID','featureID5p','featureID3p'):
                            g[field][clean(r.get(field)) or '<MISSING>']+=1
                        g['genes'][gene or '<MISSING>']+=1
                        g['roles'][role or '<MISSING>']+=1
                        g['raw_eventID'][eid or '<MISSING>']+=1
                        g['raw_variantType'][clean(r.get('variantType')) or '<MISSING>']+=1
                        if len(g['first_rows'])<2:
                            g['first_rows'].append(line)
                        if (not location or not eid) and len(mxe_examples)<example_limit:
                            mxe_examples.append(dict(sample=sample,line=line,eventID=eid,gene=gene,
                                variantType=clean(r.get('variantType')),featureID=clean(r.get('featureID')),
                                featureID5p=clean(r.get('featureID5p')),
                                featureID3p=clean(r.get('featureID3p')),
                                seqnames=clean(r.get('seqnames')),start=clean(r.get('start')),
                                end=clean(r.get('end')),strand=clean(r.get('strand')),
                                PSI=clean(r.get('PSI'))))
                    elif typ in ('A3SS','A5SS'):
                        frm=clean(r.get('from'));to=clean(r.get('to'))
                        if not (frm and to and eid):
                            c['unusable_graph_rows']+=1
                            continue
                        # Reconstruct exactly the earlier audit's graph signature,
                        # retaining per-source alternative and resolved coordinates.
                        eventkey=(sample,typ,gene,eid)
                        e=graph_events[eventkey]; e['row_count']+=1
                        triple=(role,frm,to)
                        e['alternatives'].add(triple)
                        e['coords'][triple].add(location or '<UNRESOLVED>')
                        if len(e['examples'])<2: e['examples'].append(f'{basename}:{line}')
                    elif typ in ('AFE','ALE'):
                        if not eid:
                            c['rows_without_eventID']+=1
                            continue
                        key=(sample,typ,gene,eid)
                        a=alternatives[key]
                        # Keep each source row as an *alternative candidate*;
                        # eventID alone must NOT be used as a variant key.
                        record=dict(sample=sample,event_type=typ,gene=gene,eventID=eid,line=line,
                                    junction=location,featureID=clean(r.get('featureID')),
                                    PSI=clean(r.get('PSI')),
                                    countsVariant5p=clean(r.get('countsVariant5p')),
                                    countsVariant3p=clean(r.get('countsVariant3p')),
                                    valid_PSI=psi is not None)
                        a['rows'].append(record)
    # MXE summary and evidence. These numbers diagnose exported data, not original R object.
    mxe_table=[]
    for (sample,gene,eid),g in mxe_groups.items():
        mxe_table.append(dict(sample=sample,gene=gene,eventID=eid or '<MISSING>',rows=g['rows'],
            resolved=g['resolved'],missing=g['rows']-g['resolved'],
            distinct_featureID=len(g['featureID']),
            missing_featureID=g['featureID'].get('<MISSING>',0),
            top_featureID=json.dumps(g['featureID'].most_common(3)),
            top_featureID5p=json.dumps(g['featureID5p'].most_common(3)),
            top_featureID3p=json.dumps(g['featureID3p'].most_common(3)),
            native_variant_types=json.dumps(g['raw_variantType'].most_common(3)),
            example_lines=','.join(map(str,g['first_rows']))))
    mxe_table.sort(key=lambda r:-r['rows'])
    # Graph conflict diagnostics. Compare exact role + from/to graph among samples.
    graph_by_signature=defaultdict(list)
    for (sample,typ,gene,eid),event in graph_events.items():
        graph=sorted(event['alternatives'])
        graph_sig=hash_signature([reference,typ,graph])
        alternatives_as_junctions=sorted((role,co) for (role,frm,to),values in event['coords'].items() for co in values)
        structural_sig=hash_signature([reference,typ,alternatives_as_junctions])
        graph_by_signature[(typ,graph_sig)].append(dict(sample=sample,typ=typ,gene=gene,eid=eid,
            graph=graph,graph_sig=graph_sig,structural_sig=structural_sig,
            junctions=set(alternatives_as_junctions),event=event))
    discord=[]; classification=Counter();discord_groups=Counter()
    for (typ,graph_sig),members in graph_by_signature.items():
        sigs={x['structural_sig'] for x in members}
        if len(sigs)<2: continue
        discord_groups[typ]+=1
        all_junctions=set().union(*(m['junctions'] for m in members))
        common_junctions=set.intersection(*(m['junctions'] for m in members))
        missing_by_member=[all_junctions-m['junctions'] for m in members]
        # 'subset' means at least one exporter yields a strict superset;
        # all sets may overlap without nesting in mixed cases.
        sets=[m['junctions'] for m in members]
        nested=all(a<=b or b<=a for i,a in enumerate(sets) for b in sets[i+1:])
        kind=('nested_partial_export' if nested else
              'non_nested_junction_substitutions')
        classification[(typ,kind)]+=1
        for m,missing in zip(members,missing_by_member):
            discord.append(dict(event_type=typ,graph_signature=graph_sig,
                sample=m['sample'],gene=m['gene'],eventID=m['eid'],
                discrepancy_class=kind,n_occurrences=len(members),
                n_structures=len(sigs),n_this=len(m['junctions']),
                n_union=len(all_junctions),n_common=len(common_junctions),
                n_missing=len(missing),
                missing_vs_union=';'.join(f'{role}::{coord}' for role,coord in sorted(missing)),
                present=';'.join(f'{role}::{coord}' for role,coord in sorted(m['junctions'])),
                graph=json.dumps(m['graph'],separators=(',',':')),
                example_rows=';'.join(m['event']['examples'])))
    # AFE/ALE diagnostic: no variantID, cannot safely identify alternatives.
    # Examine PSI and counts/coordinates within native event, without false conflicts.
    afeale=[]; a_stats=defaultdict(Counter)
    for (sample,typ,gene,eid),o in alternatives.items():
        rr=o['rows']; s=a_stats[typ];s['events']+=1;s['rows']+=len(rr)
        psi_vals={finite(x['PSI']) for x in rr if x['valid_PSI']}
        coords={x['junction'] for x in rr if x['junction']}
        # Compare if identical coordinates have discordant PSI.
        by_coord=defaultdict(set)
        by_psi=defaultdict(set)
        by_psi_counts=defaultdict(set)
        for x in rr:
            if x['junction'] and x['valid_PSI']:
                by_coord[x['junction']].add(finite(x['PSI']))
            if x['valid_PSI']:
                by_psi[finite(x['PSI'])].add(x['junction'])
                by_psi_counts[finite(x['PSI'])].add((x['countsVariant5p'],x['countsVariant3p']))
        conflicting_coords={x:vals for x,vals in by_coord.items() if len(vals)>1}
        # Repeated coordinate+same PSI+counts may be repeat rows from source;
        # cannot infer it is necessarily one biological variant.
        full_observation_keys={(x['junction'],x['PSI'],x['countsVariant5p'],x['countsVariant3p']) for x in rr}
        repeated_rows=len(rr)-len(full_observation_keys)
        s['events_with_multiple_psi']+=len(psi_vals)>1
        s['events_with_identical_coordinate_discordance']+=bool(conflicting_coords)
        s['events_with_repeated_identical_observations']+=repeated_rows>0
        s['repeated_identical_observation_rows']+=repeated_rows
        s['events_with_psis_on_multiple_coordinates']+=any(len(v)>1 for v in by_psi.values())
        afeale.append(dict(sample=sample,event_type=typ,gene=gene,eventID=eid,
            rows=len(rr),distinct_coordinates=len(coords),distinct_psi=len(psi_vals),
            identical_coordinate_discordances=len(conflicting_coords),
            repeated_identical_observation_rows=repeated_rows,
            psi_shared_among_junctions=sum(len(v)>1 for v in by_psi.values()),
            psi_with_mixed_variant_counts=sum(len(v)>1 for v in by_psi_counts.values()),
            coordinate_conflicts=';'.join(f'{co}:{sorted(v)}' for co,v in sorted(conflicting_coords.items())),
            source_observations=json.dumps(rr[:6],separators=(',',':'))))
    summary={'reference':reference,'samples':[r[0] for r in rows],
             'source_rows':{k:dict(v) for k,v in totals.items()},
             'missing_files':missing_files,
             'mxe':{'groups':len(mxe_groups),'source_rows':totals['MXE']['rows'],
                    'resolved_rows':totals['MXE']['with_coordinates'],
                    'unresolved_rows':totals['MXE']['rows']-totals['MXE']['with_coordinates'],
                    'largest_groups':mxe_table[:10],
                    'interpretation':'TSV-level diagnosis only; original SGSeq RDS must be examined for exporter defects'},
             'a3_a5':{'discordant_graph_groups':dict(discord_groups),
                      'discordance_types':{f'{typ}:{kind}':n for (typ,kind),n in classification.items()},
                      'discordant_occurrences':len(discord)},
             'afe_ale':{k:dict(v) for k,v in a_stats.items()},
             'limitations':['Original eventID is not a stable cross-sample identity.',
                 'AFE/ALE source exports do not retain variantID. Equal PSI is not proof of the same variant.',
                 'A3/A5 graph equivalence is a candidate match and does not prove equal PSI denominator.',
                 'MXE TSV analysis cannot establish the cause of missing fields in the upstream SGSeq RDS.',
                 'No original TSV, SQLite, or production parser is modified.']}
    output.mkdir(parents=True)
    tabulate(output/'mxe_group_diagnostics.tsv',mxe_table,
             ['sample','gene','eventID','rows','resolved','missing','distinct_featureID','missing_featureID','top_featureID','top_featureID5p','top_featureID3p','native_variant_types','example_lines'])
    tabulate(output/'mxe_unresolved_examples.tsv',mxe_examples,
             ['sample','line','eventID','gene','variantType','featureID','featureID5p','featureID3p','seqnames','start','end','strand','PSI'])
    tabulate(output/'a3_a5_graph_discordance.tsv',discord,
             ['event_type','graph_signature','sample','gene','eventID','discrepancy_class','n_occurrences','n_structures','n_this','n_union','n_common','n_missing','missing_vs_union','present','graph','example_rows'])
    tabulate(output/'afe_ale_alternative_candidates.tsv',afeale,
             ['sample','event_type','gene','eventID','rows','distinct_coordinates','distinct_psi','identical_coordinate_discordances','repeated_identical_observation_rows','psi_shared_among_junctions','psi_with_mixed_variant_counts','coordinate_conflicts','source_observations'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return summary

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',required=True,type=Path)
    ap.add_argument('--genome-version',required=True)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--example-limit',type=int,default=30)
    a=ap.parse_args()
    try:
        if not a.genome_version.strip(): raise ValueError('Missing genome version')
        summary=diagnostic(manifest_rows(a.manifest),a.genome_version,a.output,a.example_limit)
    except (OSError,ValueError) as exc:
        ap.error(str(exc))
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()