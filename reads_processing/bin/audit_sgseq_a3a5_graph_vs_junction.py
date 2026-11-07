#!/usr/bin/env python3
"""Read-only A3SS/A5SS junction reconstruction check from exported SGSeq TSVs.

Compares per-variant graph endpoints ('from','to') with exporter-assigned
junc_start/junc_end. This detects assignment mismatches, but only original
SGSeq RDS can establish which feature was actually selected upstream.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

def node(text):
    parts = str(text or '').split(':')
    if len(parts) != 4 or parts[0] not in ('D','A') or parts[3] not in ('+','-'):
        return None
    try:
        return parts[0],parts[1],int(parts[2]),parts[3]
    except ValueError:
        return None

def classify(row, typ):
    a,b=node(row.get('from')),node(row.get('to'))
    if not a or not b or a[1] != b[1] or a[3] != b[3] or a[0] != 'D' or b[0] != 'A':
        return 'invalid_graph'
    try:
        x,y=int(row['junc_start']),int(row['junc_end'])
    except (ValueError,KeyError,TypeError):
        return 'missing_exported_coordinate'
    if x>=y: return 'invalid_exported_coordinate'
    expected=(min(a[2],b[2]),max(a[2],b[2]))
    return 'matches_variant_graph' if (x,y)==expected else 'disagrees_with_variant_graph'

def run(manifest,out):
    if out.exists(): raise ValueError('Refusing to overwrite existing output directory')
    counts=defaultdict(Counter);records=[]; missing=[]
    with manifest.open(newline='') as h:
        entries=[r for r in csv.reader(h,delimiter='\t') if r]
    if not entries or any(len(x)!=3 for x in entries): raise ValueError('Expected headerless three-column manifest')
    for sample,_,directory in entries:
        for typ in ('A3SS','A5SS'):
            p=Path(directory)/('sgseq_'+typ+'.tsv')
            if not p.is_file(): missing.append(str(p));continue
            with p.open(newline='',encoding='utf-8-sig') as h:
                for idx,row in enumerate(csv.DictReader(h,delimiter='\t'),2):
                    status=classify(row,typ)
                    counts[typ][status]+=1
                    if status!='matches_variant_graph':
                        records.append(dict(sample=sample,type=typ,line=idx,eventID=row.get('eventID',''),
                            geneID=row.get('geneID',''),variantID=row.get('variantID',''),
                            variantType=row.get('variantType',''),from_node=row.get('from',''),
                            to_node=row.get('to',''),start=row.get('junc_start',''),
                            end=row.get('junc_end',''),classification=status))
    out.mkdir(parents=True)
    (out/'summary.json').write_text(json.dumps(dict(counts={k:dict(v) for k,v in counts.items()},
        missing_files=missing,unmatched_records=len(records)),indent=2)+'\n')
    cols=['sample','type','line','eventID','geneID','variantID','variantType','from_node','to_node','start','end','classification']
    with (out/'junction_graph_mismatches.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(records)
    print((out/'summary.json').read_text())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    args=p.parse_args();run(args.manifest,args.output)
