#!/usr/bin/env python3
"""Synthetic read-only regression tests for audit_sgseq_three_questions.py."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().with_name('audit_sgseq_three_questions.py')


def put(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(x.keys() for x in rows)))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        w.writeheader()
        w.writerows(rows)


class AuditTests(unittest.TestCase):
    def test_three_diagnostics_and_input_safety(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            manifest=[]
            for sample, alt in [('S1',200),('S2',210)]:
                sg=root/sample/'sg';sg.mkdir(parents=True)
                maj=root/sample/'majiq';maj.mkdir()
                manifest.append(f'{sample}\t{maj}\t{sg}\n')
                rows=[dict(eventID=1,geneID='G1',event_type='A3SS:D',
                           **{'from':'D:Chr1:100:+','to':'A:Chr1:300:+'},
                           seqnames='Chr1',strand='+',junc_start=100,junc_end=300,PSI=.4,variantID=1),
                      dict(eventID=1,geneID='G1',event_type='A3SS:P',
                           **{'from':'D:Chr1:100:+','to':'A:Chr1:200:+'},
                           seqnames='Chr1',strand='+',junc_start=100,junc_end=alt,PSI=.6,variantID=2)]
                put(sg/'sgseq_A3SS.tsv',rows)
                put(sg/'sgseq_MXE.tsv',[dict(eventID=4,geneID='G2',event_type='MXE',seqnames='',strand='',start='',end='',PSI=.5,featureID=99)]*3)
                # Inclusion PSI repeated on both inclusion junctions, not two measures
                put(sg/'sgseq_SE.tsv',[dict(eventID=3,geneID='G3',event_type='SE:I',seqnames='Chr1',strand='+',start=100,end=150,PSI=.7),
                    dict(eventID=3,geneID='G3',event_type='SE:I',seqnames='Chr1',strand='+',start=200,end=250,PSI=.7),
                    dict(eventID=3,geneID='G3',event_type='SE:S',seqnames='Chr1',strand='+',start=100,end=250,PSI=.3)])
            m=root/'manifest.tsv';m.write_text(''.join(manifest))
            out=root/'out'
            a=[sys.executable,str(SCRIPT),'--manifest',str(m),'--genome-version','Athaliana_447','--output',str(out)]
            p=subprocess.run(a,capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            import json
            j=json.loads((out/'summary.json').read_text())
            self.assertEqual(j['mxe']['event_groups'],2)
            self.assertEqual(j['mxe']['groups_with_zero_resolved_junctions'],2)
            self.assertEqual(j['a3_a5_graph_conflict_counts']['A3SS'],1)
            self.assertEqual(j['variant_psi_by_type']['SE']['extra_rows_due_to_repeated_variant_groups'],2)
            self.assertEqual(j['variant_psi_by_type']['SE'].get('groups_with_discordant_PSI',0),0)
            self.assertEqual(j['two_role_psi_sums']['SE']['approximately_sum_to_one'],2)
            again=subprocess.run(a,capture_output=True,text=True)
            self.assertNotEqual(again.returncode,0)
            self.assertTrue((root/'S1'/'sg'/'sgseq_A3SS.tsv').exists())


if __name__=='__main__':
    unittest.main()