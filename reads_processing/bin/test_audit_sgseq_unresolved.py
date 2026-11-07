#!/usr/bin/env python3
"""Small synthetic regression suite for the read-only SGSeq second-stage audit."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT=Path(__file__).with_name('audit_sgseq_unresolved.py')

def tsv(path,cols,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)

class AuditTest(unittest.TestCase):
    def test_mxe_graph_discordance_afeale_and_nonoverwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);manifest=[]
            for sample in ['S1','S2']:
                sg=root/sample/'sg';sg.mkdir(parents=True)
                ma=root/sample/'maj';ma.mkdir()
                manifest.append(f'{sample}\t{ma}\t{sg}\n')
                # Synthetic malformed MXE export with preserved eventID, missing coordinate fields.
                tsv(sg/'sgseq_MXE.tsv',['geneID','eventID','variantType','featureID','seqnames','start','end','strand','PSI'],[
                    dict(geneID='G1',eventID='1',variantType='MXE',featureID='',seqnames='',start='',end='',strand='',PSI='.4'),
                    dict(geneID='G1',eventID='1',variantType='MXE',featureID='4',seqnames='Chr1',start='100',end='200',strand='+',PSI='.4')])
                # Identical graph triples, but different junctions in second sample.
                tsv(sg/'sgseq_A3SS.tsv',['eventID','geneID','event_type','from','to','seqnames','strand','junc_start','junc_end','PSI'],[
                    dict(eventID='42',geneID='G2',event_type='A3SS:P',**{'from':'D:Chr1:1:+','to':'A:Chr1:2:+'},seqnames='Chr1',strand='+',junc_start='10',junc_end='20',PSI='.2'),
                    dict(eventID='42',geneID='G2',event_type='A3SS:D',**{'from':'D:Chr1:1:+','to':'A:Chr1:3:+'},seqnames='Chr1',strand='+',junc_start='10',junc_end='30' if sample=='S1' else '35',PSI='.8')])
                tsv(sg/'sgseq_AFE.tsv',['geneID','eventID','seqnames','strand','start','end','featureID','PSI','countsVariant3p'],[
                    dict(geneID='G3',eventID='9',seqnames='Chr1',strand='+',start='100',end='200',featureID='1',PSI='.2',countsVariant3p='8'),
                    dict(geneID='G3',eventID='9',seqnames='Chr1',strand='+',start='100',end='300',featureID='2',PSI='.8',countsVariant3p='32')])
            mp=root/'manifest.tsv';mp.write_text(''.join(manifest))
            out=root/'out';cmd=[sys.executable,str(SCRIPT),'--manifest',str(mp),'--genome-version','Athaliana_447','--output',str(out)]
            run=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['mxe']['unresolved_rows'],2)
            self.assertEqual(summary['a3_a5']['discordant_graph_groups'],{'A3SS':1})
            self.assertEqual(summary['afe_ale']['AFE']['events_with_multiple_psi'],2)
            self.assertEqual(summary['afe_ale']['AFE']['events_with_identical_coordinate_discordance'],0)
            self.assertEqual(summary['afe_ale']['AFE']['repeated_identical_observation_rows'],0)
            self.assertTrue((out/'a3_a5_graph_discordance.tsv').exists())
            self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)
            self.assertEqual(len((sg/'sgseq_AFE.tsv').read_text().splitlines()),3)

if __name__=='__main__':unittest.main()