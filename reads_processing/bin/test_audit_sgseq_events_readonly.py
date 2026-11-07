#!/usr/bin/env python3
"""Synthetic non-mutating smoke tests for audit_sgseq_events_readonly.py."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().with_name('audit_sgseq_events_readonly.py')

class AuditSGSeqTests(unittest.TestCase):
    def test_event_recurrence_and_missing_native_ids(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            manifest=[]
            for sample, ids in [('S1',['101','102']),('S2',['902'])]:
                sg=root/sample/'sg';sg.mkdir(parents=True)
                maj=root/sample/'maj';maj.mkdir()
                manifest.append(f'{sample}\t{maj}\t{sg}\n')
                with (sg/'sgseq_SE.tsv').open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=['seqnames','start','end','strand','eventID','event_type','PSI'],delimiter='\t')
                    w.writeheader()
                    for event_id in ids:
                        if event_id=='102':
                            # Separate structure
                            coords=[(300,400,'SE:S'),(300,350,'SE:I')]
                        else:
                            coords=[(100,250,'SE:S'),(100,150,'SE:I'),(200,250,'SE:I')]
                        for start,end,role in coords:
                            w.writerow(dict(seqnames='Chr1',start=start,end=end,strand='+',eventID=event_id,event_type=role,PSI='0.3'))
                    # Null event id is counted, not reconstructed as a separate event
                    w.writerow(dict(seqnames='Chr1',start=500,end=600,strand='+',eventID='',event_type='SE:S',PSI='0.5'))
            manifest_path=root/'manifest.tsv'
            manifest_path.write_text(''.join(manifest))
            out=root/'out'
            args=[sys.executable,str(SCRIPT),'--manifest',str(manifest_path),'--genome-version','Athaliana_447','--output',str(out)]
            p=subprocess.run(args,text=True,capture_output=True)
            self.assertEqual(p.returncode,0,p.stderr)
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['n_native_event_occurrences'],3)
            self.assertEqual(summary['event_type_statistics']['SE']['distinct_structural_signatures'],2)
            self.assertEqual(summary['event_type_statistics']['SE']['shared_structural_signatures'],1)
            self.assertEqual(summary['source_file_statistics']['sgseq_SE.tsv']['rows_without_native_id'],2)
            # Refuse existing directory, keep original inputs untouched
            p2=subprocess.run(args,text=True,capture_output=True)
            self.assertNotEqual(p2.returncode,0)

if __name__=='__main__':
    unittest.main()
