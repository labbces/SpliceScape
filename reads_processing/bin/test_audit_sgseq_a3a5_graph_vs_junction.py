#!/usr/bin/env python3
"""Synthetic checks for graph-to-junction assignment, no R required."""
import unittest
from audit_sgseq_a3a5_graph_vs_junction import classify

class ReconstructionTests(unittest.TestCase):
    def row(self, fr, to, st, en):
        return {'from':fr,'to':to,'junc_start':str(st),'junc_end':str(en)}
    def test_plus_correct_and_wrong_extreme(self):
        r=self.row('D:Chr1:100:+','A:Chr1:200:+',100,200)
        self.assertEqual(classify(r,'A3SS'),'matches_variant_graph')
        r['junc_end']='250'
        self.assertEqual(classify(r,'A3SS'),'disagrees_with_variant_graph')
    def test_minus(self):
        r=self.row('D:Chr2:500:-','A:Chr2:300:-',300,500)
        self.assertEqual(classify(r,'A5SS'),'matches_variant_graph')
    def test_missing_and_invalid(self):
        r=self.row('D:Chr1:100:+','A:Chr1:200:+','','')
        self.assertEqual(classify(r,'A3SS'),'missing_exported_coordinate')
        r=self.row('D:Chr1:100:+','A:Chr2:200:+',100,200)
        self.assertEqual(classify(r,'A3SS'),'invalid_graph')
if __name__=='__main__': unittest.main()
