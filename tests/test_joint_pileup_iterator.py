'''
Created on 2011-03-10

@author: andrew
'''
import unittest
from joint_snv_mix.pre_processing.bam_to_jcnt import JointPileupIterator


class Test(unittest.TestCase):    
    def test_all_positions_returned(self):        
        normal_positions = range(100)
        tumour_positions = range(0, 101, 2)
        
        real_joint_pos = set(normal_positions) & set(tumour_positions) 
        
        normal_iter = iter([PileupColumnFake(x) for x in normal_positions])        
        tumour_iter = iter([PileupColumnFake(x) for x in tumour_positions])

        joint_iter = JointPileupIterator(normal_iter, tumour_iter)        
        
        pred_joint_pos = set()
        for normal_col, tumour_col in joint_iter:
            pred_joint_pos.add(normal_col.pos)
            
        self.assertEqual(real_joint_pos, pred_joint_pos)
        
    def test_positions_match(self):   
        normal_positions = range(100)
        tumour_positions = range(0, 101, 2)

        normal_iter = iter([PileupColumnFake(x) for x in normal_positions])        
        tumour_iter = iter([PileupColumnFake(x) for x in tumour_positions])
        
        joint_iter = JointPileupIterator(normal_iter, tumour_iter)        
        
        for normal_col, tumour_col in joint_iter:            
            self.assertEqual(normal_col.pos, tumour_col.pos)
             
            
            
class PileupColumnFake:
    def __init__(self, pos):
        self.pos = pos

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
