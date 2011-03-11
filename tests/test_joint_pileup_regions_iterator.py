'''
Created on 2011-03-10

@author: andrew
'''
import unittest
from joint_snv_mix.pre_processing.bam_to_jcnt import JointPileupRegionIterator


class Test(unittest.TestCase):
    def setUp(self):
        self.start = 200
        self.end = 300
        
        self.normal_positions = range(self.start, self.end)
        self.tumour_positions = range(self.start, self.end, 2)                
        
        self.normal_iter = iter([PileupColumnFake(x) for x in self.normal_positions])        
        self.tumour_iter = iter([PileupColumnFake(x) for x in self.tumour_positions])          
     
    def test_region_contains(self):
        regions = [(self.start + 20, self.start + 50)]
        
        self.check_regions(regions)
        
    def test_multiple_regions(self):
        regions = [(self.start + 20, self.start + 30),
                   (self.start + 60, self.start + 75)]
        
        self.check_regions(regions)
        
    def test_region_on_start_boundary(self):
        regions = [(self.start, self.start + 10)]
        
        self.check_regions(regions)
        
    def test_region_on_end_boundary(self):
        regions = [(self.end - 20, self.end)]
        
        self.check_regions(regions)
     
    def test_region_before_start(self):
        regions = [(self.start - 90, self.start - 80)]
        
        self.check_regions(regions)
        
    def test_region_past_end(self):
        regions = [(self.end + 100, self.end + 200)]
        
        self.check_regions(regions)
        
    def test_positions_match(self):   
        region = (10, 20)
        joint_iter = JointPileupRegionIterator(self.normal_iter, self.tumour_iter, [region])        
        
        for normal_col, tumour_col in joint_iter:            
            self.assertEqual(normal_col.pos, tumour_col.pos)
            
    def check_regions(self, regions):
        real_joint_pos = set(self.normal_positions) & set(self.tumour_positions)
        
        real_pos = []
        for region in regions:        
            real_pos.extend([x for x in real_joint_pos if region[0] <= x <= region[1]])
                
        joint_iter = JointPileupRegionIterator(self.normal_iter, self.tumour_iter, regions)        
        
        pred_joint_pos = set()

        for normal_col, tumour_col in joint_iter:
            pred_joint_pos.add(normal_col.pos)
        
        print real_pos, pred_joint_pos
        self.assertItemsEqual(real_pos, pred_joint_pos)
             
            
            
class PileupColumnFake:
    def __init__(self, pos):
        self.pos = pos

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
