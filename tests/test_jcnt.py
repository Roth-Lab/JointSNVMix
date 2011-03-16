'''
Created on 2011-03-15

@author: andrew
'''
import os
import shutil

import unittest
from joint_snv_mix.file_formats.jcnt import JointCountsWriter, JointCountsReader


class Test(unittest.TestCase):
    def setUp(self):        
        if not os.path.exists('tmp'):
            os.mkdir('tmp')
            
        self.file_path = 'tmp/test.jcnt'
    
    def tearDown(self):        
        if os.path.exists('tmp'):
            shutil.rmtree('tmp')
    
    def test_add_row(self):
        writer = JointCountsWriter(self.file_path)
        
        chrom = '1'
        row = [12345, 'A', 'C', 'C', 100, 0, 100, 100]
        writer.add_row(chrom, row)        
        
        writer.close()
                
        reader = JointCountsReader(self.file_path)        
        size = reader.get_number_of_table_rows('1')
        reader.close()
        
        self.assertEqual(1, size)
        
    def test_add_row_using_context_manager(self):
        with JointCountsWriter(self.file_path) as writer:
            chrom = '1'
            row = [12345, 'A', 'C', 'C', 100, 0, 100, 100]
            writer.add_row(chrom, row)
        
        reader = JointCountsReader(self.file_path)        
        size = reader.get_number_of_table_rows('1')
        reader.close()
        
        self.assertEqual(1, size)
        
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
