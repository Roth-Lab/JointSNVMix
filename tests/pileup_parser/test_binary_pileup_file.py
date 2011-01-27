'''
Created on 2010-07-29

@author: Andrew Roth
'''
import unittest

from tables import openFile
from tables.table import Table

from pyleup.pileup.binary_pileup import BinaryPileupFile, EntryExistsError

class Test( unittest.TestCase ):
    def setUp( self ):
        self.file_name = 'tmp/binary_pileup_test.bpil'
        self.file = BinaryPileupFile( self.file_name, 'w' )

    def tearDown( self ):
        self.file.close()

    def test_create_entry( self ):
        chr_name = 'chr1'

        self.file.create_entry( chr_name )

        self.assertTrue( chr_name in self.file.entries )
        self.assertTrue( isinstance( self.file.get_entry( chr_name )[0], Table ) )
        self.assertTrue( isinstance( self.file.get_entry( chr_name )[1], Table ) )

    def test_create_entry_more_than_once_raises_exception( self ):
        chr_name = 'chr1'

        self.file.create_entry( chr_name )

        self.assertRaises( EntryExistsError, self.file.create_entry, chr_name )

    def test_file_structure_is_correct( self ):
        chr_name = "chr1"

        self.file.create_entry( chr_name )

        self.file.close()

        h5_file = openFile( self.file_name, 'r' )

        self.assertTrue( "pileup" in h5_file.root )

        self.assertTrue( chr_name in h5_file.root.pileup )

        chr_group = h5_file.getNode( "/pileup", chr_name )

        self.assertTrue( "index" in chr_group )
        self.assertTrue( "reads" in chr_group )

        h5_file.close()



