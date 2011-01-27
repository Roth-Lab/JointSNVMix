#!/usr/bin/env/python
import sys
from pyleup.file_formats.jcnt import convert_bpils_to_jcnt

bpil_file_name = "../tests/data/ex1.bpil"
jcnt_file_name = "../tests/data/ex1.jcnt"

convert_bpils_to_jcnt( bpil_file_name, bpil_file_name, jcnt_file_name, 5, 30 )

