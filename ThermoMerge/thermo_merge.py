###############################################################################
# MIT License
#
# Copyright (c) 2023 ArisSgouros
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

import os
import sys
import math
import argparse

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('filename', type=str, help='Path of the lammps thermo file')
parser.add_argument('header_part', type=str, help='Part of the header')
parser.add_argument('-hide_last', type=int, default=0, help='Show value from last step.')

def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

if __name__ == "__main__":
  args = parser.parse_args()
  filename = args.filename
  header_part = args.header_part
  hide_last = args.hide_last

  # Set flags
  print_flag = False
  first_header = True
  skip_current_frame = False

  # print the header
  print( header_part)

  with open(filename) as openfileobject:
    for cur_line in openfileobject:

        # Check if this line starts with string. If yes don't print.
        if cur_line.split(None, 1):
           if not is_number(cur_line.split(None, 1)[0]):
              print_flag = False
        else:
           print_flag = False

        # Print except if skip_current_frame says otherwise.
        if print_flag:
           if not skip_current_frame:
              print( cur_line, end = '')
           else:
              skip_current_frame = False

        # If this line is a header start printing
        if header_part in cur_line:
           if first_header:
              first_header = False
              skip_current_frame = False
           else:
              if hide_last:
                 skip_current_frame = True
           print_flag = True

  sys.exit(0)
