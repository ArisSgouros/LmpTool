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

import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('dump_file', type=str, help='Path of the lammps dump file')
parser.add_argument('factor', type=int, help='Decimation degree')

def IncreaseCycleIter(i,n):
  if i == n:
    return 1
  else:
    return i+1

if __name__ == "__main__":
  args = parser.parse_args()
  dump_file = args.dump_file
  factor = args.factor
  foo = open(dump_file, 'r')
  do_write = factor
  for line in foo:
    line = line.rstrip('\n')
    if 'TIMESTEP' in line:
      do_write = IncreaseCycleIter(do_write,factor)
    if do_write == 1:
      print(line)
  foo.close()
