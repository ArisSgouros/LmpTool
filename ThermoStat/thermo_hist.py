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
import math as m
import numpy as np

parser = argparse.ArgumentParser(description='Statistics of thermo quantities')
parser.add_argument('thermo_file', type=str, help='Path of the lammps thermo file')
parser.add_argument('header', type=str, help='header')
parser.add_argument('end_flag', type=str, help='end_flag')
parser.add_argument('-skip', type=int, default=0, help='skip lines')
parser.add_argument('-nbin', type=int, default=10, help='histogram bins')


args = parser.parse_args()
FILENAME          = args.thermo_file
HEADER            = args.header
THERMO_QUANTITIES = ["val1", "val2", "val3"]
END_FLAG          = args.end_flag
N_HIST_BINS       = args.nbin
N_SKIP_LINES      = args.skip

#sys.exit()

#
# Set the parameters
#

#FILENAME          = "o.sample_log.0"
#HEADER            = "Step val1 val2 val3"
#THERMO_QUANTITIES = ["val1", "val2", "val3"]
#END_FLAG          = "END FLAG"
#N_HIST_BINS       = 100
#N_SKIP_LINES      = 0


THERMO_DATA = {}
for thermo in THERMO_QUANTITIES:
   THERMO_DATA[thermo] = []

log_data = {}

for ifile in range(1):

   FILENAME_FULL = FILENAME

   if not os.path.exists(FILENAME_FULL):
       break
   #
   foo = open(FILENAME_FULL, 'r')

   #Skip the first lines untill reaching the header
   iline = 0
   read_flag = False
   while True:
       line = foo.readline()

       if line == '':
          break
       if END_FLAG in line:
          read_flag = False

       if read_flag:
          iline += 1
          if iline > N_SKIP_LINES:
              thermo_line = line.split()
              for ithermo in THERMO_QUANTITIES:
                  icol = col_of_thermo[ithermo]
                  THERMO_DATA[ithermo].append(float(thermo_line[icol]))

       if not read_flag and HEADER in line:
          read_flag = True
          header = line.split()
          col_of_thermo = {}

          icol = 0
          for ithermo in header:
              for jthermo in THERMO_QUANTITIES:
                  if ithermo == jthermo:
                      col_of_thermo[ithermo] = icol
              icol += 1
foo.close()

print("%-15s " % ("thermo"), end='')
for ithermo in THERMO_DATA:
   print("%-15s " % ithermo, end='')
print()

print("%-15s " % ("mean"), end='')
for data in THERMO_DATA.values():
   print("%-15f " % np.average(data), end='')
print()

print("%-15s " % ("std"), end='')
for data in THERMO_DATA.values():
   print("%-15f " % np.std(data, ddof=1), end='')
print()

print("%-15s " % ("err"), end='')
for data in THERMO_DATA.values():
   print("%-15f " % (np.std(data, ddof=1)/np.sqrt(len(data))), end='')
print()

# export histograms

for ithermo in THERMO_DATA:
   data    = THERMO_DATA[ithermo]
   min_val = min(data)
   max_val = max(data)
   nbin    = N_HIST_BINS
   step    = (max_val - min_val) / float(N_HIST_BINS)
   bins    = [min_val + ii*step for ii in range(nbin + 1)]

   foo = open("o.hist."+ithermo, 'w')

   hist, bin_edge = np.histogram(data, bins=bins, density=False)

   for ii in range(nbin):
      foo.write("%-15d %-15f %-15f\n" % (ii, bin_edge[ii], hist[ii]))
   foo.close()
