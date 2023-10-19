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
import argparse
import os
import math as m

parser = argparse.ArgumentParser(description='Statistics of Lammps profiles')
parser.add_argument('filename', type=str, help='Path of the profile file')
parser.add_argument('-sum', type=int, default=0, help='Multiply the values with the weights')
parser.add_argument('-every', type=int, default=1, help='Set the reading frequency')
parser.add_argument('-density', type=int, default=1, help='Calculate the density of each bin')
parser.add_argument('-mean', type=int, default=1, help='Calculate the mean')
parser.add_argument('-std', type=int, default=1, help='Calculate the std')
parser.add_argument('-err', type=int, default=1, help='Calculate the standard error')


#if __name__ == "__main__":
args = parser.parse_args()
FILENAME = args.filename
READ_EVERY = args.every
print_DENS = bool(args.density)
print_MEAN = bool(args.mean)
print_STD = bool(args.std)
print_ERR = bool(args.err)
SUM = bool(args.sum)

# Read the header and get the number of Columns and number of Bins

data_file = open(FILENAME, 'r')

data_file.readline()
data_file.readline()

# Read the number of Columns
header = data_file.readline().split()
nCol = len(header) - 4

# Read the number of Bins
nBin = int(data_file.readline().split()[1])

# Read the coordinate of each bin
coord = [0.0 for ii in range(nBin)]
for ii in range(len(coord)):
    coord[ii] = float(data_file.readline().split()[1])

print( 'number of columns: ' ,nCol)
print( 'number of Bins: '    ,nBin)
print( 'Bin width: '    ,str(coord[1]-coord[0]))

data_file.close()

gM1 = [[0.0 for x in range(nBin)] for y in range(nCol)]
gM2 = [[0.0 for x in range(nBin)] for y in range(nCol)]
nDens = [0.0 for x in range(nBin)]
bin_count = [0 for x in range(nBin)]

n_frame=0

with open(sys.argv[1], 'r') as data_file:

   # Skip the first three lines
   data_file.readline()
   data_file.readline()
   data_file.readline()

#   for kk in range(1):
   while True:

       # Skip frames
       if (n_frame % READ_EVERY != 0):
           for bb in range(nBin+1):
               data_file.readline(),

       # Check if we reached end of file
       current_line = data_file.readline()
       if current_line == '': break

       n_frame += 1

       for bb in range(nBin):

           current_line = data_file.readline().split()
           n_atom = int(current_line[2])

           # Append the atoms per bin to the density profile
           nDens[bb] += n_atom

           if SUM:
              bin_count[bb] += 1
              for cc in range(nCol):
                 x = n_atom*float(current_line[cc+3])
                 delta = x - gM1[cc][bb]
                 gM1[cc][bb] += delta/bin_count[bb]
                 delta2 = x - gM1[cc][bb]
                 gM2[cc][bb] += delta*delta2

           if not SUM:
              for ii in range(n_atom):
                 bin_count[bb] += 1
                 for cc in range(nCol):
                    #if multiply_with_beads_per_bin:
                    #    x = n_atom*float(current_line[cc+3])
                    #else:
                    x = float(current_line[cc+3])
                    delta = x - gM1[cc][bb]
                    gM1[cc][bb] += delta/bin_count[bb]
                    delta2 = x - gM1[cc][bb]
                    gM2[cc][bb] += delta*delta2


print("%-16s" %("coord"), end='')
if print_DENS:
    print( "%-16s" %("DENS"), end='')
if print_MEAN:
    for jj in range(nCol):
        print( "AVE_%-12s " % (header[jj+4]), end='')
if print_STD and n_frame > 1:
    for jj in range(nCol):
        print( "STD_%-12s " % (header[jj+4]), end='')
if print_ERR and n_frame > 1:
    for jj in range(nCol):
        print( "ERR_%-12s " % (header[jj+4]), end='')
print('\n', end=' ')

for ii in range(nBin):
    print( "%-16f " %(coord[ii]), end='')
    if print_DENS:
        print( "%-16f " %(nDens[ii]/n_frame), end='')
    # Print the Averages
    if print_MEAN:
        for jj in range(nCol):
            print( "%-16f " % (gM1[jj][ii]), end='')
    # Print the STD
    if print_STD and n_frame > 1:
        for jj in range(nCol):
            print( "%-16f " % (m.sqrt(gM2[jj][ii]/(bin_count[ii]-1))), end='')
    # Print the Error
    if print_ERR and n_frame > 1:
        for jj in range(nCol):
            print( "%-16f " %(m.sqrt(gM2[jj][ii]/(bin_count[ii]-1)/bin_count[ii])), end='')

    print( '\n', end=' ')
