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
parser.add_argument('-read_every', type=int, default=1, help='Set the reading frequency')
parser.add_argument('-density', type=int, default=1, help='Calculate the density of each bin')
parser.add_argument('-mean', type=int, default=1, help='Calculate the mean')
parser.add_argument('-std', type=int, default=1, help='Calculate the std')
parser.add_argument('-err', type=int, default=1, help='Calculate the standard error')


if __name__ == "__main__":
   args = parser.parse_args()
   filename      = args.filename
   read_every    = args.read_every
   print_density = bool(args.density)
   print_mean    = bool(args.mean)
   print_std     = bool(args.std)
   print_err     = bool(args.err)
   MultWeight    = bool(args.sum)

   print( 'File name     : ', filename)
   print( 'Read every    : ', read_every)
   print( 'Print density : ', print_density)
   print( 'Print mean    : ', print_mean)
   print( 'Print std     : ', print_std)
   print( 'Print err     : ', print_err)
   print( 'Sum           : ', print_err)

   # Read the header and get the number of Columns and number of Bins
   data_file = open(filename, 'r')
   data_file.readline()
   data_file.readline()

   # Read the number of Columns
   header = data_file.readline().split()
   n_col = len(header) - 4

   # Read the number of Bins
   n_bin = int(data_file.readline().split()[1])

   # Read the coordinate of each bin
   coord = [0.0 for ii in range(n_bin)]
   for ii in range(len(coord)):
      coord[ii] = float(data_file.readline().split()[1])

   print( 'number of columns   : ', n_col)
   print( 'number of Bins      : ', n_bin)
   print( 'Bin width           : ',str(coord[1]-coord[0]))

   data_file.close()

   prof_ave = [[0.0 for x in range(n_bin)] for y in range(n_col)]
   prof_var = [[0.0 for x in range(n_bin)] for y in range(n_col)]
   prof_dens = [0.0 for x in range(n_bin)]
   prof_weight = [0 for x in range(n_bin)]

   n_frame=0

   with open(sys.argv[1], 'r') as data_file:
      # Skip the first three lines
      data_file.readline()
      data_file.readline()
      data_file.readline()

      while True:
         # Skip frames
         if (n_frame % read_every != 0):
            for mm in range(n_bin+1):
               data_file.readline(),

         # Check if we reached end of file
         line = data_file.readline()
         if line == '': break

         n_frame += 1
         for mm in range(n_bin):
            line = data_file.readline().split()
            n_count = int(line[2])

            # Append the instances per bin to the density profile
            prof_dens[mm] += n_count

            if MultWeight:
               prof_weight[mm] += 1
               for cc in range(n_col):
                  x                 = n_count*float(line[cc+3])
                  delta             = x - prof_ave[cc][mm]
                  prof_ave[cc][mm] += delta/prof_weight[mm]
                  delta2            = x - prof_ave[cc][mm]
                  prof_var[cc][mm] += delta*delta2

            if not MultWeight:
               for ii in range(n_count):
                  prof_weight[mm] += 1
                  for cc in range(n_col):
                     x                 = float(line[cc+3])
                     delta             = x - prof_ave[cc][mm]
                     prof_ave[cc][mm] += delta/prof_weight[mm]
                     delta2            = x - prof_ave[cc][mm]
                     prof_var[cc][mm] += delta*delta2


   print("%-16s " %("coord"), end='')
   print( "%-16s " %("weight"), end='')
   if print_density:
      print( "%-16s " %("density"), end='')
   if print_mean:
      for jj in range(n_col):
         print( "%-16s " % ("ave_"+header[jj+4]), end='')
   if print_std and n_frame > 1:
      for jj in range(n_col):
         print( "%-16s " % ("std_"+header[jj+4]), end='')
   if print_err and n_frame > 1:
      for jj in range(n_col):
         print( "%-16s " % ("err_"+header[jj+4]), end='')
   print()

   for ii in range(n_bin):
      print( "%-16f " %(coord[ii]), end='')
      print( "%-16f " %(float(prof_weight[ii])), end='')
      if print_density:
         print( "%-16f " %(prof_dens[ii]/n_frame), end='')
      # Print the Averages
      if print_mean:
         for jj in range(n_col):
            print( "%-16f " % (prof_ave[jj][ii]), end='')
      # Print the STD
      if print_std and n_frame > 1:
         for jj in range(n_col):
            print( "%-16f " % (m.sqrt(prof_var[jj][ii]/(prof_weight[ii]-1))), end='')
      # Print the Error
      if print_err and n_frame > 1:
         for jj in range(n_col):
            print( "%-16f " %(m.sqrt(prof_var[jj][ii]/(prof_weight[ii]-1)/prof_weight[ii])), end='')

      print()
print()
