import sys
import os
import math as m


if len(sys.argv) != 8:
   print()
   print( "This script needs at least two arguments to run")
   print()
   print( "The required format is the following:")
   print()
   print( "python lmp_prof_ave.py FILENAME NEVERY MULT DENS MEAN STD ERR")
   print()
   print( "EVERY (INT): Calculare averages every NEVERY frames")
   print( "MULT (BOOL): Before averaging multiply the values with")
   print( "             beads per bin. Thats the proper way to do")
   print( "             it in the number of beads per bin. Thats the" )
   print( "             proper way to do it in stress calculations")
   print( "DENS (BOOL): Compute the density profile?")
   print( "MEAN (BOOL): Compute the averages?")
   print( "STD  (BOOL): Compute the Standard deviations?")
   print( "ERR  (BOOL): Compute the error of the mean?")
   print()
   print( "Example:")
   print()
   print( "python lmp_prof_ave.py profile.stresses 1 1 1 1 1 1")
   print()
   print( "This does everything mentioned above")
   print( "If any of the aformentioned funcionalities are not desired")
   print( "set the boolean to zero for those entries")
   print()
   print( "exiting..")
   sys.exit()



READ_EVERY = int(sys.argv[2])
#multiply_with_beads_per_bin = int(sys.argv[3])
print_DENS = int(sys.argv[4])
print_MEAN = int(sys.argv[5])
print_STD = int(sys.argv[6])
print_ERR = int(sys.argv[7])

apply_weight = True

#print "start"
#print "re",READ_EVERY
#print multiply_with_beads_per_bin
#print print_DENS
#print print_MEAN
#print print_STD
#print print_ERR
#print "end"

# Parameter section
#READ_EVERY = 1
#multiply_with_beads_per_bin = True
#print_DENS = 1
#print_MEAN = True
#print_STD = True
#print_ERR = True

# Read the header and get the number of Columns and number of Bins

data_file = open(sys.argv[1], 'r')

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

#print n_frame, n_frame % READ_EVERY

       n_frame += 1

       for bb in range(nBin):

           current_line = data_file.readline().split()
           n_atom = int(current_line[2])
           print("n_atom: ", n_atom)

           # Append the atoms per bin to the density profile
           nDens[bb] += n_atom

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
