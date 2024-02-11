import sys

class c_group:
   def __init__(self, Id, nmol, segs):
      self.nmol = nmol
      self.segs = segs
      self.segs_per_mol = len(segs)
      self.Id = Id

      if self.segs_per_mol == 1:
         self.type = "SINGLE"
      else:
         self.type = "COPOL"

class c_atom:
   def __init__(self):
      self.Id_local = -1
      self.head_rmv = False
      self.tail_rmv = False
      self.head_change = False
      self.tail_change = False
      self.opls_name = "noname"
      self.opls_class = "noclass"
      self.opls_element = "noelement"
      self.comment = "nocoment"

class c_seg:
   def __init__(self):
      self.first = False
      self.last = False
      self.my_name = "noname"
      self.next_name = "void"
      self.prev_name = "void"

def read_SMART(seg_filename):

   try:
      with open(seg_filename) as goo:
         lines = goo.readlines()
   except:
      print("ERROR: FILE " + seg_filename + " NOT FOUND!")
      quit()

   atom_list = {}
   for line in lines:
      linesplit = line.split()
      iat = c_atom()
      iat.Id_local = int(linesplit[0][:-1])
      iat.opls_element = linesplit[1]
      iat.opls_name = linesplit[3]
      iat.opls_class = linesplit[4]

      try:
          iat.comment = line.split("<===")[1][:-2]
      except:
          iat.comment = " nocomment"

      if "HEAD_REMOVE" in line:
         iat.head_rmv = True
      if "TAIL_REMOVE" in line:
         iat.tail_rmv = True
      if "HEAD_CHANGE" in line:
         iat.head_change = True
      if "TAIL_CHANGE" in line:
         iat.tail_change = True

      atom_list[iat.Id_local] = iat

   return atom_list


# This function checks whether the group has been read by
# checking that the next arg is an integer (start of a new
# group), or that iarg+1 does not exist (end of last group)
def group_end (a, i):
   try:
      a[i+1]
      try:
         int(a[i+1])
         return True
      except:
         return False
   except:
      return True

#
# This section parses the command line arguments
#
if len(sys.argv) == 1:
   print("The format of the command arguments is the following:\n")
   print("python lmp_gen_SMART_MOL.py PATH_OF_SMART_FILES mol_A_number SEG_A1 SEG_A2 .. SEG_An mol_B_number SEG_B1 SEG_B2 .. SEG_Bm .. \n")
   print("Example commands:\n")
   print("SYSTEM with 1 SPE molecule")
   print("python lmp_gen_SMART_MOL.py ../SMART_FILES 1 SPE\n")
   print("SYSTEM with 1 SPE-VPA copolymer")
   print("python lmp_gen_SMART_MOL.py ../SMART_FILES 1 SPE VPA\n")
   print("SYSTEM with 10 SPCE water molecules")
   print("python lmp_gen_SMART_MOL.py ../SMART_FILES 10 SPCE\n")
   print("SYSTEM with 300 SPCE water molecules and 2 SPE-VPA copolymers")
   print("python lmp_gen_SMART_MOL.py ../SMART_FILES 10 SPCE 2 SPE VPA\n")
   print("exiting...")
   exit()


SMART_FILES_PATH = sys.argv[1] + "/"


sys.argv.pop(0)
sys.argv.pop(0)
args = sys.argv

groups = []
segs_list_aux = []
seg = c_seg()

print("\nParsing command line arguments and building connectivity..\n")

groupId = 1
for iarg in range(len(args)):
   try:
      nmol = int(args[iarg])
   except:
      segs_list_aux.append(args[iarg])

   if group_end (args, iarg):
      # 1. generate seg objects
      # 2. assess their order in group
      # 3. place them in segs list
      segs = []
      segId = 1
      nseg = len(segs_list_aux)
      for iseg in range(nseg):
         seg = c_seg()

         seg.my_name = segs_list_aux[iseg]
         seg.Id = segId

         # Assess whether these elements are at the ends of the group
         if iseg == 0:
            seg.first = True
         if iseg is nseg-1:
            seg.last = True

         # Find the previous and next segs
         if iseg < nseg-1:
            seg.next_name = segs_list_aux[iseg+1]
         if iseg > 0:
            seg.prev_name = segs_list_aux[iseg-1]
         segs.append(seg)

         segId += 1

      groups.append(c_group(groupId, nmol, segs))

      segs_list_aux = []
      groupId += 1

# Print output characteristics and generate the output filename

print("%-15s %-15s %-15s %-15s %-15s" % ("Group", "nmol", "type", "segs/mol", "mols"))

output_SMART_NAME = "SMART"
for ig in groups:

   output_SMART_NAME += "_" + str(ig.nmol) + "-"
   print( "%-15d %-15d %-15s %-15d" % (ig.Id, ig.nmol, ig.type, ig.segs_per_mol)),

   for seg in ig.segs:
      output_SMART_NAME += seg.my_name
      #print( "-15s" % (seg.my_name)),
      print( "%-15s" % seg.my_name),
   print("")

output_SMART_NAME += ".dat"

# Print the connectivity
for ig in groups:
   print( "\nConnectivity of Group %d.." % (ig.Id))
   print( "%-15s %-15s %-15s %-15s" % ("seg_ID", "prev_seg", "cur_seg", "next_seg"))
   for seg in ig.segs:
      print( "%-15d %-15s %-15s %-15s" % (seg.Id, seg.prev_name, seg.my_name, seg.next_name))

print("\nOutput file name: " + output_SMART_NAME + "\n")

foo = open(output_SMART_NAME, 'w')

# Start doing the thingy..
Id_global = 1
for ig in groups:
   print("Generating SMART file from Group %d" % (ig.Id))
   full_mol = []
   for seg in ig.segs:

      # Read the patch-with-prev SMART file if not first
      if not seg.first:
         #print "stitch with prev:"
         patch_filename = SMART_FILES_PATH + "SMART_PATCH_" + seg.my_name + "-" + seg.prev_name + ".dat"
         prev_patch_atom_list = read_SMART(patch_filename)
      else:
         prev_patch_atom_list = {}

      # Read the patch-with-prev SMART file if not last
      if not seg.last:
         #print "stitch with next:"
         patch_filename = SMART_FILES_PATH + "SMART_PATCH_" + seg.my_name + "-" + seg.next_name + ".dat"
         next_patch_atom_list = read_SMART(patch_filename)
      else:
         next_patch_atom_list = {}

      # Read the SMART file of the present seg
      seg_filename = SMART_FILES_PATH + "SMART_" + seg.my_name + ".dat"
      seg_atom_list = read_SMART(seg_filename)

      # Replace the OPLS classes if needed to patch the segs properly
      for key in seg_atom_list:
         iat = seg_atom_list[key]

         # remove the head/tail atom in case this is not the first/last seg
         if iat.head_rmv and not seg.first:
            continue
         if iat.tail_rmv and not seg.last:
            continue

         # replace the head_change/tail_chain atom in case this is not the first/last seg
         if iat.tail_change and not seg.last:
            iat = next_patch_atom_list[iat.Id_local]
         if iat.head_change and not seg.first:
            iat = prev_patch_atom_list[iat.Id_local]

         full_mol.append(iat)

   for imol in range(ig.nmol):
      for ii in range(len(full_mol)):
         iat = full_mol[ii]
         foo.write('%6d: %-5s==> %-19s %-6s <=== %s\n' % (Id_global, iat.opls_element, iat.opls_name, iat.opls_class, iat.comment))
         Id_global += 1

print("\nSUCCESS!")
