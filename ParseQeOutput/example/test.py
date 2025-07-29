import os
import sys

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.insert(0, parent_dir)

from parse_qe_output import write_xyz_file, parse_qe_output

for [case_, file_qe] in [["01", "i.MoSSe_3x3_01_vcrelax.out"],
                         ["02", "i.MoSSe_3x3_02_relax.out"],
                         ["03", "i.MoSSe_3x3_03_scf.out"],
                         ["04", "i.MoSSe_3x3_04_scf_no_force_no_stress.out"]]:

   print()
   print('####################################')
   print('Parsing file {:s}'.format(file_qe))
   frame_def, status = parse_qe_output(file_qe, xyzfmt = 'ase', debug = False)

   print()
   print('Display parsed properties:')
   # Check whether the properties were parsed
   for key in frame_def['cell_properties']:
      print(key, frame_def['cell_properties'][key])

   if ('forces' in frame_def['atoms'][0] or 'force' in frame_def['atoms'][0]):
      print('found forces')

   if 'virial' in frame_def['cell_properties']:
      print('found virial')

   # set the default weight
   frame_def['cell_properties']['weight'] = 1.0
   frame_def['cell_properties']['config_type'] = case_

   file_xyz = "o.{:s}_ase.xyz".format(case_)
   write_xyz_file(file_xyz, [frame_def], status='w', xyzfmt = 'ase')
   file_xyz = "o.{:s}_nep.xyz".format(case_)
   write_xyz_file(file_xyz, [frame_def], status='w', xyzfmt = 'nep')
