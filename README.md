# LmpTool
A suite of processing scripts for Lammps files

# Author
- Dr. Aristotelis P. Sgouros (arissgouros@gmail.com)

# Organization
The folder includes the following files and directories:
 - README               -> current file
 - LICENSE              -> MIT LICENSE
 - CountDumpFrame       -> Count frames of .lammpstrj and .xyz files
 - DatafileMerge        -> Merge Lammps data files
 - DatafileReplCoeff    -> Replace the force-field coefficients
 - DatafileRmvDuplCoeff -> Remove types with identical coefficients
 - DataToDump           -> Convert Lammps data to .lammpstrj or .xyz files
 - DumpConstAtom        -> Convert traj w/ nonconst atoms to vdf friendly format
 - DumpDecimator        -> Reduce the frame frequency of lammps dump files
 - DumpSortCol          -> Sort integer columns of dump files
 - DumpTypeStrip        -> Remove atoms with specific types from dump files
 - ProfStat             -> Export statistics of Lammps profiles
 - ThermoMerge          -> Merge multiple thermo files
 - ThermoStat           -> Export statistics of Lammps thermo outputs
 - MsdMol               -> Calculate the mean square displacement of molecules
 - Dev                  -> Miscellaneous code under development
