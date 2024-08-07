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
 - DatafileSwapType     -> Swap atom/bond/angle/dihedral/improper types
 - DatafileRmvDuplCoeff -> Remove types with identical coefficients
 - DatafileTranspose    -> Transpose dimensions of lammps datafiles
 - DataToDump           -> Convert Lammps data to .lammpstrj or .xyz files
 - DumpConstAtom        -> Convert traj w/ nonconst atoms to vdf friendly format
 - DumpDecimator        -> Reduce the frame frequency of lammps dump files
 - DumpSortCol          -> Sort integer columns of dump files
 - DumpTypeStrip        -> Remove atoms with specific types from dump files
 - ProfStat             -> Export statistics of Lammps profiles
 - SheetToCnt           -> Fold 2D sheets to CNTs
 - ThermoMerge          -> Merge multiple thermo files
 - ThermoStat           -> Export statistics of Lammps thermo outputs
 - MsdMol               -> Calculate the mean square displacement of molecules
 - OplsDeployer         -> Deploy OPLS force field
 - Dev                  -> Miscellaneous code under development
