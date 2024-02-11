# OplsDeployer
Deploy OPLS force field

# Author
- Dr. Aristotelis P. Sgouros (arissgouros@gmail.com)

# Description
This script facilitates the deployment of the OPLS force field.
To utilize it, the user needs to follow these steps:

   Provide a LAMMPS data file
   Generate a SMART file to specify atom types.
   Provide a force field file.

While the script has been validated independently, it is recommended to leverage the MoSDeF framework [https://mosdef.org/] for its well-established testing and increased stability.


# Organization
The folder includes the following files and directories:
 - README             -> current file
 - LICENSE            -> MIT LICENSE
 - examples/          -> directory containing an indicative example
 - opls_deployer.py   -> python script deploying the OPLS force field
 - smart_generator.py -> python script generating smart atom types
