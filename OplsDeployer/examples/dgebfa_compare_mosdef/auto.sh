#!/bin/bash
path_smart="../../smart_generator.py"
path_smart_files="smart_files"
path_opls="../../opls_deployer.py"
path_ff='../../force_field/ff_oplsaa_foyer_2023_mod.xml'
name='dgebfa'
path_rmv_dupl_coeff="../../../DatafileRmvDuplCoeff/rmv_dupl_coeff.py"

# apply opls force field method from MOSDEF
python mosdef/mosdef_opls.py $name $path_ff > o.log_mosdef

# remove duplicate coeffs
python $path_rmv_dupl_coeff "$name"_mosdef.data "$name"_mosdef_no_dupl.data  -unique_pair_coeff=0

# apply opls force field method from OplsDeployer
python $path_smart $path_smart_files 30 DGEBF_A > o.log_smart
python $path_opls "$name"_mosdef.data SMART_30-DGEBF_A.dat $path_ff "$name"_opls_depl.data > o.log_opls_deployer
