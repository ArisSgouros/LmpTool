#!/bin/bash
path_smart="../../smart_generator.py"
path_smart_files="smart_files"
path_opls="../../opls_deployer.py"
path_ff="../../force_field/ff_oplsaa_foyer_2020_mod_SPE.xml"
python $path_smart $path_smart_files 200 SPCE 1 SPE 200 SPCE 1 SPE 200 SPCE 1 SPE 200 SPCE 1 SPE > o.log_smart
python $path_opls S1_4-W800_dreiding.data SMART_200-SPCE_1-SPE_200-SPCE_1-SPE_200-SPCE_1-SPE_200-SPCE_1-SPE.dat $path_ff S1_4-W800_opls.data > o.log_apply_opls
