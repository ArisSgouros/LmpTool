#!/bin/bash
ase2nep="../../ase2nep.py"
nep2ase="../../nep2ase.py"
python $ase2nep i.ase.xyz > o.ase_nep.xyz
python $nep2ase i.ase.xyz > o.ase_ase.xyz
python $ase2nep i.nep.xyz > o.nep_nep.xyz
python $nep2ase i.nep.xyz > o.nep_ase.xyz
