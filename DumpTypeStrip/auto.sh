DATAFILE="1_POS.dat"
DUMP_FILE="12_prod.lammpstrj"
NFRAME=750
EVFRAME=1

W_TYPE=6

source activate python3.7
#python lmp_MD2CG.py $DATAFILE $DUMP_FILE $NFRAME $EVFRAME
conda deactivate

source activate python2
python lmp_angle_dist_v2.py CG.dat [1] $NFRAME CG.lammpstrj
python lmp_angle_dist_v2.py CG.dat [2] $NFRAME CG.lammpstrj
python lmp_angle_dist_v2.py CG.dat [3] $NFRAME CG.lammpstrj
python lmp_angle_dist_v2.py CG.dat [4] $NFRAME CG.lammpstrj
#python lmp_type_strip.py CG.lammpstrj 3 "[$W_TYPE]" > CG.strip_$W_TYPE.lammpstrj
python lmp_bond_dist_v2.py CG.dat [1] $NFRAME CG.lammpstrj
python lmp_bond_dist_v2.py CG.dat [2] $NFRAME CG.lammpstrj
python lmp_bond_dist_v2.py CG.dat [3] $NFRAME CG.lammpstrj
python lmp_bond_dist_v2.py CG.dat [4] $NFRAME CG.lammpstrj
python lmp_bond_dist_v2.py CG.dat [5] $NFRAME CG.lammpstrj
python lmp_bond_dist_v2.py CG.dat [6] $NFRAME CG.lammpstrj
conda deactivate
