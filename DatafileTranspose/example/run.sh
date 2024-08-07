#!/bin/bash
python ../datafile_transpose.py graphene x y z -atomtype=full > graphene_xyz
python ../datafile_transpose.py graphene y x z -atomtype=full > graphene_yxz
python ../datafile_transpose.py graphene x z y -atomtype=full > graphene_xzy
python ../datafile_transpose.py graphene z x y -atomtype=full > graphene_zxy
