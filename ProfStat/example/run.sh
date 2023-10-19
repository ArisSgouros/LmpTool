#!/bin/bash
echo -e 'Test: process file profile_sample\n'
python ../prof_stat_v4.py profile_sample

echo -e 'Test: process file profile_sample (values are multiplied w/ weights)\n'
python ../prof_stat_v4.py profile_sample -sum=1

echo -e 'Test: process file profile_sample (export mean only)\n'
python ../prof_stat_v4.py profile_sample -sum=0 -std=0 -err=0 -density=0
