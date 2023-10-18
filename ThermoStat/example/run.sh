#!/bin/bash
python ../thermo_hist.py in.sample_log 'Step val1 val2 val3' 'END FLAG' -nbin 20 -skip 0 > o.log
