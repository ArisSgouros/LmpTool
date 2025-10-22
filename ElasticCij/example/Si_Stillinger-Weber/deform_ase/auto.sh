#!/bin/bash
up="1e-6"
python deform.py --epsilon $up
python emin.py
python ../../../src/stress2cij.py --epsilon $up
mv o.stiffness.json o.stiffness_this_$up.json
