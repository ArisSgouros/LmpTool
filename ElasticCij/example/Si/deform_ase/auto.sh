#!/bin/bash
up="1.0e-6"
python deform.py --epsilon $up
python ../../../src/stress2cij.py --epsilon $up
mv o.stiffness.json o.stiffness_this_$up.json
