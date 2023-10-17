#!/bin/bash
echo -e "Testing thermo_merge:"
echo -e "  header string : 'Step KinEng aux1'"
echo -e "  output        : 'o.log_step_kineng_aux1'"
python ../thermo_merge.py o.sample_log "Step KinEng aux1" > o.log_step_kineng_aux1

echo -e "Testing thermo_merge:"
echo -e "  header string : 'Step KinEng aux1'"
echo -e "  hide last     : True"
echo -e "  output        : 'o.log_step_kineng_aux1_hide_last'"
python ../thermo_merge.py o.sample_log "Step KinEng aux1" -hide_last 1 > o.log_step_kineng_aux1_hide_last

echo -e "Testing thermo_merge:"
echo -e "  header string : 'Step KinEng aux1 aux2'"
echo -e "  output        : 'o.log_step_kineng_aux1_aux2'"
python ../thermo_merge.py o.sample_log "Step KinEng aux1 aux2" > o.log_step_kineng_aux1_aux2

echo -e "Testing thermo_merge:"
echo -e "  header string : 'Step KinEng aux1 aux2'"
echo -e "  hide last     : True"
echo -e "  output        : 'o.log_step_kineng_aux1_aux2_hide_last'"
python ../thermo_merge.py o.sample_log "Step KinEng aux1 aux2" -hide_last 1 > o.log_step_kineng_aux1_aux2_hide_last
