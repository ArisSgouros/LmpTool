import json
import sys

json_in = sys.argv[1]
with open(json_in) as f:
    d = json.load(f)

info_attr_1, info_val_1 = sys.argv[2], sys.argv[3]
info_attr_2, info_val_2 = sys.argv[4], sys.argv[5]
d['info'][info_attr_1] = info_val_1
d['info'][info_attr_2] = info_val_2
print(json.dumps(d, indent=2))
