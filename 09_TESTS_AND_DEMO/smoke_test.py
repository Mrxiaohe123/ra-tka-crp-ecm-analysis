#!/usr/bin/env python3
import csv, math
from pathlib import Path
p=Path(__file__).with_name('clinical_synthetic_example.csv')
rows=list(csv.DictReader(p.open()))
assert len(rows)==8 and all(r['patient'].startswith('SYNTHETIC_') for r in rows)
a=sum(math.log(float(r['test_pre'])) for r in rows)/8
b=sum(math.log(float(r['test_m6'])) for r in rows)/8
print({'synthetic_rows':8,'synthetic_log_m6_pre_difference':b-a})
