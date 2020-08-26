#!/usr/bin/env python
"""
Compile a HTML report from a variant .maf file
"""
import os
import sys
import csv
from jinja2 import FileSystemLoader, Environment, select_autoescape

template_dir = os.path.join(os.path.dirname(__file__), "templates")
loader = FileSystemLoader(template_dir)
environment = Environment(
    loader = loader,
    autoescape = select_autoescape(['html'])
    )
template = environment.get_template('report.html')

input_file = "analyst_file.new.maf"
with open(input_file) as fin:
    reader = csv.DictReader(fin, delimiter = '\t')
    row_list = [ row for row in reader ]

parsed = template.render(row_list = row_list)

print(parsed, file = sys.stdout)
