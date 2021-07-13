#!/usr/bin/env python

from os import system

log = f'to_summary.log'

# Defaults follow each option as a comment after '#'. 
I = ' -i ../Data/' # -i ../Data/ # input directory
O = ' -o ../Output/' # -o ../Output/ # output directory   
C = ' -c 0.95' # -c 0.95 # confidence limit
B = ' -b ../Data/Bib/bib.csv' # -b ../Data/Bib/bib.csv # input file with bibliographic metadata  
S = ' -s' # analyze controls and treatment arms # -s analyzes control arms only 

system( f'python to_summary.py {I} {O} {C} {B} {S} > {log}' )
