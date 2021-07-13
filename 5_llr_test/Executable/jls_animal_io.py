#!/usr/bin/env python
"""
Provides I/O for animal trials.
"""

import pandas as pd

def _cols():
    return [ 'Bib', 'First Author', 'Year', 'PMID', 'Source', 
             'Virus', 'Macaque', 'Route', 'Weeks' ]

# Returns a dictionary with contents of 4_nested_test/Data/bib.csv.    
def to_uid2dict( bib_fn ):
    uid2dict = {}
    df = pd.read_csv( bib_fn )
    k2v = df.to_dict( 'list' )
    count = len( k2v.get('PMID') )
    for i in range( count ):
        key = str(k2v.get('PMID')[i]) + '_' + str(k2v.get('Source')[i]) # uid
        value = {}
        end_note = str(k2v.get('Bib')[i])
        value[_cols()[0]] = end_note
        s = end_note.split(',')
        value[_cols()[1]] = s[0][1:]          # First Author
        value[_cols()[2]] = int(s[1][1:].split()[0])                # Year
        for j in range(3, 9):
            value[_cols()[j]] = k2v.get(_cols()[j])[i] # PMID
        uid2dict[ key ] = value
        print( key,':',value, flush=True )
    return uid2dict

def main(): 
    pass
    
if __name__ == "__main__":
    main()
