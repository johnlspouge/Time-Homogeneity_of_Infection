#!/usr/bin/env python
"""
Reads in properly named digitization CSVs from https://apps.automeris.io/wpd/ and writes a summary CSV.
"""

import argparse
from sys import exit
from os.path import exists, isfile
from os import listdir, mkdir
import math
import pandas as pd

from jls_animal_model import ConstantHazard, ArithmeticPriming, GeometricPriming, StepPriming, BetaFrailty, DeltaFrailty
from jls_animal_io import to_uid2dict

CSV = '.csv'
SEP = '_'

abbrs = 'AP,GP,SP1,SP2,SP3,BF,DF'.split(',')
abbr_two2model = {
   'AP': {'model':ArithmeticPriming, 'param':'eps'},
   'GP': {'model':GeometricPriming, 'param':'r'},
   'SP': {'model':StepPriming, 'param':'p2'},
   'BF': {'model':BetaFrailty, 'param':'p_var'},
   'DF': {'model':DeltaFrailty, 'param':'delta'} }

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    df_summary = to_df_summary( argument )    
    df_summary.to_csv( f'{argument.odir}summary.csv', index=False )

def to_df_summary( argument ):
    print( cols(), flush=True )
    df_summary = pd.DataFrame( columns=cols() )
    uid2dict = to_uid2dict( argument.bib )
    fns_all = sorted( listdir( argument.idir ) )
    fns_csv = [ fn for fn in fns_all if fn.endswith( CSV ) ]
    for fn in fns_csv:
        df_trial = pd.read_csv( f'{argument.idir}{fn}' )
        arms = to_arms( argument, df_trial )
        fbn = fn[:-4]
        if fbn not in uid2dict:
            raise Exception( fn )
        d = uid2dict.get( fbn )
        print( d.get('PMID'), d.get('Source'), arms, flush=True )
        for arm in arms:
            row = [ d.get(cols()[i]) for i in range(8) ]
            row.append( arm ) # [...arm] 
            survivors = [ int( survivor ) for survivor in df_trial[ arm ].tolist() if not math.isnan( survivor ) ]
            print( survivors, flush=True )
            row.append( survivors[0] )
            row.append( survivors[-1] )
            row.append( len( survivors) - 1 )
            # reduced model
            model = ConstantHazard( survivors )
            mle = model.mle()
            row.append( float(mle[0]) )
            row.append( model.ln_likelihood( mle ) )
            row.append( float(model.fisher_information( mle )[0][0]) )
            row.extend( model.lr_interval( argument.confidence ) )
            # full models
            for abbr in abbrs:
                l_step = None
                if len( abbr ) == 3:
                    l_step = int( abbr[2] )
                    abbr = abbr[:2]
                if len( abbr ) != 2:
                    raise Exception('The model abbreviation is neither 2 nor 3 characters')
                row = to_append( row, survivors, abbr_two2model[ abbr ]['model'], l_step )
            print( row, flush=True )
            df_summary.loc[len(df_summary.index)] = row
    return df_summary

# Returns row after appending results for Model( survivors ).
def to_append( row, survivors, Model, l_step=None ):
    if l_step is None:
        model = Model( survivors )
    else:
        model = Model( survivors, l_step )
    mle = model.mle()
    fisher_information = model.fisher_information( mle )
    row.append( model.llr_pvalue() )
    row.extend( mle )
    row.append( model.ln_likelihood( mle ) )
    row.append( fisher_information[0][0] )
    row.append( fisher_information[1][0] )
    row.append( fisher_information[1][1] )
    return row

def cols():
    cols = [
        'First Author', 'Year', 'PMID', 'Source', 'Virus', 'Macaque', 'Route', 'Weeks',
        'Arm', 'n', 'n[-1]', 'm', 'p_RM', 'l_RM', 'I_RM', 'p_lo_RM', 'p_hi_RM' ]
    for abbr in abbrs:
        cols.extend( [ 'pval_' + abbr, 
                       'p_' + abbr, abbr_two2model[ abbr[:2] ]['param'] + '_' + abbr, 
                       'l_' + abbr, 'I11_' + abbr, 'I12_' + abbr, 'I22_' + abbr ] )
    return cols

def to_arms( argument, df_trial ):
    assert 'Control' in df_trial.columns 
    if argument.sham:
        arms = ['Control']
    else:
        arms = list( df_trial ) # headers
        arms.remove('Challenge')
    return arms

# Check and fixes arguments if possible.    
def check( argument ): 
    if not exists( f'{argument.idir}' ):
        print( f'Error: a valid INPUT_DIRECTORY "{argument.idir}" is required.' )
        exit()
    if not argument.idir.endswith('/'):
        argument.idir += '/'
    if not exists( f'{argument.odir}' ):
        mkdir( f'{argument.odir}' )
    if not argument.odir.endswith('/'):
        argument.odir += '/'
    if not isfile( f'{argument.bib}' ):
        raise Exception( f'Error: a valid BIB "{argument.bib}" is required.' )
        
def getArguments():
    parser = argparse.ArgumentParser(description='Reads in digitization CSVs with formatted filename for each trial, and writes a summary CSV.\n')
    parser.add_argument("-i", "--idir", dest="idir", default="../Data/", # input directory
                        help="INPUT_DIRECTORY", metavar="INPUT_DIRECTORY")
    parser.add_argument("-o", "--odir", dest="odir", default="../Output/", # input directory
                        help="OUTPUT_DIRECTORY", metavar="OUTPUT_DIRECTORY")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # confidence interval for restricted model
                        help="CONFIDENCE", metavar="CONFIDENCE")
    parser.add_argument("-b", "--bib", dest="bib", default='../Data/Bib/bib.csv', # bibliography
                        help="BIB", metavar="BIB")
    parser.add_argument("-s", "--sham", dest="sham", default=False, action="store_true", # ? controls only ?
                        help="SHAM")
    return parser
    
if __name__ == "__main__":
    main()
