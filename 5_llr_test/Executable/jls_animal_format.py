#!/usr/bin/env python
"""
Provides subroutines for general data manipulation in animal trials.
"""

import math

#from scipy.stats import norm
import pandas as pd

def to_df( fn ):
    return pd.read_csv( fn, usecols = ['Challenge', 'Uninfected'] )

def exact_stats(): # sufficient exact statistics
    return [ 
        "# of maximum challenges for uninfected animals", # None, if all animals are infected.
        "# animals", 
        "# uninfected animals", 
        "# non-infecting challenge doses total given to infected animals" 
    ]
    
# Returns a dictionary with keys=names of exact stats 2 values.
def to_exact_stats( canonical_npts ): # canonical_npts = df.values.tolist()
    n = canonical_npts[0][1]  # animals
    u = canonical_npts[-1][1] # uninfected animals
    if u == 0:
        m = None # no uinfected animals
    else:
        m = canonical_npts[-1][0] # maximum # of challenges for uninfected animals
    v = 0 # total # non-infecting challenge doses given to infected animals
    for npt in canonical_npts:
        if npt is canonical_npts[0]:
            assert( npt[0] == 0 )
            npt0 = npt
            continue
        v += (npt[0] - npt0[0]) * (npt0[1] - u)
        npt0 = npt
    v -= n - u
    return dict( zip( exact_stats(), [m,n,u,v] ) )

def _to_exact_stats_test():
    npts2stats = {
        ( (0,8), (4,5), (5,4), (6,3), (10,3) ) : [10,8,3,18], 
        ( (0,8), (1,4), (2,2), (4,1), (10,1) ) : [10,8,1,5],
        ( (0,8), (1,4), (2,2), (4,0) )         : [None,8,0,8]
    }
    keys = exact_stats()
    for (k, v) in npts2stats.items():
        assert( list( to_exact_stats( k ).keys() ) == keys )
        assert( list( to_exact_stats( k ).values() ) == v )

# Returns a list whose elements are n[t+1] = #(uninfected survivors after challenge t) for t=0...t[max].
def to_survivors( canonical_npts, m=None ): # canonical_npts = df.values.tolist(), m = # maximum challenges
    survivors = []
    for npt in canonical_npts:
        if npt is canonical_npts[0]:
            assert( npt[0] == 0 )
            npt0 = npt
            continue
        for c in range( npt0[0], npt[0] ):
            survivors.append( npt0[1] )
        npt0 = npt
    survivors.append( canonical_npts[-1][1] )
    return survivors

def _to_survivors_test():
    npts2survivors = {
        ( (0,5), (2,3), (4,2), (5,2) ) : [5,5,3,3,2,2], # repeated final value
        ( (0,5), (2,3), (4,2), (5,1) ) : [5,5,3,3,2,1],
        ( (0,5), (2,3), (4,2), (5,0) ) : [5,5,3,3,2,0],
        ( (0,5), (2,3), (4,0), (5,0) ) : [5,5,3,3,0,0]  # robust against repeated 0
    }
    for (k, v) in npts2survivors.items():
        assert( to_survivors( k ) == v )

# Returns a list whose elements are (t[k], n[t[k]] = #(uninfected survivors after challenge t[k]),
#   where t[k] = min{ t > t[k - 1] : n[t] < n[t - 1] } .
def to_canonical_npts( survivors ): # the inverse function of to_survivors( canonical_npts )
    canonical_npts = [ [ 0, survivors[0] ] ]
    tks = [ i != j for i,j in list(zip( survivors[1:], survivors)) ]
    canonical_npts.extend( [ [ i + 1, survivors[i + 1] ] for i in range( len(survivors) - 2 ) if tks[i] ] )
    canonical_npts.append( [ len(survivors) - 1, survivors[-1] ] )
    return canonical_npts

def _to_canonical_npts_test():
    npts2survivors = {
        ( (0,5), (2,3), (4,2), (5,2) ) : [5,5,3,3,2,2], # repeated final value
        ( (0,5), (2,3), (4,2), (5,1) ) : [5,5,3,3,2,1],
        ( (0,5), (2,3), (4,2), (5,0) ) : [5,5,3,3,2,0],
        ( (0,5), (2,3), (4,0), (5,0) ) : [5,5,3,3,0,0], # robust against repeated 0
        ( (0,8), (4,5), (5,4), (6,3), (10,3) ) : [8,8,8,8,5,4,3,3,3,3,3], 
        ( (0,8), (1,4), (2,2), (4,1), (10,1) ) : [8,4,2,2,1,1,1,1,1,1,1],
        ( (0,8), (1,4), (2,2), (4,0) )         : [8,4,2,2,0]
    }
    for (k, v) in npts2survivors.items():
        assert( to_canonical_npts( v ) == list(map(list,k)) )

# Returns a list counting the challenges before each animal is infected, or m+1 if animal is uninfected
def to_challenges( survivors ): 
    infections = [ i - j for i,j in zip( survivors, survivors[1:] ) ]
    challenges = []
    for i in range( len( infections ) ):
        challenges.extend( [ i + 1 ] * infections[i]  )
    challenges.extend( [ len( infections ) + 1 ] * survivors[-1] )
    return challenges

def _to_challenges_test():
    survivors2challenges = {
        (5,5,3,3,2,2) : [2,2,4,6,6], # repeated final value
        (5,5,3,3,2,0) : [2,2,4,5,5], 
        (5,5,3,3,0,0) : [2,2,4,4,4], 
        (8,8,8,8,5,4,3,3,3,3,3) : [4,4,4,5,6,11,11,11], 
        (8,4,2,2,1,1,1,1,1,1,1) : [1,1,1,1,2,2,4,11], 
        (8,4,2,2,0) : [1,1,1,1,2,2,4,4] 
    }
    for (k, v) in survivors2challenges.items():
        assert( to_challenges( k ) == v )

# Pads a survivors with final 0 out to a list of indicated length.
#   Does nothing if final element != 0 or length is None or len( survivors ) == length.
#   Raises exception if length < len( survivors ).
def pad_survivors( survivors, length=None ): 
    if length is None or len( survivors ) == length:
        return survivors
    if length < len( survivors ):
        raise Exception( 'survivors {} is already longer than length {}.', survivors, length )
    if survivors[-1] != 0:
        pad = float('NaN')
    else:
        pad = 0
    survivors = list( survivors )
    for i in range( len( survivors ), length ):
        survivors.append( pad )
    return survivors

def _pad_survivors_test():
    NaN = float('NaN')
    npts2survivors = {
        (5,5,3,3,2,0) : [5,5,3,3,2,0,0,0],
        (5,5,3,3,0,0) : [5,5,3,3,0,0,0,0],  # robust against repeated 0
        (5,5,3,3,2,2) : [5,5,3,3,2,2,NaN,NaN]  # robust against repeated 0
        }
    for (k, v) in npts2survivors.items():
        assert( pad_survivors( k, None ) == k )
    for (k, v) in npts2survivors.items():
        assert( pad_survivors( k, 6 ) == k )
    for (k, v) in npts2survivors.items():
        v0 = pad_survivors( k, 8 )
        assert( all( [ i == j or (math.isnan( i ) and math.isnan( j )) for i,j in zip(v,v0) ] ) )

# Removes NaNs and ensures the result is a list of integers.
def rstrip_nan( survivors ): 
    survivors_no_nans = [ int( survivor ) for survivor in survivors if not math.isnan( survivor )]
    if not all( [ i == j for i,j in zip(survivors_no_nans, survivors) ] ):
        raise Exception('The survivor list has internal NaNs.')
    return survivors_no_nans

def _rstrip_nan_test():
    NaN = float('NaN')
    npts2survivors = {
        (5,5,3,3,2,0) : [5,5,3,3,2,0],
        (5,5,3,3,0,0,NaN) : [5,5,3,3,0,0],
        (5,5,3,3,2,2,NaN,NaN) : [5,5,3,3,2,2]  # robust against repeated 0
        }
    for (k, v) in npts2survivors.items():
        assert( rstrip_nan( k ) == v )

# Returns True if survivors contains decreasing non-negative integers, not all the same.
def is_survivors( survivors ): # # uninfected animals after Challenge i.
    if not all( isinstance( n, int ) and n >= 0 for n in survivors ):
        return False
    if len( survivors ) <= 1:
        return False
    successive_pairs = list( zip( survivors, survivors[1:] ) )
    if not all( j <= i for i,j in successive_pairs ):
        return False
    return True

def _is_survivors_test():
    assert( is_survivors( [5,5,3,3,2,2] ) )
    assert( is_survivors( [5,5,3,3,2,1] ) )
    assert( is_survivors( [5,5,3,3,2,0] ) )
    assert( is_survivors( [5,5,3,3,0,0] ) )
    assert( is_survivors( [5,5] ) )
    assert( not is_survivors( [5] ) )
    assert( not is_survivors( [5,5,6,3,0,0] ) )

# Returns True if infection occurs.
#   Assumes is_survivors( survivors ).
def is_infection( survivors ): # # uninfected animals after Challenge i.
    return survivors[-1] < survivors[0]

def _is_infection_test():
    assert( is_infection( [5,5,3,3,2,2] ) )
    assert( is_infection( [5,5,3,3,2,1] ) )
    assert( is_infection( [5,5,3,3,2,0] ) )
    assert( is_infection( [5,5,3,3,0,0] ) )
    assert( not is_infection( [5] ) )
    assert( not is_infection( [5,5] ) )
    assert( not is_infection( [5,5,5,5,5,5] ) )

# Converts survivors n[t+1] to a list [n[t+1],n[t]-n[t+1]] of 2-lists [at risk, infected] at Challenge t, shortened by 1.
#   The columns are for zipping into tables for logrank_statistic.
def _to_columns( survivors ):
    columns = []
    npt0 = None
    for npt in survivors:
        if npt0 is None:
            npt0 = npt
            continue
        columns.append( [npt0, npt0 - npt] )
        npt0 = npt
    return columns

def _to_columns_test():
    npts2survivors = {
        (5,5,3,3,2,2):[[5,0],[5,2],[3,0],[3,1],[2,0]], # repeated final value
        (5,5,3,3,2,1):[[5,0],[5,2],[3,0],[3,1],[2,1]],
        (5,5,3,3,2,0):[[5,0],[5,2],[3,0],[3,1],[2,2]],
        (5,5,3,3,0,0):[[5,0],[5,2],[3,0],[3,3],[0,0]]  # robust against repeated 0
    }
    for (k, v) in npts2survivors.items():
        assert( _to_columns( k ) == v )

# Converts control and vaccine survivors, challenge by challenge, to a list of tables for logrank_statistic.
def to_tables( survivors_control, survivors_treated ):
    return list( zip( _to_columns( survivors_control ), _to_columns( survivors_treated ) ) )

def _to_tables_test():
    survivors2tables = {
        ((5,5,3,3,2,2),(5,5,3,3,2,1)) : [([5,0],[5,0]),([5,2],[5,2]),([3,0],[3,0]),([3,1],[3,1]),([2,0],[2,1])], # repeated final value
        ((5,5,3,3,2,0),(5,5,3,3,0,0)) : [([5,0],[5,0]),([5,2],[5,2]),([3,0],[3,0]),([3,1],[3,3]),([2,2],[0,0])],
        ((5,5,3,3,2,0),(5,5,3,3,0))   : [([5,0],[5,0]),([5,2],[5,2]),([3,0],[3,0]),([3,1],[3,3])]  # different lengths
    }
    for (k, v) in survivors2tables.items():
        assert( to_tables( k[0], k[1] ) == v )

# Converts survivors in a row of a Fisher 2x2 table.
#   The whole table is [ to_fisher_row( survivors_control ),to_fisher_row( survivors_treated ) ].
def to_fisher_row( survivors ):
    return [survivors[0] - survivors[-1], survivors[-1]]

def _to_fisher_row_test():
    survivors2row = {
        (5,5,3,3,2,2):[3,2], 
        (5,5,3,3,2,1):[4,1], 
        (5,5,3,3,2,0):[5,0],
        (5,5,3,3,0,0):[5,0],
        (5,5,3,2):[3,2],
        (6,5,3,3):[3,3]
    }
    for (k, v) in survivors2row.items():
        assert( to_fisher_row( k ) == v )

def main():
    _to_exact_stats_test()
    _to_survivors_test()
    _to_canonical_npts_test()
    _to_challenges_test()
    _is_survivors_test()
    _is_infection_test()
    _pad_survivors_test()
    _rstrip_nan_test()
    _to_columns_test()
    _to_tables_test()
    _to_fisher_row_test()
    
if __name__ == "__main__":
    main()
