#!/usr/bin/env python
"""
Provides minus ln-likelihoods for the models in animal trials.
    jls_animal_nested_constant_hazard.py contains the regressive tests.
    The full models require survivors be set (as a global variable for minimizations).
    All routines presume jls_animal_format.is_survivors( survivors ).

Derived Classes of Abstract Base Class AnimalModel:
    ConstantHazard
Derived Classes of Base Class ConstantHazard: (Nested Models)
    ArithmeticPriming
    GeometricPriming
    StepPriming
    BetaFrailty
    DeltaFrailty    
"""

import warnings

import math
import numpy as np
import scipy.optimize as opt
from scipy.stats import chi2
import numdifftools as ndt

from abc import ABC, abstractmethod

from jls_animal_format import is_survivors, rstrip_nan, is_infection
 
class AnimalModel(ABC):
    def __init__(self, survivors ):
        survivors = rstrip_nan( survivors )
        if not is_survivors( survivors ):
            raise Exception('invalid list of surccessive survivor counts')
        self.survivors = survivors # uninfected animals after Challenge i.
        self.pts = AnimalModel._to_pts( self.survivors )
        super().__init__()
    # Returns model name.
    def name(self):
        return type(self).__name__
    # Returns probability of infection corresponding to i = pt[0] = t-1.
    @abstractmethod
    def p_infection(self, i, x ):
        pass
    # Returns False if x violates bounds.
    @staticmethod
    @abstractmethod
    def is_in_bounds(x):
        pass
    # Returns True if the null model is on the boundary of the model parameter space.
    # Returns None if the model is the null model.
    @staticmethod
    def is_null_on_boundary(x):
        pass
    # Returns mle for constant hazard of infection as a scalar.
    @staticmethod
    @abstractmethod
    def x0( survivors ):
        return (survivors[ 0 ] - survivors[ -1 ]) / sum( survivors[ :-1 ] )
    #  Depends on individual models to calculate p_infection, probability of infection.
    #  pt = [ t-1, ns[t], ds[t] ], where ns count the challenged animals; and ds, the deaths.
    def ln_likelihood(self, x ):
        if not is_infection(self.survivors ) and np.allclose( x, self.x0(self.survivors ) ):
            return 0.0
        if not self.is_in_bounds(x):
            return -math.inf
        ln_likelihood = 0.0
        for pt in self.pts:
            ln_likelihood += AnimalModel._add_ln_likelihood( pt, self.p_infection( pt[0], x ) )
        return ln_likelihood 
    #  -self.ln_likelihood( x ) for minimization in scipy.opt.
    def _minus_ln_likelihood(self, x ):
        return -self.ln_likelihood( x ) 
    # Returns the maximum likelihood estimator as an array, even in one dimension.
    def mle(self, method='Basinhopping' ):
        #print(method)
        x0 = self.x0( self.survivors )
        if not is_infection( self.survivors ):
            return x0
        with warnings.catch_warnings():
            warnings.filterwarnings( "ignore", category=RuntimeWarning )
            _EPS = 1.0e-06
            if method == 'Nelder-Mead':
                optimum = opt.minimize( self._minus_ln_likelihood, x0, method='Nelder-Mead', 
                             bounds=None, tol=None, callback=None, 
                             options={'xatol': _EPS, 'fatol': _EPS, 'maxiter': None, 'maxfev': None, 'disp': False, 'return_all': False, 'adaptive': True})
            elif method == 'Powell':
                optimum = opt.minimize( self._minus_ln_likelihood, x0, method='Powell', 
                             bounds=None, tol=None, callback=None, 
                             options={'xtol': _EPS, 'ftol': _EPS, 'maxiter': None, 'maxfev': None, 'disp': False, 'return_all': False})
                if len( x0 ) == 1: # Converts Powell optimum to list for consistency.
                    optimum.x = [optimum.get('x').tolist()]
                    #print(optimum.x)
            elif method == 'Basinhopping':
                optimum = opt.basinhopping( self._minus_ln_likelihood, x0 )
            else:
                raise Exception('unknown optimization method')
        return optimum.get('x')
    # Returns arrays of NaN if is_in_bounds(x) but on the boundary.
    def fisher_information(self, x): # usually the maximum likelihood estimator
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return ndt.Hessian( self._minus_ln_likelihood )(x)

##########################################################################################
# private routines
##########################################################################################

    # Returns pts[t] = [ t-1, ns[t], ds[t] ], where ns count the challenged animals; and ds, the deaths.
    @staticmethod
    def _to_pts(survivors):
        ns = survivors[:-1]
        ds =  [ i - j for i,j in list( zip( survivors, survivors[1:] ) ) if not math.isnan( j ) ]
        assert( len(ns) == len(ds) )
        return list( zip( range( len(ns) ), ns, ds ) ) # (t-1, ns[t], ds[t])
    # Returns the increments to ln_likelihood from Challenge t.
    #   Depends on individual models to calculate p_infection, probability of infection.
    #   pt = [ t-1, ns[t], ds[t] ], where ns count the challenged animals; and ds, the deaths.
    @staticmethod
    def _add_ln_likelihood( pt, p_infection ):
        p_infection = min( 1.0, max( 0.0, p_infection ) ) # sentinels
        if p_infection == 0.0 and pt[2] == 0:
            return 0.0
        elif p_infection == 1.0 and pt[2] == pt[1]:
            return 0.0
        elif p_infection == 0.0 or p_infection == 1.0: # impossibility
            return -math.inf
        ln_p = math.log( p_infection, math.e ) # ln probability of deaths
        ln_q = math.log( 1.0 - p_infection, math.e ) # ln probability of non-deaths
        return pt[2] * ln_p + ( pt[1] - pt[2] ) * ln_q

def _test_AnimalModel():
    survivors = [ 64, 32, 16, 8, 4, 2, 1 ]
    assert( AnimalModel.x0( survivors ) == 0.5 )
    survivors = [ 64, 16, 4, 1 ]
    assert( AnimalModel.x0( survivors ) == 0.75 )
    
##########################################################################################
# derived classes
##########################################################################################

_EPS = 0.003 # accuracy for numerical tests
_METHODS = [ 'Nelder-Mead', 'Powell', 'Basinhopping' ]

# R.R. Regoes et al. (2005) Preclinical assessment of HIV vaccines and microbicides by repeated low-dose virus challenges. PLoS Med 2: e249.
class ConstantHazard( AnimalModel ): # p # constant probability p of infection on Challenge t
    def __init__(self, survivors):
        super().__init__(survivors)
    # Returns probability of infection corresponding to pt[0] = t-1.
    def p_infection(self, i, x):
        return x[0]
    # Returns one-dimensional list as MLE for reduced model.
    def lr_interval(self, confidence ):
        DF = 1
        p_hat = AnimalModel.x0( self.survivors )
        chi = chi2.ppf( confidence, DF )
        def diff( x ):
            return self.ln_likelihood( [x] ) - self.ln_likelihood( [p_hat] ) + 0.5 * chi
        if p_hat == 0.0:
            lo = 0.0
        else:
            lo = opt.brentq( diff, 0.0, p_hat )
        if p_hat == 1.0:
            hi = 1.0
        else:
            hi = opt.brentq( diff, p_hat, 1.0 )
        return [lo, hi]

    # Returns 2.0 * deviation of full from reduced model.
    def chisquared_fct(self):
        constantHazard = ConstantHazard( self.survivors )
        return 2.0 * (self.ln_likelihood( self.mle() ) - constantHazard.ln_likelihood( constantHazard.x0( self.survivors ) ))

    # Returns False if x violates bounds.
    @staticmethod
    def is_in_bounds(x):
        return 0.0 < x[0] <= 1.0
    # Returns True if x is on the boundary of the model parameter space.
    @staticmethod
    def is_null_on_boundary(x):
        return None # ConstantHazard is the null model.
    # Returns one-dimensional list as MLE for reduced model.
    @staticmethod
    def x0( survivors ):
        return [AnimalModel.x0( survivors )]

def _test_ConstantHazard():
    #print('constant_hazard')
    data = {
        ( 64, 32, 16, 8, 4, 2, 1 ):{'x':[0.5],'fun':-87.34,
            'fisher_information':[[504.]]}, # 504. checked by hand.
        ( 64, 16, 4, 1 ):{'x':[0.75],'fun':-47.24,
            'fisher_information':[[448.]]}, 
    }  
    for survivors,optimize_result0 in data.items():
        #print(survivors)
        x_hat0 = optimize_result0.get('x')
        fun0 = optimize_result0.get('fun')
        information0 = np.asarray(optimize_result0.get('fisher_information'))
        
        model = ConstantHazard( survivors )
        assert( model.name() == 'ConstantHazard' )
        x = [0.2]
        [ p ] = x
        for i in range(10):
            assert( model.p_infection(i, x) == p )
        
        assert( all( [ model.p_infection( i, x_hat0 ) == x_hat0[0] for i in range(10) ] ) )
        for method in _METHODS:
            #print(method)
            x_hat = model.mle( method )
            fun = model.ln_likelihood( x_hat )
            information = model.fisher_information( x_hat )
            assert( all( math.isclose( i, j, abs_tol=_EPS ) for i,j in zip( x_hat, x_hat0 ) ) )
            assert( math.isclose( information, information0, rel_tol=_EPS ) )
            assert( math.isclose( fun, fun0, rel_tol=_EPS ) )

class ConstantHazardFullModel( ConstantHazard ):
    # Returns 2.0 * deviation of full from reduced model.
    def chisquared_fct(self):
        constantHazard = ConstantHazard( self.survivors )
        return 2.0 * (self.ln_likelihood( self.mle() ) - constantHazard.ln_likelihood( constantHazard.x0( self.survivors ) ))
    # Returns p-value corresponding to the chisquared_fct.
    def df(self):
        return len( self.x0( self.survivors ) ) - len( super().x0( self.survivors ) )
    # Returns p-value corresponding to the chisquared_fct.
    def llr_pvalue(self):
        return chi2.sf(self.chisquared_fct(), self.df() )

# R.R. Regoes (2012) The role of exposure history on HIV acquisition: insights from repeated low-dose challenge studies. PLoS Comput Biol. 8: p. e1002767.
class ArithmeticPriming( ConstantHazardFullModel ): # p_infection = p + (t - 1) * eps on Challenge t
    def __init__(self, survivors):
        super().__init__(survivors)
    # Returns probability of infection corresponding to i = pt[0] = t-1.
    def p_infection(self, i, x):
        [ p, eps ] = x
        return p + i * eps
    # Returns False if x violates bounds.
    @staticmethod
    def is_in_bounds(x):
        return 0.0 < x[0] <= 1.0
    # Returns True if x is on the boundary of the model parameter space.
    @staticmethod
    def is_null_on_boundary(x):
        return False
    # Returns one-dimensional list as MLE for reduced model.
    @staticmethod
    def x0(survivors):
        return [AnimalModel.x0( survivors ), 0.0]

def _test_ArithmeticPriming():
    #print('arithmetic_priming')
    data = {
        ( 64, 32, 16, 8, 4, 2, 1 ):{'x':[0.5,0.0],'fun':-87.34,'llr_pvalue':1.0},
        ( 64, 16, 4, 1 ):{'x':[0.75,0.0],'fun':-47.24,'llr_pvalue':1.0}, 
    }  
    #print('ArithmeticPriming')
    for survivors,optimize_result0 in data.items():
        #print(survivors)
        x_hat0 = optimize_result0.get('x')
        fun0 = optimize_result0.get('fun')
        
        model = ArithmeticPriming( survivors )
        assert( math.isclose(  model.llr_pvalue(), optimize_result0.get('llr_pvalue'), abs_tol=_EPS ) )
        assert( model.name() == 'ArithmeticPriming' )
        x = [0.2, 0.1]
        [ p, eps ] = x
        for i in range(10):
            assert( model.p_infection(i, x) == p + i * eps )

        for method in _METHODS:
            #print(method)
            x_hat = model.mle( method )
            fun = model.ln_likelihood( x_hat )
            assert( all( math.isclose( i, j, abs_tol=_EPS ) for i,j in zip( x_hat, x_hat0 ) ) )
            assert( math.isclose( fun, fun0, rel_tol=_EPS ) )

class GeometricPriming( ConstantHazardFullModel ): # p_infection = p * r**(t - 1) on Challenge t
    def __init__(self, survivors):
        super().__init__( survivors )
    # Returns probability of infection corresponding to i = pt[0] = t-1.
    def p_infection(self, i, x):
        [ p, r ] = x
        return p * r ** i
    # Returns False if x violates bounds.
    @staticmethod
    def is_in_bounds(x):
        return 0.0 < x[0] <= 1.0 and 0.0 < x[1]
    # Returns True if x is on the boundary of the model parameter space.
    @staticmethod
    def is_null_on_boundary(x):
        return False
    # Returns one-dimensional list as MLE for reduced model.
    @staticmethod
    def x0( survivors ):
        return [AnimalModel.x0( survivors ), 1.0]

def _test_GeometricPriming():
    data = {
        ( 64, 32, 16, 8, 4, 2, 1 ):{'x':[0.5, 1.0],'fun':-87.34,
            'fisher_information':[[ 504.,228.],[228.,282.]],'llr_pvalue':1.0}, # 504. checked by hand.
        ( 64, 16, 4, 1 ):{'x':[0.75, 1.0],'fun':-47.24,
            'fisher_information':[[ 448.,96.],[96.,96.]],'llr_pvalue':1.0}, 
        ( 16384, 12288, 10752, 10080, 9765 ):{'x':[0.25, 0.5],'fun':-17758.51,
            'fisher_information':[[ 132139.4,33316.08],[33316.08,30196.32]],'llr_pvalue':0.0}, 
        ( 16, 12, 10, 10, 10 ):{'x':[0.2746, 0.3388],'fun':-15.18, # Nelder-Mead minimization
            'fisher_information':[[ 103.9577,22.89840],[22.89840,30.11120]],'llr_pvalue':0.01586106}, 
    }  
    #print('GeometricPriming')
    for survivors,optimize_result0 in data.items():
        #print(survivors)
        x_hat0 = optimize_result0.get('x')
        fun0 = optimize_result0.get('fun')
        information0 = np.asarray(optimize_result0.get('fisher_information'))

        model = GeometricPriming( survivors )
        assert( math.isclose(  model.llr_pvalue(), optimize_result0.get('llr_pvalue'), abs_tol=_EPS ) )
        #assert( math.isclose( model.llr_pvalue(), 0.0, abs_tol=_EPS ) )
        assert( model.name() == 'GeometricPriming' )
        x = [0.2, 0.1]
        [ p, r ] = x
        for i in range(10):
            assert( model.p_infection(i, x) == p * r ** i )
        
        for method in _METHODS:
            #print(method)
            x_hat = model.mle( method )
            fun = model.ln_likelihood(x_hat)
            information = model.fisher_information( x_hat )
            #print(fisher_information)
            assert( all( math.isclose( i, j, abs_tol=_EPS ) for i,j in zip( x_hat, x_hat0 ) ) )
            assert( np.allclose(information, information0, rtol=_EPS ) )
            assert( math.isclose( fun, fun0, rel_tol=_EPS ) )

# R.R. Regoes (2012) The role of exposure history on HIV acquisition: insights from repeated low-dose challenge studies. PLoS Comput Biol. 8: p. e1002767.
class StepPriming( ConstantHazardFullModel ): # p_infection  = p_1, but switches to p_2 strictly after Challenge l_step
    def __init__(self, survivors, l_step): # l_step is the time t at which p_2 starts to pertain.
        assert( isinstance( l_step, int ) and 0 < l_step )
        if len( survivors ) <= l_step:
            raise Exception('The change-point occurs after the end of challenges.')
        self.l_step = l_step
        super().__init__(survivors)
    # Returns probability of infection corresponding to i = pt[0] = t-1.
    def p_infection(self, i, x):
        [ p_1, p_2 ] = x
        if i < self.l_step:
            return p_1
        return p_2
    # Returns False if x violates bounds.
    @staticmethod
    def is_in_bounds(x):
        return 0.0 < x[0] <= 1.0 and 0.0 <= x[1] <= 1.0
    # Returns True if x is on the boundary of the model parameter space.
    @staticmethod
    def is_null_on_boundary(x):
        return False
    # Returns one-dimensional list as MLE for reduced model.
    @staticmethod
    def x0( survivors ):
        return [AnimalModel.x0( survivors ), AnimalModel.x0( survivors )]

def _test_StepPriming():
    data = {
        ( 64, 32, 16, 8, 4, 2, 1 ):{'x':[0.5,0.5],'fun':-87.34,'llr_pvalue':1.0},
        ( 64, 16, 4, 1 ):{'x':[0.75,0.75],'fun':-47.24,'llr_pvalue':1.0}, 
    }  
    #print('StepPriming')
    for survivors,optimize_result0 in data.items():
        #print(survivors)
        x_hat0 = optimize_result0.get('x')
        fun0 = optimize_result0.get('fun')
        
        for l_step in range(1,3):
            model = StepPriming( survivors, l_step )
            assert( math.isclose(  model.llr_pvalue(), optimize_result0['llr_pvalue'], abs_tol=_EPS ) )
            assert( model.name() == 'StepPriming' )
            x = [0.2, 0.1]
            [ p_1, p_2 ] = x
            for i in range(10):
                if i < l_step:
                    assert( model.p_infection(i, x) == p_1 )
                else:
                    assert( model.p_infection(i, x) == p_2 )
            for method in _METHODS:
                #print(method)
                x_hat = model.mle( method )
                fun = model.ln_likelihood(x_hat)
                assert( all( math.isclose( i, j, abs_tol=_EPS ) for i,j in zip( x_hat, x_hat0 ) ) )
                assert( math.isclose( fun, fun0, rel_tol=_EPS ) )

# R.R. Regoes (2012) The role of exposure history on HIV acquisition: insights from repeated low-dose challenge studies. PLoS Comput Biol. 8: p. e1002767.
class BetaFrailty( ConstantHazardFullModel ): # p # constant probability p of infection on Challenge t
    def __init__(self, survivors): 
        super().__init__(survivors)
    # Returns probability of infection corresponding to i = pt[0] = t-1.
    def p_infection(self, i, x): # x = [ p_mean, p_var ]
        [ a, b ] = BetaFrailty._to_beta_params( x )
        p_infection = a / (a + b + i)
        return p_infection
    # Returns False if x violates bounds.
    @staticmethod
    def is_in_bounds(x):
        return 0.0 < x[0] <= 1.0 and 0.0 <= x[1] <= x[0] * (1 - x[0])
    # Returns True if x is on the boundary of the model parameter space.
    @staticmethod
    def is_null_on_boundary(x):
        return True
    # Returns the first two centered moments for the beta distribution for the "beta_frailty" full model.
    @staticmethod
    def _to_moments( beta_params ): # [a,b] = beta_params
        [ a, b ] = beta_params 
        p_mean = a / (a + b)
        p_var = (a / (a + b)) * (b / (a + b)) / (a + b + 1.0)
        return [ p_mean, p_var ]
    # Returns [a,b] = beta_params for the beta distribution for the "beta_frailty" full model.
    @staticmethod
    def _to_beta_params( moments ): # [a,b] = beta_params
        [ p_mean, p_var ] = moments
        TOL = 1.0e-12
        s = p_mean * (1.0 - p_mean) / max( TOL, p_var ) - 1.0
        a = p_mean * s
        b = (1.0 - p_mean) * s
        return [ a, b ]
    # Returns one-dimensional list as MLE for reduced model.
    @staticmethod
    def x0( survivors ):
        p_mean0 = AnimalModel.x0( survivors )
        p_var0 = p_mean0 * (1.0 - p_mean0) * 0.1
        return [p_mean0, p_var0]

def _test_BetaFrailty():
    # test reparametrization
    [ a0, b0 ] = [ 3.0, 4.0 ]
    [ a, b ] = BetaFrailty._to_beta_params( BetaFrailty._to_moments( [ a0, b0 ] ) )
    assert ( abs( a / a0 - 1.0 ) < _EPS )
    assert ( abs( b / b0 - 1.0 ) < _EPS )
    data = { # Nelder-Mead minimization
        ( 64, 32, 16, 8, 4, 2, 1 ):{'x':[0.5, 0.0],'fun':-87.34,'llr_pvalue':1.0},
        ( 64, 16, 4, 1 ):{'x':[0.75, 0.0],'fun':-47.24,'llr_pvalue':1.0},
        ( 16384, 12288, 10752, 10080, 9765 ):{'x':[0.2534, 0.1114],'fun':-17821.39,
            'fisher_information':[[ 269904.7,-331621.3],[-331621.3,607597.8]],'llr_pvalue':0.0},
        ( 16, 12, 10, 10, 10 ):{'x':[0.2593, 0.1303],'fun':-15.71,
            'fisher_information':[[ 273.6344,-358.8308],[-358.8308,691.0599]],'llr_pvalue':0.02930025}
    }
    #print('BetaFrailty')
    for survivors,optimize_result0 in data.items():
        #print(survivors)
        x_hat0 = optimize_result0.get('x')
        [ p_mean0, p_var0 ] = x_hat0
        fun0 = optimize_result0.get('fun')

        model = BetaFrailty( survivors )
        assert( math.isclose(  model.llr_pvalue(), optimize_result0['llr_pvalue'], abs_tol=_EPS ) )
        assert( model.name() == 'BetaFrailty' )
        x = [0.2, 0.1]
        [ a, b ] = BetaFrailty._to_beta_params( x )
        for i in range(10):
            assert( model.p_infection(i, x) == a / (a + b + i) )
        
        for method in _METHODS:
            #print(method)
            x_hat = model.mle( method )
            fun = model.ln_likelihood(x_hat)
            [ p_mean, p_var ] = x_hat
            if p_var0 < _EPS: # boundary value
                assert( math.isclose( p_mean, p_mean0, rel_tol=_EPS ) )
                assert( math.isclose( p_var, p_var0, abs_tol=_EPS ) )
            else:
                information = model.fisher_information( x_hat )
                information0 = np.asarray(optimize_result0.get('fisher_information'))
                assert( all( math.isclose( i, j, rel_tol=_EPS ) for i,j in zip( x_hat, x_hat0 ) ) )
                assert( math.isclose( fun, fun0, rel_tol=_EPS ) )
                assert( np.allclose(information, information0, rtol=_EPS ) )

# M.G. Hudgens, et al. (2009) Power to detect the effects of HIV vaccination in repeated low-dose challenge experiments. J Infect Dis. 200: p. 609-13.
class DeltaFrailty( ConstantHazardFullModel ): # p # constant probability p of infection on Challenge t
    def __init__(self, survivors): 
        super().__init__( survivors )
    # Returns probability of infection corresponding to i = pt[0] = t-1.
    def p_infection(self, i, x): # x = [ p_mean, p_var ]
        [ p, theta ] = x
        p_infection = ((1.0 - theta) * p * (1.0 - p)**i) / (theta + (1.0 - theta) * (1.0 - p)**i)
        return p_infection
    # Returns False if x violates bounds.
    @staticmethod
    def is_in_bounds(x):
        return 0.0 < x[0] <= 1.0 and 0.0 <= x[1] < 1.0
    # Returns True if x is on the boundary of the model parameter space.
    @staticmethod
    def is_null_on_boundary(x):
        return True
    # Returns one-dimensional list as MLE for reduced model.
    @staticmethod
    def x0( survivors ):
        if not is_infection( survivors ):
            return [ 1.0, 1.0 ]
        survivor_count = survivors[-1]
        survivors0 = [ i - survivor_count for i in survivors ]
        p0 = AnimalModel.x0( survivors0 )
        theta0 = survivor_count / survivors[0]
        return [ p0, theta0 ]

def _test_DeltaFrailty():
    data = {
        ( 64, 32, 16, 8, 4, 2, 1 ):{'x':[0.5,0.0],'fun':-87.34,'llr_pvalue':1.0},
        ( 64, 16, 4, 1 ):{'x':[0.75,0.0],'fun':-47.24,'llr_pvalue':1.0}, 
        ( 16384, 12288, 10752, 10080, 9765 ):{'x':[0.5904, 0.5843],'fun':-17765.62,
            'fisher_information':[[28437.1,-7555.1],[-7555.1,64268.9]],'llr_pvalue':0.0},
        ( 16, 12, 10, 10, 10 ):{'x':[0.7397, 0.6232],'fun':-15.06,
            'fisher_information':[[ 35.61855,-1.804427],[-1.804427,67.64198697]],'llr_pvalue':0.01388016}
    }  
    #print('DeltaFrailty')
    for survivors,optimize_result0 in data.items():
        #print(survivors)
        x_hat0 = optimize_result0.get('x')
        fun0 = optimize_result0.get('fun')
        
        model = DeltaFrailty( survivors )
        assert( math.isclose(  model.llr_pvalue(), optimize_result0['llr_pvalue'], abs_tol=_EPS ) )
        assert( model.name() == 'DeltaFrailty' )
        x = [0.2, 0.1]
        [ p, theta ] = x
        for i in range(10):
            assert( model.p_infection(i, x) == ((1.0 - theta) * p * (1.0 - p)**i) / (theta + (1.0 - theta) * (1.0 - p)**i) )
        for method in _METHODS:
            #print(method)
            x_hat = model.mle( method )
            fun = model.ln_likelihood(x_hat)
            information = model.fisher_information( x_hat )
            assert( all( math.isclose( i, j, abs_tol=_EPS ) for i,j in zip( x_hat, x_hat0 ) ) )
            assert( math.isclose( fun, fun0, rel_tol=_EPS ) )
            if x_hat0[1] == 0.0: # The mle of full model is on the boundary.
                assert( np.all( np.isnan(information) ) )
            else:
                information0 = np.asarray( optimize_result0.get( 'fisher_information' ) )
                assert( np.allclose( information, information0, rtol=_EPS ) )

def main(): 
    _test_AnimalModel()
    _test_ConstantHazard()
    _test_ArithmeticPriming()
    _test_GeometricPriming()
    _test_StepPriming()
    _test_BetaFrailty()
    _test_DeltaFrailty()
    
if __name__ == "__main__":
    main()
