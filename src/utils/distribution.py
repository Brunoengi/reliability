# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 05:17:45 2024

@author: BrunoTeixeira
"""

from src.distribution import Normal, Beta, Frechet, Gama, Gumbell, LogNormal, Uniform, Weibull

DISTRIBUTION_MAP = {
  'gauss': {
    'names': ['norm', 'normal', 'gauss'],
    'distribution': Normal
      },
  'uniform': {
    'names': ['uniform', 'uniforme', 'const'],
    'distribution': Uniform
  },
  'lognorm': { 
    'names': ['lognormal', 'lognorm', 'log'],
    'distribution': LogNormal,
  },
  'gumbel': {
    'names': ['gumbel', 'extvalue1', 'evt1max'],
    'distribution': Gumbell  
  },
  'frechet': {
    'names': ['frechet', 'extvalue2', 'evt2max'],
    'distribution': Frechet
  },
  'weibull': {
    'names': ['weibull', 'extvalue3', 'evt3min'],
    'distribution': Weibull
  },
  'beta': {
    'names': ['beta', 'beta_dist'],
    'distribution': Beta
  },
  'gamma': {
    'names': ['gamma', 'gama'],
    'distribution': Gama
  }
}

def renameVariableDistribuition(oldName):
        
    lowerOldName = oldName.lower()

    for new_name, synonyms in DISTRIBUTION_MAP.items():
        if lowerOldName in synonyms['names']:
            return new_name 

    raise ValueError(f"{oldName} is not a valid name")

def createDistribution(name):
    new_name = renameVariableDistribuition(name)
    distribution_class = DISTRIBUTION_MAP[new_name]['distribution']
    return distribution_class