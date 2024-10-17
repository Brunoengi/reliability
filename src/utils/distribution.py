# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 05:17:45 2024

@author: BrunoTeixeira
"""

from distribution.Beta import Beta
from distribution.Frechet import Frechet
from distribution.Gama import Gama
from distribution.Gumbel import Gumbel
from distribution.LogNormal import LogNormal
from distribution.Normal import Normal
from distribution.Uniform import Uniform
from distribution.Weibull import Weibull

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
    'distribution': Gumbel
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

def renameVariableDistribution(oldName):
        
    lowerOldName = oldName.lower()

    for new_name, synonyms in DISTRIBUTION_MAP.items():
        if lowerOldName in synonyms['names']:
            return new_name 

    raise ValueError(f"{oldName} is not a valid name")

def createDistribution(name, dictionaryInput):
    # Rename the distribution using the renameVariableDistribuition function
    new_name = renameVariableDistribution(name)
    
    # Get the corresponding distribution class from the dictionary
    distribution_class = DISTRIBUTION_MAP[new_name]['distribution']

    # Instantiate the distribution with the parameters
    return distribution_class(dictionaryInput)