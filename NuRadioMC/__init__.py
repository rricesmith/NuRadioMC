""" NuRadioMC: Simulating the radio emission of neutrinos from interaction to detector"""
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(dir_path, '..'))
import NuRadioReco
__version__ = NuRadioReco.__version__
