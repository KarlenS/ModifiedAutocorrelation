'''Module for calculating the modified autocorrelation function described in `Li (2009)`_. This method calculates temporal properties, e.g, power density, coherence, and time lag in the time domain, without utilizing Fourier transforms. The method is optimimal for studying rapid variablity on short time scales.

.. _Li (2009): http://adsabs.harvard.edu/abs/2001ChJAA...1..313L
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def readData(filename):
    return np.genfromtxt(filename,names=['t','e'])
    #return pd.readcsv(filename, delim_whitespace=True)

def calculateMACF(lc):
    return calculateMCCF(lc,lc)

def calculateMCCF(lc1, lc2):


    v1 = np.subtract(lc1,np.mean(lc1))
    v2 = np.subtract(lc2,np.mean(lc2))

    v1*v2 / sigma1 / sigma2 #equation 20

    return 


def constructLC(tte,dt):
    lc_start = np.min(tte)
    lc_end = np.max(tte)

    full_lc = np.arange(np.min(tte),np.max(tte),dt)

    return v

def plotMACF(macf):

    plt.show()

def main():

    parser = argparse.ArgumentParser(description='Calculates the modified autocorrelation function (MACF), provided a file with time-tagged events.')
    parser.add_argument('-f',required=True,help='File containing the time-tagged event list.')
    args = parser.parse_args()


    tte = readData(args.f)
    lc = prepLC(tte)

    macf = calculateMACF()


if __name__ == '__main__':
    main()
