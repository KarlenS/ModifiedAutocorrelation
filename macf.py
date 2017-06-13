'''Module for calculating the modified autocorrelation function described in `Li (2009)`_. This method calculates temporal properties, e.g, power density, coherence, and time lag in the time domain, without utilizing Fourier transforms. The method is optimimal for studying rapid variablity on short time scales.

.. _Li (2009): http://adsabs.harvard.edu/abs/2001ChJAA...1..313L
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import time

def readData(filename):
    return pd.read_csv(filename, delim_whitespace=True,usecols=[0],names=['t'],header=0) 

def calculateMACF(lc):
    return calculateMCCF(lc,lc)


def delayLC(lc):

    print lc.index*0.1
    #return delayed_lc

def calculateMCCF(lc, delayed_lc):

    v1 = np.subtract(lc,np.mean(lc))
    v2 = np.subtract(delayed_lc,np.mean(delayed_lc))

    sigma1 = np.sum(np.square(v1))
    sigma2 = np.sum(np.square(v2))

    return v1*v2 / sigma1 / sigma2 #equation 20


def constructLCs(tte,dt, nbins = 20,delay=0):
    
    # force series to start at 0
    tte = tte.subtract(tte['t'].min())

    # introduce delay
    tte_delayed = tte.add(delay)

    # divisor to chop up light curves into chunks
    chunk_bounds = nbins * dt

    # create bins column following the required dt
    bins = dt*np.floor(tte['t']/dt)
    tte['bin'] = bins
    tte_delayed['bin'] = bins

    chunks = chunk_bounds*np.floor(tte['bin']/chunk_bounds)
    tte['chunk'] = chunks
    tte_delayed['chunk'] = chunks

    tte = tte.set_index('chunk')
    tte_delayed = tte_delayed.set_index('chunk')

    lc_list = []
    delayed_lc_list = []
    #handling the case where single-row chunks get turned series...
    for chunk in np.unique(tte.index):
        try:
            tte.ix[chunk].groupby('bin').size()
            tte_delayed.ix[chunk].groupby('bin').size()
            lc_list.append(tte.ix[chunk].groupby('bin').size())
            delayed_lc_list.append(tte_delayed.ix[chunk].groupby('bin').size())
        except ValueError:
            if isinstance(tte.ix[chunk]['bin'],np.float64):
                lc_list.append(tte.ix[chunk].to_frame().transpose().groupby('bin').size())
            elif isinstance(tte_delayed.ix[chunk]['bin'],np.float64):
                delayed_lc_list.append(tte_delayed.ix[chunk].to_frame().transpose().groupby('bin').size())
            else:
                raise ValueError('Single-valued series-to-frame conversion failed.')


    return [lc_list,delayed_lc_list]


    #return [ tte.ix[chunk].groupby('bin').size() for chunk in np.unique(tte.index) ]

    # create a light curve by counting up events in each bin
    #lc = tte.groupby('bin').size()
    

    # forcing uniform binning (othewise no bins for 0 counts) - this part might not be necessary... not sure
    #lc_bins = np.arange(tte['t'].min(),tte['t'].max(),dt)
    #zeros = np.zeros_like(lc_bins)
   
    # creating uniform bin dataframe with 0 counts and merging with the light curve
    # null_lc = pd.DataFrame(zeros,index = lc_bins,columns=['zeros'])
    # data_lc = pd.DataFrame(lc.values,index = lc.index.values,columns=['counts'])

    #print data_lc

    #full_lc = pd.concat([data_lc,null_lc],axis=1)
    #full_lc = full_lc.fillna(0)


    #print np.array_split( full_lc['counts'] , chunks ) 
    #print np.array_split(np.array(full_lc['counts']),3)

    #print np.array_split(full_lc['counts'],3)
    #returning light curve split into chunks
    #return np.array_split( full_lc['counts'] , chunks )


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description='Calculates the modified autocorrelation function (MACF), provided a file with time-tagged events.')
    parser.add_argument('-f',required=True,help='File containing the time-tagged event list.')
    args = parser.parse_args()

    tte = readData(args.f)


    # have light curves divided into chunks and binned for each dt
    # outer list is divided by dt, inner divided by chunks
    # 0 count bins are not included

    #all_lcs = [ constructLCs(tte,dt) for dt in np.arange(1,300,1.) ]
    #delayed_lcs = [ constructLCs(tte,dt,delay=10) for dt in np.arange(1,300,1.) ]

    
    for delay in np.arange(1.,300.,0.02):
        b,z = constructLCs(tte,1,delay=delay)
        if not (np.size(z[0].index) == np.size(b[0].index)):
            print 'fuck', delay, z[0].index, b[0].index

    #delayLC(all_lcs[0][0])

    print "|o|o|o|o| %s seconds |o|o|o|o|" % (time.time() - start_time)
     

    

    #lcs = constructLCs(tte,1.)

    #macf = calculateMACF()


if __name__ == '__main__':
    main()
