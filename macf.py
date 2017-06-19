#!/local/gammasoft/anaconda2/bin/python
import time
import numpy as np
import pandas as pd
import argparse

class TimeDomainAnalysis(object):

    '''Module for calculating the modified autocorrelation function described in `Li (2009)`_. This method calculates temporal properties, e.g, power density, coherence, and time lag in the time domain, without utilizing Fourier transforms. The method is optimimal for studying rapid variablity on short time scales.
    
    .. _Li (2009): http://adsabs.harvard.edu/abs/2001ChJAA...1..313L
    '''


    def __init__(self,filename):
        self.filename = filename
        self.mccf_outfile = None

    def readData(self):
        return pd.read_csv(self.filename, delim_whitespace=True,usecols=[0],names=['t'],header=0) 

    def writeMCCF(self,mccf):
        np.savetxt(self.macf_outfile,np.transpose(mccf),fmt='%.4e')
    
    def calculateMACF(self,base_lc,del_lc):
    
        v = base_lc.subtract(base_lc.mean())
        v_del = del_lc.subtract(del_lc.mean())
        var = v.apply(lambda x: x**2).sum()
        macf = v.multiply(v_del)/var
    
        return macf.sum()
    
    def calculateMCCF(self,base_lc, del_lc):
    
        v = base_lc.subtract(base_lc.mean())
        v_del = del_lc.subtract(del_lc.mean())
        sigma1 = np.sqrt(v.apply(lambda x: x**2).sum())
        sigma2 = np.sqrt(v_del.apply(lambda x: x**2).sum())
        mccf = v.multiply(v_del)/sigma1/sigma2
    
        return mccf.sum() 
    
    
    def processLCs(self,tte,dt, nbins = 20,delay_res=0.02):

        self.macf_outfile = 'macf/macf_%.2f.dat' %dt
        
        # force series to start at 0
        tte = tte.subtract(tte['t'].min())

        # create bins column following the required dt
        tte['bin'] = dt*np.floor(tte['t']/dt)

        # get longer individual light curves
        if np.abs(dt) < 10:
            nbins = 100
        chunk_bounds = nbins * dt

        chunk = chunk_bounds*np.floor(tte['bin']/chunk_bounds)
        tte['chunk'] = chunk 
    
        tte = tte.set_index('chunk')

        delays = np.arange(0,dt*nbins,delay_res*dt) - dt*nbins
       
        ccfs = np.zeros_like(delays)

        for count,delay in enumerate(delays):
            print 'on delay %s' %delay
            # introduce delay
            tte['t_delayed'] = tte['t'].add(delay)
    
            # divisor to chop up light curves into chunks
            # higher dt will have longer but fewer chunks
    
            tte['delayed_bin'] = dt*np.floor_divide(tte['t_delayed'],dt)
    
            mean_ccf = 0.
            achunks = np.unique(tte.index)
            nchunks = np.size(achunks)
            #handling the case where single-row chunks get turned series...
            for chunk in achunks:
                try:
                    base_lc = tte.ix[chunk].groupby('bin').size()
                    del_lc = tte.ix[chunk].groupby('delayed_bin').size()
                    #running average of ccfs
                    mean_ccf += self.calculateMACF(base_lc,del_lc) / nchunks
                except ValueError:
                    if isinstance(tte.ix[chunk]['bin'],np.float64):
                        base_lc = tte.ix[chunk].to_frame().transpose().groupby('bin').size()
                    elif isinstance(tte.ix[chunk]['delayed_bin'],np.float64):
                        del_lc = tte.ix[chunk].to_frame().transpose().groupby('delayed_bin').size()
                    else:
                        raise ValueError('Single-valued series-to-frame conversion failed.')
                    mean_ccf += self.calculateMACF(base_lc,del_lc) / nchunks

            ccfs[count] = mean_ccf
    
        return np.array([delays,ccfs])

def main():

    # have light curves divided into chunks and binned for each dt
    # outer list is divided by dt, inner divided by chunks
    # 0 count bins are not included


    parser = argparse.ArgumentParser(description='Calculates the modified autocorrelation function (MACF), provided a file with time-tagged events.')
    parser.add_argument('-f',required=True,help='File containing the time-tagged event list.')
    parser.add_argument('-t',type=float,required=True,help='Time scale to search for variability.')
    args = parser.parse_args()

    start_time = time.time()

    tda = TimeDomainAnalysis(args.f)

    tte = tda.readData()

    nbins = 20
    delay_res = 0.02

    mccf = tda.processLCs( tte,args.t,delay_res=delay_res,nbins=nbins )
    tda.writeMCCF(mccf)

    print 'Finished running in %s seconds.' %(time.time() - start_time)



if __name__ == '__main__':
    main()
