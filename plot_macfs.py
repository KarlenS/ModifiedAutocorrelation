import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

def gauss(x, *p):
    a, mu, sigma = p
    #a, mu, sigma, d = p
    return a * np.exp( -(x-mu)**2 / (2.*sigma**2)) 
    #return a * np.exp( -(x-mu)**2 / (2.*sigma**2)) + d

def fitGauss(x,y,*p):

    pars,cov = curve_fit(gauss, x, y, p0=p)
    perr = np.sqrt(np.diag(cov))
    return [pars[2],perr[2]]


def main():

    #tscales = np.concatenate([np.arange(1,21,1),np.arange(25,125,5)])
    tscales = np.concatenate([np.arange(1,121,1),np.arange(125,305,5)])

    fig, ax = plt.subplots(5, 8, figsize=(25, 25))
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05,top=0.95)

    col = 1

    sns.set_style("ticks")
    fig2, ax2 = plt.subplots(1,1)
    fwhm_norm = []
    times = []
    #for countz,tscale in enumerate(zip(ax.flatten(),tscales)):
    for tscale in tscales: 
        #ax_ts = tscale[0]
        #ax_ts.annotate('t=%s' %str(tscale),xy=(2,1),xytext=(120,0.95),fontsize=10)
        #print tscale[1]

        macfs = 'macfs/macf_%.3i.dat' %tscale
        data = np.genfromtxt(macfs,names=['t','macf'],comments='M')
        #ax_ts.plot(data['t'],data['macf'])
        #ax_ts.tick_params(axis='both', which='major', labelsize=6)

        p0 = [0.99, 0.0001, 3.]
        pars = fitGauss(data['t'],data['macf'],*p0)
        fwhm_norm.append(pars[0]*2.355/tscale)
        times.append(tscale)


    ax2.plot(times,fwhm_norm)

    ax2.legend()
    plt.show()


if __name__ == '__main__':
    main()
