import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

def gauss(x, *p):
    a, mu, sigma = p
    #a, mu, sigma, d = p
    return a * np.exp( -(x-mu)**2 / (2.*sigma**2)) 
    #return a * np.exp( -(x-mu)**2 / (2.*sigma**2)) + d


#def getFWHM():

def fitGauss(x,y,*p):

    pars,cov = curve_fit(gauss, x, y, p0=p)
    perr = np.sqrt(np.diag(cov))
    return [pars[2],perr[2]]


def main():

    tscales = np.concatenate([np.arange(1,21,1),np.arange(25,125,5)])
    fig, ax = plt.subplots(5, 8, figsize=(25, 25))
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05,top=0.95)

    #fig.text(0.5, 0.015, 'Time Delay', ha='center', va='center',fontsize=16)
    #fig.text(0.015, 0.5, 'MACF', ha='center', va='center', rotation='vertical',fontsize=16)
    col = 1

    sns.set_style("ticks")
    fig2, ax2 = plt.subplots(1,1)
    for count,lvl in enumerate([0.0]):
    #for count,lvl in enumerate(np.arange(0,0.03,2E-3)):
        print "Processed lvl: ",lvl
        #print "Processed timescale: ",tscale[1]
        fwhm_norm = []
        times = []
        for countz,tscale in enumerate(zip(ax.flatten(),tscales)):
            ax_ts = tscale[0]
            ax_ts.annotate('t=%s' %str(tscale[1]),xy=(2,1),xytext=(120,0.95),fontsize=10)

            macfs = 'macfs/macf_%.3i.dat' %tscale[1]
            #macfs = 'macfs/macf_bglvl%s_%s_iter101.dat' %(str(lvl), str(tscale[1]))
            #print count%5
            data = np.genfromtxt(macfs,names=['t','macf'],comments='M')
            #print count/col
            #if count/col > 4:
            #    col += 1
            #ax_ts.set_xlabel(r'Time delay (s)',fontsize=14)
            #ax_ts.set_ylabel(r'MACF',fontsize=14)
            #ax_ts.set_ylim(-0.2,1.2)
            ax_ts.plot(data['t'],data['macf'])
            ax_ts.tick_params(axis='both', which='major', labelsize=6)

            p0 = [0.99, 0.0001, 3.]
            #p0 = [0.99, 0.0001, tscale[1]]
            pars = fitGauss(data['t'],data['macf'],*p0)
            fwhm_norm.append(pars[0]*2.355/tscale[1])
            times.append(tscale[1])
            #try:
            #    pars = fitGauss(data['t'],data['macf'],*p0)
            #    fwhm_norm.append(pars[0]*2.355/tscale[1])
            #    times.append(tscale[1])
            #except:
            #    continue
            #dt = (np.max(data['t']) - np.min(data['t'])) / 1000.
            #t_rs = np.arange(np.min(data['t']), np.max(data['t']), dt)


        if(count%3 == 0):
            ax2.plot(times,fwhm_norm,label = str(lvl))

    ax2.legend()
    #ax2.set_ylim(0.7,1.6)
    #ax2.savefig("macf_fits.png")
    #ax.savefig("macf_sims.png")
    plt.show()


if __name__ == '__main__':
    main()
