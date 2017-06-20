import subprocess
import os
import time
import errno
import numpy as np
import argparse



def submit_job(macf_func = 'macf.py',macfDir='macfs',logDir='logs',args='',tscale=1.):

    #create directory for logs, if it doesn't already exist
    try:
        os.makedirs(logDir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(logDir):
            pass
        else:
            raise IOError

    #create directory for logs, if it doesn't already exist
    try:
        os.makedirs(macfDir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(logDir):
            pass
        else:
            raise IOError

    logFile = '%s/macf_%.3i' % (logDir,tscale)

    submitOut = subprocess.Popen(['condor_submit', 'executable=%s' %macf_func, 
                                    'arguments=%s' %args,'log=%s.condor.log' %logFile, 
                                    'output=%s.out' %logFile,'error=%s.out' %logFile, 
                                    'getenv=true', '/data/gammaray/users/shahin/Mrk421_autocorr/macfsignaloptimizationmc/vegas.submit'],stdout=subprocess.PIPE)

    job_submit, error = submitOut.communicate()
    #return [job_submit,error]


def main():

    parser = argparse.ArgumentParser(description='Wrapped for submitting MACFs to Lucifer clusters.')
    parser.add_argument('-f',required=True,help='File containing the time-tagged event list.')
    args = parser.parse_args()


    tscales = np.concatenate([np.arange(1,121,1),np.arange(125,305,5)])
    #tscales = np.concatenate([np.arange(1,21,1),np.arange(25,125,5)])

    for dt in tscales:
        macf_args = '-f %s -t %s' %(args.f,dt)
        submit_job(args=macf_args,tscale = dt)

    print '%s jobs submitted.' %np.size(tscales)

        #print "|o|o|o|o| delta_t %s done in %s seconds |o|o|o|o|" % (dt,time.time() - start_time)

    #print "|o|o|o|o| Final length: %s seconds |o|o|o|o|" % (time.time() - start_time)




if __name__ == '__main__':
    main()
