''' Curve fitting, using curve_fit from scipy.optimize.
    The fitting function is defined in func(xx,*par).
'''
#__version__ = 'v03 2016-10-13' #  
#__version__ = 'v04 2016-10-14' # gPeaks is flat array now.
#__version__ = 'v05 2016-10-17' # Some cleanup. peak_finder replaces find3Peaks.
#__version__ = 'v06 2016-10-23' # Major cleanup. Testing moved to separate file.
#__version__ = 'v07 2016-10-26' # user defined function: amplitude2par added.
__version__ = 'v04 2016-10-27' # gaussian() replaced with peak_shape(), 


import numpy as np
from scipy.optimize import curve_fit


def peak_shape(xx, halfWidth):
    ''' Shape of the peak. The height is not important as it will be 
    calibrated using amplitude2par() function. 
    halfWidth is measure of the half width (sigma in the case of gaussian).'''
    #return 1./(np.sqrt(2.*np.pi)*halfWidth)*np.exp(-np.power(xx/halfWidth, 2.)/2) # full-fledged gaussian
    #return 1./(np.sqrt(2.*np.pi)*halfWidth)*np.exp(-((xx/halfWidth)**2)/2)
    # note: 1./(halfWidth)*np.exp(-(xx/halfWidth)**2) would be twice faster
    # to prove: timeit.timeit('import numpy as np; halfWidth,xx=1,1; 1./(halfWidth)*np.exp(-(xx/halfWidth)**2)', number=10000)
    return 1./(halfWidth)*np.exp(-(xx/halfWidth)**2) # fitting time dropped from 12 ms to 7 ms!
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Fitting function, user modified.

def func(xx, *par):
    ''' Fitting function'''
    nx = len(xx)
    nPeaks = 3
    sumy = np.zeros(nx)
    for ii in range(nPeaks):
        yy = peak_shape(xx-par[nPeaks*ii+1],par[nPeaks*ii+2])*par[nPeaks*ii+3]
        sumy += yy
    sumy += par[0]
    return sumy
    
def amplitude2par(yy, sigma): #convert amplitude to function parameters of the single peak
    return yy/peak_shape(0,sigma)
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Helper functions
# 
def db2w(dB): 
    ''' Decibel to linear conversion. '''
    return np.power(10.,np.array(dB,dtype='double')/10.)

def w2db(w):
    ''' Linear to decibel conversion '''
    return 10.*np.log10(np.array(w,dtype='double'))

def ind2par(indexes,nPoints):
    ''' Convert indexes in range(nPoints) to function's argument space (0,1)'''
    #print index, np.array(index), np.array(index,dtype='float')/nPoints*2.
    return np.array(indexes,dtype='float')/nPoints

def par2ind(parameters,nPoints):
    ''' Convert parameters in np.arange(0,1) to index range(nPoints)'''
    return (parameters)*nPoints

def signalPower(yy):
    return np.sqrt(np.sum(yy*yy))/len(yy)

def rms(yy):
    return np.sqrt(np.sum(yy*yy)/len(yy))
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Curve fitting
#
def curfit(xx,yy,noise_rms,guess=[]):
    #gNIterations = 0
    popt, pcov = curve_fit(func, xx, yy, p0=guess)    
    degrees_of_freedom = len(yy)-len(guess)
    reduced_chi_squared = np.sum(((func(xx,*popt) - yy)/noise_rms)**2)/degrees_of_freedom
    if np.abs(reduced_chi_squared - 1) > 0.5:
        print('WARNING! ######### Bad fit! chisq:'+str('%0.2g' % reduced_chi_squared)+
          '. Is noise estimate correct? I suggest '+
          str('%0.2g' % (w2db(noise_rms) + w2db(reduced_chi_squared)/2.))+' dBm')
    stdevs = np.sqrt(np.diag(pcov))
    return reduced_chi_squared, popt, stdevs, func(xx,*popt)
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Peak finding
#
def peak_finder(yy,minWidth=10,maxNumberOfPeaks=3):
    ''' Finds positions of peaks of minimal width.'''
    # Returns all local maxima of yy array, the peak is discarded if it is 
    # within minWidth distance from a bigger peak.
    # For better peak detection the yy array should be match-filtered with 
    # desired template using scipy.signal.correlate(yy,template,mode='same').
    pp=[]
    previ,prevp = -minWidth,-9.e99
    #nn = 0
    cur,next = prevp, yy[0] # the cycling works twice faster
    for ii in range(0,len(yy)-1):
        #if (yy[ii]>yy[ii-1]) & (yy[ii]>=yy[ii+1]):
        prev,cur,next = cur,next,yy[ii+1]                    
        if (prev<cur) & (cur>=next):
            #print 'local peak ', yy[ii], ' at ',ii, nn
            #if (ii - pp[nn][0]) > minWidth:
            if (ii - previ) >= minWidth:
                 pp.append((ii,cur))
                 #nn +=1
                 previ,prevp = (ii,cur)
                 #print 'inserted ',pp[nn-1][0]
            else:
                 #print 'too close to prev ', pp[nn][0]
                 #if yy[ii] > pp[nn][1]:
                 if cur > prevp:
                     #print 'but it is higher, replace the old one'
                     try:
                        pp[len(pp)-1] = (ii,cur)
                     except:
                        #print('peak_finder trying to replace [',+str(len(pp))+'], ii:'+str(ii))
                        print 'ERROR peak_finder trying to replace ',pp,ii,minWidth
                        print prev, cur, next, previ, prevp
                     previ,prevp = (ii,cur)
    #print 'peakPositions:', pp 
    pp.sort(key=lambda x: x[1],reverse=True)
    #print 'sorted:',pp
    peakPositions = map(lambda x: x[0],pp)[:maxNumberOfPeaks]
    #print peakPositions
    peakPositions.sort()
    return peakPositions
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

