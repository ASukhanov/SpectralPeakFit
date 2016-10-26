''' Peaks Generator. 
'''
__version__ = 'v01 2016-10-20' # adopted from curfit.py

import numpy as np
from scipy.optimize import curve_fit

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Peak-generating functions.
#
gPeaks = [] # flat array of peak parameters

def gaussian(xx, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((xx - mu)/sig, 2.)/2)
    
# Signal generator.
def generate_peaks(xx, pars = None):
    ''' Generate samples.'''
    global gPeaks
    if pars: 
        #print('peakgen setting peaks to: '+', '.join(['%0.2g' %i for i in pars]))
        gPeaks = pars
    #yy = sigFunction(xx)
    yy = np.zeros(len(xx))
    nPeaks = len(gPeaks)/ParsPerPeak
    for ii in range(nPeaks):
        scale = gPeaks[ii*ParsPerPeak + 2]/gaussian(0,0,gPeaks[ii*ParsPerPeak+1]) # to match the desired height
        yy += gaussian(xx,gPeaks[ii*ParsPerPeak],gPeaks[ii*ParsPerPeak+1])*scale
    yy += gBaseLine
    #print 'generate_peaks:',pars
    return yy

# Internal helper function.
def white_noise_psd(nSamples=201, samplingFrequency=0, powerDensity=1.e-6):
    ''' Generates samples of power spectrum density of white noise.'''
    # It probably could be simplified to just quadratic noise.
    # Note: sampling frequency only defines scale of the output frequencies
    if samplingFrequency == 0:
        samplingFrequency = nSamples*2
    noise_power = powerDensity * samplingFrequency / 2
    xx = np.random.normal(scale=np.sqrt(noise_power), size=nSamples*2-1)
    frequencies, powerDensities = signal.periodogram(xx, samplingFrequency)
    return frequencies, powerDensities

#  Smeared signal generator.
def smearSignal(signal,noiseLevel):
    ''' Signal corrupted with noise. '''
    #nP = len(signal)
    # Option 1: linear noise
    #return signal + np.random.normal(0,noiseLevel,len(signal)) # with linear noise
    #
    # Option 2: quadratic noise. 
    # Note the resulting noise will be not gaussian, the chisquare - underestimated
    #return signal + np.power(np.random.normal(0,noiseLevel,len(signal)),2)
    #    
    # Option 3: white noise power spectrum density
    f, psd = white_noise_psd(nSamples=len(signal), powerDensity=noiseLevel)
    # Bug?: the psd[0] is unrealisticly small (-400). Clone it from psd[1]:
    psd[0] = psd[1]
    return signal + psd
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def ind2par(indexes,nPoints):
    ''' convert indexes in range(nPoints) to function's argument space (0,1)'''
    #print index, np.array(index), np.array(index,dtype='float')/nPoints*2.
    return np.array(indexes,dtype='float')/nPoints

def par2ind(parameters,nPoints):
    ''' convert parameters in np.arange(0,1) to index range(nPoints)'''
    return (parameters)*nPoints
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def signalPower(yy):
    return np.sqrt(np.sum(np.power(yy,2)))/len(yy)

def rms(yy):
    return np.sqrt(np.mean(yy*yy))
    
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# The code below is for standalone testing
#
gNPoints = 1000
peakSigma = 0.01 # sigma of the main peak
sigAmp = 1.e-8 # amplitude of the main peak
gFloor = 1e-10
print('peakGen default signal amplitude:'+str(sigAmp)+', noise floor:'+str(gFloor))
gBaseLine = 0. # gBaseline have no sense for power spectral density plots as it is strictly 0.

#'' set for fitting of power spectral density plots
#sf = sigAmp/gaussian(0,0,peakSigma) 
ParsPerPeak = 3
#flat arrays of peak parameters
gPeaks = [0.25,peakSigma,sigAmp, 0.5,peakSigma,sigAmp, 0.75,peakSigma,sigAmp] 

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# The main function for debugging.
#
import matplotlib.pyplot as plt
from scipy import signal

# Y-scaling of plots. Uncomment one of 2 options:
# Option1: no scaling
#def scaleY(yy): return yy # no scaling
#yLabel='PSD[V**2/Hz]'
# option2: logarithmic scaling
def w2db(w):  return 10.*np.log10(np.array(w,dtype='double'))
def scaleY(yy): return w2db(yy) # decibels
yLabel='dBm'

if __name__ == "__main__":
# The code below is skipped if file is imported
    ix = np.arange(0, gNPoints, 1) 
    xx = ind2par(ix,gNPoints)
    print('Generating signal.')
    noiseFreeSignal = generate_peaks(xx)
    sigAmplitude = np.max(noiseFreeSignal)-gBaseLine
    print('Peaks with max amplitude generated:'+str(gPeaks))
    print('Signal amplitude:'+str(sigAmplitude)+', power:'+str(signalPower(noiseFreeSignal)))
    noise = smearSignal(np.zeros(gNPoints),gFloor)
    noiseRMS = np.std(noise)
    plt.plot(scaleY(noise),label='noise')
    print('NoiseFloor:'+str(gFloor)+', noise RMS:'+str(noiseRMS)+', power:'+str(signalPower(noise)))
    print('SNR:'+str(sigAmplitude/noiseRMS))
    while 1:
        plt.plot(scaleY(noiseFreeSignal),label='Original')
        #plot_original_samples()
        yy = smearSignal(noiseFreeSignal,gFloor)
        plt.plot(scaleY(yy),label='+ noise')
        legend = plt.legend(loc='upper right', shadow=True)
        plt.ylabel(yLabel)
        plt.xlabel('Frequency')
        if yLabel == 'dBm': plt.ylim([-110, -60])
        plt.show()


