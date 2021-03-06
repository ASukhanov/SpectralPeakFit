''' Curve fitting test program
'''
#__version__ = 'v01 2016-10-23' # Extracted from curfit.py
#__version__ = 'v02 2016-10-25' # gen_peaks parameters, print helpers, ipython example.
#__version__ = 'v03 2016-10-26' # scaleY: scale parameter added, default:'dBm'.
__version__ = 'v05 2016-10-27' # cf.gaussian replaced with cf.peak_shape, removed mu parameter in peak_shape

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import signal
from  timeit import default_timer as timer

import curfit as cf
import peakgen as pg

print('curfit_test version:'+str(__version__))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# This section is for standalone testing
#
gNPoints = 1000
gSNR = 100 # desired signal/noise of the generated signal
peakSigma = 0.05 # sigma of the main peak
sigAmp = 1.e-8 # amplitude of the main peak
gFloor = sigAmp/gSNR
print('curfit_test default signal amplitude:'+str(sigAmp)+', noise floor:'+str(gFloor))
gBaseLine = 0. # gBaseline have no sense for power spectral density plots as it is strictly 0.

#'' set for fitting of power spectral density plots
sf = sigAmp/cf.peak_shape(0,peakSigma) 
ParsPerPeak = 3
gGuess = [gBaseLine, 0.2,peakSigma*0.8,0.5*sf, 0.45,peakSigma*.8,sf, 0.75,peakSigma*0.8,0.5*sf]
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
class cfTest():

    def set_debug(self,vv): self.dbg = vv

    def generate_peaks(self, pars=None): # generate peaks
        # to change the peaks:
        #if pars : pg.gPeaks = pars
        if self.dbg & 0x1: print('Generating signal['+str(self.nX)+'].')
        self.noiseFreeSignal = pg.generate_peaks(self.xx,pars)
        self.sigAmplitude = np.max(self.noiseFreeSignal)-gBaseLine
        if self.dbg & 0x1: 
            print('Peaks with max amplitude generated:'+str(self.sigAmplitude))
            print('Signal amplitude:'+str(self.sigAmplitude)+' power:'+str(cf.signalPower(self.noiseFreeSignal))+' RMS:'+str(cf.rms(self.noiseFreeSignal)))
        
    def print_peakgen_params(self):
        print('peakgen parameters:['+self.dlm.join(['%0.2g' %i for i in pg.gPeaks])+']')
    def print_guess_params(self):
        print('guess parameters:['+self.dlm.join(['%0.2g' %i for i in self.guess])+']')
    def print_fitted_params(self):
        print('guess parameters:['+self.dlm.join(['%0.2g' %i for i in self.popt])+']')
        
    def new_guess(self, pars=gGuess):  # Change all guess prameters.
        self.guess = pars
        print('guess:'+', '.join(['%0.2g' %i for i in self.guess]))
        
    def change_guess_sigmas(self, vv): # Change guess sigmas
        self.guess[2], self.guess[5], self.guess[8] = vv,vv,vv        
        
    def change_filterWidth(self, fw=10): #
        self.filterWidth = cf.ind2par(fw,self.nX) - cf.ind2par(0,self.nX)
        self.change_guess_sigmas(self.filterWidth)
        print 'filterWidth and guess sigmas changed to ', self.filterWidth
        
    def scaleY(self,yy, scale='dBm'): # logarithmic Y scale for plotting
    #def scaleY(self,yy, scale=''): # linear Y scale for plotting
        ''' Y-scaling of plots, uncomment one of two options: linear or logarithmic.'''
        if scale == 'dBm':
            self.yLabel='dBm'
            yy = cf.w2db(yy)
        else:
            self.yLabel='PSD[V**2/Hz]'
        return yy
           
    def run(self,block=True):
        plt.close()
        plt.plot(self.scaleY(self.noiseFreeSignal),label='Original')
        
        # corrupt the signal with noise
        yy = pg.smearSignal(self.noiseFreeSignal,gFloor)
        plt.plot(self.scaleY(yy),label='+ noise')

        # filter the signal using cross-correlation with the template
        #plt.plot(self.template+gBaseLine, label='self.template')
        start = timer()
        filtered = signal.correlate(yy-gBaseLine,self.template,mode='same')*self.filterScale
        end = timer()
        if self.dbg & 0x1: print('filtering time '+str(end-start)+' s.') # ~0.08ms
        #print('filtered: '+' power:'+str(cf.signalPower(filtered))+' RMS:'+str(rms(filtered)))
        plt.plot(self.ix,self.scaleY(filtered+gBaseLine), label='Filtered')

        # find peak positions and use it as a seed for the curfit
        start = timer()
        peaks = cf.peak_finder(filtered,self.filterWidth)
        end = timer()
        if self.dbg & 0x1: 
            print('peak finding time:'+str(end-start)+' s.') # ~0.5-1 ms 
            print('Peaks found at x:'+str(cf.ind2par(peaks,self.nX))+', or samples:'+str(peaks)+
                ', amplitudes:['+ ', '.join(['%0.2g' % i for i in filtered[peaks[:3]]])+']')

        # adjust guessed positions
        self.guess[1], self.guess[4], self.guess[7] = cf.ind2par(peaks[:3],self.nX)
        # adjust guessed amplitudes
        for ii in range(3):
            self.guess[1+ii*3+2] = cf.amplitude2par(filtered[peaks[ii]],self.guess[1+ii*3+1])
        print self.guess

        # fitting
        try:
            start = timer()
            reduced_chi_squared, self.popt, stdevs, yf = cf.curfit(self.xx,yy,self.noiseSTD,guess=self.guess)
            end = timer()
            if self.dbg & 0x1: 
                print('fitting time:'+str(end-start)) # ~15 ms
                print('chi_sq/dof:'+str('%0.2g' % reduced_chi_squared))
                print('params:'+str(['%0.2g' % i for i in self.popt]))
                print('param stdevs:'+str(["%0.2g" % i for i in stdevs[1:]]))
            plt.plot(self.ix,self.scaleY(yf),label='Fitted')    
        except Exception, e:
            print('EXCEPTION curfit failed')
            print(str(e))
            yf = np.zeros(self.nX)
            reduced_chi_squared = 999999.

        # plot everything
        plt.plot(self.ix,self.scaleY(cf.func(self.xx, *self.guess)),label='Guess')
        if block == True: print('Close the plot to get the next one.')
        legend = plt.legend(loc='upper right', shadow=True)
        plt.ylabel(self.yLabel)
        plt.xlabel('Frequency')
        if self.yLabel == 'dBm': plt.ylim([-110, -60])
        plt.show(block)
        
    def __init__(self):
        # The code below is skipped if file is imported
        self.dbg = 1
        self.dlm = ', ' # delimiter for printing
        self.ix = np.arange(0, gNPoints, 1) 
        self.nX = len(self.ix)
        self.xx = cf.ind2par(self.ix,self.nX)
        self.generate_peaks()
        noise = pg.smearSignal(np.zeros(self.nX),gFloor)
        self.noiseSTD = np.std(noise)
        print('Noise: floor:'+str(gFloor)+' power:'+str(cf.signalPower(noise))+' RMS:'+str(cf.rms(noise)))
        print('SNR:'+str(self.sigAmplitude/self.noiseSTD))

        self.new_guess()
        # Create filter. Crop out 10% tails to speed up cross-correlation.
        #fullFilter = np.ones(100)/100 # flat filter
        self.change_filterWidth()
        fullFilter = cf.peak_shape(self.xx-0.5,self.filterWidth)
        croppedArea = 0.1 * np.sum(fullFilter)
        fl = len(fullFilter)
        fHW = 0 # filter halfwidth
        while croppedArea > 0:
          croppedArea -= fullFilter[fHW] + fullFilter[fl-1-fHW]
          fHW += 1
        fHW = fl/2 - fHW
        template = fullFilter[fl/2-fHW:fl/2+fHW]
        self.template = fullFilter[fl/2-fHW:fl/2+fHW]
        print('Filter: sigma='+str(self.filterWidth)+' or '+
          str(cf.par2ind(self.filterWidth,self.nX))+' samples, sum='+str(np.sum(fullFilter))+
          ', template['+str(2*fHW)+'] sum='+str(np.sum(self.template)))
        self.filterScale = 1./(np.sum(self.template))
        print('filter: scale:'+str(self.filterScale)+' power:'+str(cf.signalPower(self.template))+' RMS:'+str(cf.rms(self.template)))

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# The main function for debugging.
#
if __name__ == "__main__":
    # The code below is skipped if file is imported.
    # The best way is to run it in ipython.
    # ipython example:
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import curfit_test
    
    cft = curfit_test.cfTest()
    cft.run(block=False)

    # one may want to look and change the peak parametes:
    cft.print_peakgen_params()

    # quite impressive peak fitting:
    cft.generate_peaks([0.25, 0.01, 10e-9, 0.30, 0.01, 1e-9, 0.75, 0.01, 10e-9])
    cft.run(block=False)

    # loop
    cft.set_debug(0)
    plt.ion()
    while 1:
        rPos = np.random.random()*0.4+0.3
        rAmp = 1e-9*(0.1+5*np.random.random())
        rs2 = 0.005 + 0.025*np.random.random()
        cft.generate_peaks([0.25, 0.01,1e-9, rPos,rs2,rAmp, 0.75, 0.01, 1e-9])   
        cft.run(block=False)
        plt.pause(1)
        
    '''
    # for pure python:
    cft = cfTest()
    while 1:
        rPos = np.random.random()*0.4+0.3
        rAmp = 1e-9*(0.1+5*np.random.random())
        rs2 = 0.005 + 0.025*np.random.random()
        #print rPos,rAmp
        cft.generate_peaks([0.25, 0.01,1e-9, rPos,rs2,rAmp, 0.75, 0.01, 1e-9])   
        cft.run()

