# SpectralPeakFit
Spectral density peak finding and fitting using python scipy.

The peak finding and fitting functions in curfit.py are generic, can be used for any data.
The peak generator peakgen generates 3 peaks and adds the power spectral density, generated from the white noise (psd-transformed noise). 

--- Prerequisites:

numpy, scipy - both are provided through Anaconda Distribution.
python2.7 - it looks like scipy.optimize.curve_fit() is not available in 2.6.

--- Performance of cuftit_test on aclinec using ipython:
filtering time: 0.060 ms.
peak finding time: 0.3 ms.
fitting time: 7.0 ms.

Typical result of processing 3 peaks:

![Alt text](curfit_result.png?raw=true "Result")

Example of testing using ipython:

    import matplotlib.pyplot as plt
    import numpy as np
    import curfit_test
    
    cft = curfit_test.cfTest()
    cft.run(block=False)

    # one may want to look and change the peak parameters:
    cft.print_peakgen_params()

    # loop with randomly changing one of the peaks 
    cft.set_debug(0)
    plt.ion()
    while 1:
        rPos = np.random.random()*0.4+0.3
        rAmp = 1e-9*(0.1+5*np.random.random())
        rs2 = 0.005 + 0.025*np.random.random()
        cft.generate_peaks([0.25, 0.01,1e-9, rPos,rs2,rAmp, 0.75, 0.01, 1e-9])   
        cft.run(block=False)
        plt.pause(1)


