# SpectralPeakFit
Spectral density peak finding and fitting using python scipy.

The peak finding and fitting functions in curfit.py are generic, can be used for any data.
The peak generator peakgen generates 3 peaks and adds the power spectrum density, generated from the white noise (psd-transformed noise). 

Example of testing using ipython:

    import curfit_test
    import matplotlib.pyplot as plt
    import numpy as np
    
    cft = curfit_test.cfTest()
    
    # process default peaks
    cft.run(block=False)

    # one may want to look and change the peak parametes:
    cft.print_peakgen_params()
    
    # change peaks and process them
    cft.generate_peaks([0.25, 0.01, 10e-9, 0.30, 0.01, 1e-9, 0.75, 0.01, 10e-9])
    cft.run(block=False)

    # loop with changing one peak
    cft.set_debug(0)
    plt.ion()
    while 1:
        rPos = np.random.random()*0.4+0.3
        rAmp = 1e-9*(0.1+5*np.random.random())
        rs2 = 0.005 + 0.025*np.random.random()
        cft.generate_peaks([0.25, 0.01,1e-9, rPos,rs2,rAmp, 0.75, 0.01, 1e-9])   
        cft.run(block=False)
        plt.pause(1)

