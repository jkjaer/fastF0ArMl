# fastF0ArMl - Fast Maximum Likelihood Estimation of the Fundamental Frequency in AR noise
This git repository contains MATLAB-implementations of fast algorithms for fundamental frequency estimation in autoregressive (AR) noise. The details of the algorithms can be found in

- Fast Algorithms for Fundamental Frequency Estimation in Autoregressive Noise. Quinn, Barry Gerard; Nielsen, Jesper Kjær; Mads Græsbøll; Signal Processing, submitted, 2020

Moreover, the git repository also contains all the code and data needed for reproducing the results in the paper. Finally, an example is also given to show how the estimator is used.

## How to reproduce the results in the paper
First, go to the `code/experiments` folder. There, you will find five folders, each corresponding to a section in the experimental part of the paper. In each folder, you can find a number of files. The files starting with a lower-case `s` are scripts and should be run in order to generate the results. The files starting with a lower-case `sPlot` are used for plotting the results.

Please note that the it can take several days to recompute the computation times and estimation accuracy results, depending on the hardware you have access to.

## Using the estimator on your own data
Before you use the estimator on your own data, please do the following:
1. Make the example in the `example` folder run on your own computer.
2. Read the documentation in the class-file `code/estimator/@FastF0ArMl/FastF0ArMl.m`.

## Scope and limitations
Note that the estimator can be refined in a number of ways. First of all, a fundamental frequency tracker should be added on top of the estimator for segment-by-segment processing. The tracker is there to minimise the risk of, e.g., octave errors (aka pitch halving) by smoothing out the estimates produced by the estimator which analyse the data on a segment-by-segment basis. For the best performance, we recommend that the tracker by either Shi et al. or Tabrikian et al. is used. See more in

- Robust Bayesian Pitch Tracking Based on the Harmonic Model. Shi, L.; Nielsen, J. K.; Jensen, J. R.; Little, M. A.; Christensen, M. G. In: IEEE/ACM Transactions on Audio, Speech, and Language Processing. 27 (11), 2019, pp. 1737 - 1751.
- Maximum a-posteriori probability pitch tracking in noisy environments using harmonic model, Tabrikian, J.; Dubnov, S.; Dickalov, Y. In: IEEE Transactions on Speech and Audio Processing 12 (1), 2004, pp. 76-87.

Alternatively, and much simpler, median smoothing or dynamic programming can also be used. See more in
- Postprocessing techniques for voice pitch trackers,  Secrest, B.; Doddington, G.; in Proc. IEEE Int. Conf. Acoust., Speech, Signal Process. IEEE, 1982, vol. 7, pp. 172--175.
- Pitch and voicing determination of speech with an extension toward music signals, Hess, W. J., Springer handbook of speech processing. Springer, Berlin, Heidelberg, 2008. 181-212.

Although the estimator does not have many user parameters, it should be set-up for the data that are going to be processed. Most importantly the range of valid fundamental frequencies must be specified (a class property called `pitchBounds`). Also, the valid number of harmonic components must be tuned to the data (a class property called `maxPitchOrder`).

Please note that we have only focused on making the computation of the residual mean square for all candidate sinusoidal and AR-orders fast. However, the linear parameters for the estimated model and fundamental frequency are just computed naïvely. These can also be computed more efficiently in the applications where they are needed.

## License

The code is published under a BSD 3-clause "New" or "Revised" License. This is a very permissive license that enables the code to be used free of charge in both academia and industry. Besides honouring the license, we therefore also expect that you 

- cite our work if you use the code academically, and
- notify us if you use the code commercially.

Constructive feedback in any form is also very much appreciated. This applies to everything from bugs to data where the estimator fails.

The code has been written and is maintained by

- Assoc. Prof. Jesper Kjær Nielsen, Audio Analysis Lab, CREATE, Aalborg University, Denmark
