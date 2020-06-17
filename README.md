# Continuous Media Monte Carlo Transporter

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3897410.svg)](https://doi.org/10.5281/zenodo.3897410)
[![License](https://img.shields.io/badge/License-CeCILL%20v2.1-brightgreen)](https://github.com/HunterBelanger/continuous-media-monte-carlo-transporter/blob/master/LICENSE)

This repository contains a brief code which aims to compare different transport
methods for solving the Boltzmann transport equation when the total cross
section is spatially non-uniform.

Please refer to the following article for more information:
"H. Belanger et al., submitted to EPJ+ (2020)"

The example output files in the ```data``` directory contain the data sets used
to write this article. A python script is also provided there to plot the
collision density and FOM for all problems.

## Problem Description
The code uses as one-dimensional rod geometry. The media is located in within
the domain 0 <= x <= 2, with vacuum boundary conditions imposed at x=0 and x=2.
Particles all enter the rod at x=0, moving in the positive direction. Despite
the total cross section varying as a function of position within the rod,
individual reaction channel probabilities are kept constant at all locations.
The probability of absorption is 0.3, and that of scatter is 0.7. Scattering is
isotropic, with equal probability of scattering in the forward or backwards
direction.

Implicit capture and Russian Roulette are used for variance reduction. The
survival weight is 1 while the roulette cutoff weight is 0.6. Particle
spliting is also used as a variance reduction technique for the transport
methods which augment particle weights. Particles with a weight w greater than
or equal to 2 are split into N identical particles, each with weight w/N,
where N is the rounded value of the original weight.

## Transport Methods
The transport methods evaulated in this code are based on the following works:

* Direct Sampling : F.B. Brown, W.R. Martin, *Direct Sampling of Monte Carlo
  Flight Paths in Media with Continuously Varying Cross-Sections*, in *ANS
  Mathematics & Computation Topical Meeting* (2003), Vol. 2

* Delta Tracking : E.R. Woodcock, T. Murphy, P.J. Hemmings, T.C. Longworth,
  Technical report. ANL-7050, ArgonneNational Laboratory, USA (1965) 

* Negative Weighted Delta Tracking : D. Legrady, B. Molnar, M. Klausz, T. Major,
  Ann. Nucl. Energy **102**, 116 (2017)

* Carter, Cashwell, and Taylor Tracking : L.L. Carter, E.D. Cashwell, W.M. Taylor,
  Nucl. Sci. Eng. **48**, 403 (1972) 

## License
This software is distributed under the terms and conditions of the CeCILL v2.1
license. The English version is provided in the LICENSE file, and the equally
valid and authenic French version is provided in LICENSE-FRENCH. This is a
copyleft license, very simillar to the GNU GPL, and is compatible with both
the GPLv2 and GPLv3. If you would like more information about the CeCILL
license beyond that provided in the two license files, please refer to
the [CeCILL website](http://cecill.info/index.en.html).

## Install
To build and install cmmct, run the following commands:
```
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/for/binary ..
$ make
```
This will build a ```cmmct``` executable, with OpenMP support automatically
enabled. If you would like to dissable OpenMP, run cmake with the 
```-Dopenmp=off``` option.
