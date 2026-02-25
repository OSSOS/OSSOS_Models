# OSSOS_Models
Nominal models of the Outer Solar System populations derived from OSSOS++ sample

---
## OSSOS Models 1.0 Documentation

OSSOS Models provide calibrated sets of orbits and absolute magnitudes for the various populations of the small bodies of the Outer Solar System:
Classical, Detached, Distant resonances, Inner, Plutinos, Scattering, Twotinos. These models are in agreement with the detections in the OSSOS++
sample (_Bannister et al (2018) ApJS, 236, 18_).

---
### Licence
This software is released under the terms of the European Union Public
Licence v1.1 (EUPL v.1.1). See files eupl1.1.-en_0_0.pdf (Preamble) and
eupl1.1.-licence-en_0.pdf (detailed description of the licence).

This source code is provided as is, with no warranty of any kind. The user takes full responsiblity for any damage to system, and
for any scientific conclusion drawn.


### Contact
The primary contacts for the OSSOS models are:  
* Jean-Marc Petit: Jean-Marc.Petit@normalesup.org for the Fortran models (Classical, Detached, Inner, Plutinos, Scatteing, Twotinos)
* JJ Kavelaars: JJ.Kavelaars@nrc-cnrc.gc.ca, and Lowell Peltier: lowell.peltier@gmail.com for the Python models (distant resonances)


### Acknowledgement
  

If you make use of the Survey characterizations and detections please cite
the appropriate survey paper:  
* CFEPS:
   * _Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)_
* OSSOS:
   * _Bannister et al (2016) AJ, 152, 70_  
   * _Bannister et al (2018) ApJS, 236, 18_
* HiLat:
   * _Petit et al (2017), AJ, 153, 236_
* MA Survey:
   * _Alexandersen et al (2016), AJ, 152, 111_

---
## Overview
This package of provides programmes and data aiming at generating absolutely calibrated models of the various component of the small bodies of the Outer Solar System.

These models can serve as comparison for the outcome of Solar System formation simulations, or as input for running a Survey Simulator, such as the [OSSOS Survey Simulator](https://github.com/OSSOS/OSSOS).
