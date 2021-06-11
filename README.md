# Beamfit

- Offline and Online analysis for laser beams, models include
  - least square elliptical rotated Gaussian fit (D4σ specification)
  - least square trepanning symmetrical Gaussian fit (D4σ specification)
  - second-order moments / ISO11146 beam width for elliptic beams (D4σ specification)
  
- Image processing
  - built-in functionality includes various DC-offset removal / noise removal techniques including 
  [noisecomp][kovesi] from Peter Kovesi and [TV-L1 denoising][tvl1] from Manolis Lourakis.
  - Why is this important? Second-order moment determinateion of beam diameter/radius is particularly prone to measurement noise
  - Additionally / alternatively an offset background image may also be provided.
  
- Offline analysis supports image and video processing
  - examples and results for usage with example_general_purpose.m are provided in \examples\
- Online analysis supports IDS uEye cameras through uEyeDotNet.dll


## Features

- Horrible mess of procedural code, held in place by duct-tape.
  If I were to do it today, I would do it properly and most likely in Python.
- That being said, the code works and has proven to be useful in practice. Especially when compared to many proprietary / expensive beam measurement software/hardware bundles.
- Usage is intended with Camera sensors WITHOUT lens. If you have an imaging system, you need to take the specification of your imaging setup into account.
If this means nothing to you, a commercial system may be better suited for your purposes.


[kovesi]: <https://www.peterkovesi.com/matlabfns/>
[tvl1]: <https://de.mathworks.com/matlabcentral/fileexchange/57604-tv-l1-image-denoising-algorithm/>
