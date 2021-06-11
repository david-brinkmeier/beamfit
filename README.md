# Beamfit

![Alt text](docs/logo.png?raw=true "Title")


## Offline and Online analysis for laser beams, models include
  - Least square elliptical rotated Gaussian fit (D4σ specification)
    - cf. beamfit/docs/
  - Least square trepanning symmetrical Gaussian fit (D4σ specification)
  - Second-order moments / ISO11146 beam width for elliptic beams (D4σ specification)
    - Based on this [implementation][imagemoments] by Raphaël Candelier with minor adjustments to fit the ISO11146 spec.
  
## Image processing
  - Built-in functionality includes various DC-offset removal / noise removal techniques including 
  [noisecomp][kovesi] from Peter Kovesi and [TV-L1 denoising][tvl1] from Manolis Lourakis.
  - Automated and/or GUI-cropping / pre-scaling of input for faster fitting etc.
  - Why image cleanup? Second-order moment determination of beam diameter/radius is particularly prone to measurement noise
    - Warning: Due to the definition of beam diameter through the second order moments this means it's also possible to effectively falsify the ISO11146 measurement through usage of excessive denoising.
  - Additionally / alternatively an offset background image which is subtracted from the input may be provided.
  
- Offline analysis supports image and video processing
  - Examples and results for usage with example_general_purpose.m are provided in \examples\
  - example_synthetic.m allows for generation and analysis of noisy rotated offset elliptical Gaussians as well as noisy trepanning offset symmetrical Gaussians for testing purposes.
- Online analysis supports IDS uEye cameras through uEyeDotNet.dll
  - Specification of uEye cameras (Pixel pitch etc.) are auto-detected.
  - Minor bugfixes were applied to the [uEye-dotnet Matlab library][ueye_lib] from [Dr. Adam Wyatt][adamwyatt]


## Features

- Horrible mess of procedural code, held in place by duct-tape.
  If I were to do it today, I would do it properly and most likely in Python.
- That being said, the code works and has proven to be useful in practice. Especially when compared to many proprietary / expensive beam measurement software/hardware bundles.
- Usage is intended with Camera sensors WITHOUT lens. If you have an imaging system, you need to take the specification of your imaging setup into account.
If this means nothing to you, a commercial system may be better suited for your purposes.

## Disclaimer

- As stated in the license...code provided "as is". If you need certified measurements, use certified measurement tools and software!

[imagemoments]: <http://raphael.candelier.fr/?blog=Image%20Moments>
[kovesi]: <https://www.peterkovesi.com/matlabfns/>
[tvl1]: <https://de.mathworks.com/matlabcentral/fileexchange/57604-tv-l1-image-denoising-algorithm/>
[ueye_lib]: <http://matlabtidbits.blogspot.com/2016/12/ueye-camera-interface-in-matlab-net.html>
[adamwyatt]: <https://www.clf.stfc.ac.uk/Pages/Adam-Wyatt.aspx>
