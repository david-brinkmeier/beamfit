function output_image = ShotAndReadNoise(input_image,settings)
% Source: Modified from Photon_Gaussian_Noise.m
% https://www.researchgate.net/profile/Milo_Hyde
% 
%
% quantum_well_depth in photoelectrons, smaller values = more noise
% sigma_read in photoelectrons RMS
%
%
% Readout Noise and dark current: http://homepage.physics.uiowa.edu/~pkaaret/2013f_29c137/Lab03_noise.html
% http://homepage.physics.uiowa.edu/~pkaaret/2013f_29c137/Lab03_noise.html
% The CCD in the ST-402 can make accurate, but not perfect, measurements of the charge
% accumulated in each CCD pixel.  One important limitation is the noise associated with
% the electronics that amplifies and digitizes the charge signal in the CCD readout.
% Even if the exact same charge is placed in a given pixel in two different images,
% this noise will produce fluctuations in the number of ADU (analog to digital units) recorded.
% We can measure the noise by taking repeated images with the same amount of charge in each pixel.
% Since the amount of charge in a pixel depends on the amount of light entering that pixel,
% the easiest way to get the same amount of charge is to have no light enter the pixel.
% Thus, for these images, we block any light from entering the camera by simply not opening the shutter.
%
% Because the analog to digital converter (ADC) in the CCD reports only positive values,
% while the noise fluctuations can be positive or negative, a constant offset,
% called a 'bias', is added to the ADC value.  Thus, these images obtained with
% no light entering the camera and with no exposure length are called 'bias' frames.
% The CCD noise is the fluctuation in the ADU counts around the constant bias value.
% Even in the absence of light, charge will accumulate in each pixel. This is because
% the CCD is at a temperature above absolute zero and thermal fluctuations in the
% silicon can release electrons. This accumulation of charge is called 'dark current'.
% Dark current depends on the CCD temperature (as you will see in the analysis section).

% ----- Pixel Well Depth ----- %
% What is Pixel Well Depth?
% May 25, 2017 | Detectors, Glossary of Spectroscopy Terms, Uncategorized
% https://www.stellarnet.us/pixel-well-depth/
% Pixels work by converting photons into electrons.
% Most detectors store electrons in the pixels while they are exposed to the light,
% then convert the electrons into a signal after exposure. The pixel well depth,
% also called the full well capacity, determines how many electrons the pixel
% can store before it saturates and starts “leaking” electrons to
% neighboring pixels, which distorts the signal.

%%

input_image_photoelectrons = round(input_image*settings.quantum_well_depth/2^8);

% ----- Shot Noise ----- %
image_shot_noise = poissrnd(input_image_photoelectrons);
corrfact = max(abs(image_shot_noise(:)))/max(input_image(:));
image_shot_noise = image_shot_noise./corrfact;

% ----- Read Noise ----- %
% The signal-to-noise ratio is an important factor from the viewpoint of
% accurate measurements. Here, the signal-to-noise ratio is defined as the
% ratio of the mean value of the signal count rate to the fluctuations of
% the counted signal and noise pulses (expressed in standard deviation or
% root mean square)

image_read_noise = settings.sigma_read*randn(size(image_shot_noise));
read_bias = -min(image_read_noise(:));
image_read_noise = image_read_noise+read_bias;

% ----- Combined Noise Normalized ----- %
image_shot_noise_read_noise = image_shot_noise + image_read_noise;
corrfact = max(abs(image_shot_noise_read_noise(:)))/max(input_image(:));
image_shot_noise_read_noise = image_shot_noise_read_noise./corrfact;

% ----- Output Image ----- %
output_image = image_shot_noise_read_noise;

%%

if settings.plot == true
    [noise_eval,figdidexist] = genORselectfigbyname('Noise evaluation');
    if ~figdidexist
        screensize = get(0,'ScreenSize');
        set(noise_eval,'Position',[screensize(3)-800-8 46 800 640]);
    end
    
    [ha, ~] = tight_subplot(2,2,[.08 .05],[.05 .05],[.15 .1]);
    axes(ha(1))
    imagesc(input_image); colormap(gray); axis image; axis off; set(gca,'YDir','normal')
    title('Original Image');
    axes(ha(2))
    imagesc(image_read_noise), colormap(gray); axis image; axis off; set(gca,'YDir','normal')
    title('Read Noise');
    caxis([min(input_image(:)) max(input_image(:))])
    axes(ha(3))
    imagesc(image_shot_noise); colormap(gray); axis image; axis off; set(gca,'YDir','normal')
    title('Shot Noise');
    axes(ha(4))
    imagesc(image_shot_noise_read_noise); colormap(gray); axis image; axis off; set(gca,'YDir','normal')
    title('Shot Noise + Read Noise');
end

if settings.saveimages == true
    name = strcat('len_',num2str(settings.length),'_','w0_',num2str(settings.w0a),...
        '-',num2str(settings.w0b),'_x0y0_',num2str(settings.x0),'-',num2str(settings.y0),'_',...
        num2str(settings.angle),'deg_pitch_',num2str(settings.pixelpitch));
    
    imwrite(uint8(input_image),strcat(name,'.png'))
    imwrite(uint8(image_shot_noise_read_noise),strcat(name,'_noisy.png'))
end

end

