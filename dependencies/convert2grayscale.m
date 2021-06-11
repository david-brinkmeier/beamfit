function output_image = convert2grayscale(input_image)
% converts input rgb/grayscale uint image to double grayscale

if ndims(input_image) == 3 %if image is RGB convert to grayscale
    output_image = double(rgb2gray(input_image));
else
    output_image = double(input_image);
end

end

