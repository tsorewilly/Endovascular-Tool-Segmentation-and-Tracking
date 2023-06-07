function disp_pixel_val = MinMaxVal(im, min, max)
    nNumPixels = size(im, 1) * size(im, 2);
    pixel_val = im(:);
    for i=1 : nNumPixels
		disp_pixel_val = (pixel_val - min)*255.0/(double)*(max - min);
    end
    disp_pixel_val = reshape(disp_pixel_val, size(im, 1),size(im, 2));
end
