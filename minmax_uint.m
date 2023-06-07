function disp_pixel_val = minmax_unit(im)
    %https://blog.csdn.net/songzitea/article/details/8505469
    %In the display, there are usually only 8 bits (0~255), while the data has 12~16 bits.
    %If data is converted btw dynamic min and max range) to 8-bit 0~255 the process is lossy, and the resulting image often highlights some noise. 
    %We put forward some three requirements and a window-leveling method based on requirements to solve this problem; that it:
        % (1) : Make full use of the effective range of display between 0-255; 
        % (2) : Minimize the loss caused by value range compression; 
        % (3) : Don't lose the organizational part that should be prominent.

    %The algorithm: Proposed window-leveling normalization
    
    nNumPixels = size(im, 1) * size(im, 2);
    window_center
    window_width
    min_im = (2*window_center - window_width)/2.0 + 0.5;
    max_im = (2*window_center + window_width)/2.0 + 0.5;
    dFactor = 255.0/(max_im - min_im);

    for i = 1:nNumPixels % != nNumPixels; i++){
        disp_pixel_val = (pixel_val - min_im)*255.0/(max_im - min_im);

        %Non-linear transformation
        disp_pixel_val = 255.0 * pow(disp_pixel_val/(max_im-min_im), 1.0/gamma);

        %A few more points; check that values ??outside effective gray scale are process carefully
        nPixelVal = ((pixel_val - min_im)*dFactor);
        if (nPixelVal < 0) 
            disp_pixel_val = 0;
        elseif (nPixelVal > 255) 
            disp_pixel_val = 255;
        else
            disp_pixel_val = nPixelVal;
        end
    end
    disp_pixel_val = reshape(size(im, 1), size(im, 2));
    %Window-level: processing of raw data outside the min and max of a window
    % dFactor = 255.0/(double)(max - min);
    % for (i = 0; i < nNumPixels; i++){
    %     if (pixel_val < min){
    %         disp_pixel_val = 0;
    %         continue;
    %     }
    %     if (pixel_val > max){
    %         disp_pixel_val = 255;
    %         continue;
    %     }
    % 
    %     nPixelVal = (int)((pixel_val - min)*dFactor);
    % 
    %     if (nPixelVal < 0) disp_pixel_log = 0;
    %     else if (nPixelVal > 255) disp_pixel_log = 255;
    %     else disp_pixel_log = nPixelVal;
    % }
end