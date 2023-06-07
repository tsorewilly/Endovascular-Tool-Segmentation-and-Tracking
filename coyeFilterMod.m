function outImage = coyeFilterMod(imInput)
    % Convert RGB to Gray via PCA
    B = imInput;
    lab = rgb2lab(imInput);
    fnVr = 0;
    wlab = reshape(bsxfun(@times,cat(3, 1-fnVr, fnVr/2, fnVr/2), lab),[],3);
    [C,S] = pca(wlab);
    S = reshape(S,size(lab));
    S = S(:,:,1);
    gray = (S-min(S(:)))./(max(S(:))-min(S(:)));
    
    J = adapthisteq(gray,'numTiles',[8 8],'nBins',128); % Contrast Enhancment of gray image using CLAHE
    
    h = fspecial('average', [9 9]);                     % Background Exclusion, Apply Average Filter
    JF = imfilter(J, h);
    Z = imsubtract(JF, J);                              % Take the difference between the gray image and Average Filter                
    %level=isodata(Z);                                   % Threshold using the IsoData Method
    level = graythresh(Z);                             % isodata was proposed in the study for thresholding

    BW = im2bw(Z, level);                               % Convert to Binary            
    BW = imbinarize(Z);
    BW2 = bwareaopen(BW, 100);                          % Remove small pixels
    outImage = imoverlay(B, imcomplement(BW2),[0 0 0]);   % Overlay and output
end