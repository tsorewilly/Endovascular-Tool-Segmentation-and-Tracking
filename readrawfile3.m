function [img,K] = readrawfile3(image_file,imagesize)
    %fid = fopen(image_file,'r','b');
    fid = fopen(image_file,'r');
    clear img
    img = fread(fid,'uint8');
    fclose(fid);
    K = length(img)/(imagesize(1)*imagesize(2));
    img = reshape(img,imagesize(1),imagesize(2),K);
end