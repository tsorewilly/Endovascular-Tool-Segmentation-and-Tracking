function xImage = readrawfile(fname)

    fid=fopen(fname, 'r');
    Image=fread(fid,inf,'uint16','l');
    try
        xImage=reshape(Image,[1560 1440]);
        xImage=flip(imrotate(xImage,-90),2);
    catch e
        fclose(fid);
        fid=fopen(fname, 'r');
        Image=fread(fid,inf,'uint16','n');
        %disp(['RAW File is read as 512 x 256']);
        xImage=reshape(Image,[512 256]);
        xImage=xImage(1:256,:);
    end    
    fclose(fid);
    figure;imshow(xImage,[]);
end

