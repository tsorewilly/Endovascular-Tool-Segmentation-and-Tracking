%function ProcessData(FolderDir)
clear; clc;
warning off;
    FolderDir='E:\2020.04.21\Trial-01\';
    list = ls(FolderDir);
    Output = [];
    tracker=1;
    
    se = strel('disk',12);
    splitThresh1 = uint16(30000);
    splitThresh2 = uint16(40000);
    %scanBox = [235 725 300 300];
    scanBox = [50 50 minus(size(imread(deblank([FolderDir, list(3,:)]))), [50 50])];
    n = 2; %ceil(6.*rand(1,1)) % n = neighborhood

    for f = 1 :1: size(list,1) 
        fileName = deblank([FolderDir, list(f,:)]);

        if(length(fileName)>3)
            ext = extractAfter(deblank(fileName), length(deblank(fileName))-3);        
            if(strcmp(ext, 'tif') || strcmp(ext, 'png') || strcmp(ext, 'psd'))
                thisframe = imread(fileName);
                
                %thisframe = rgb2gray(imcrop(thisframe,scanBox));
                %thisframe = imcrop(thisframe,scanBox);
                figure(101); clf; imshow(thisframe); title('Original Image');
                
                %Calculate the standard deviation and the image mean for stretching.                                  
                Idouble = im2double(thisframe);        avg_fg = mean2(Idouble);              sigma_fg = std2(Idouble);
                low_fg = avg_fg-n*sigma_fg; if (low_fg<0); low_fg=0; end
                high_fg = avg_fg+n*sigma_fg; if (high_fg>1); high_fg=1; end
                
                %Adjust the contrast based on the standard deviation. %use medfilt2 or imguidedfilter for filtering
                ContAdjDSframe = imadjust(imguidedfilter(thisframe),[low_fg high_fg],[]); 
                figure(102);clf; imshow(ContAdjDSframe);    title('Current Adjusted Fore Image (Guided Filtered)');
                ContAdjDSframe(ContAdjDSframe>10) = max(max(ContAdjDSframe(:,1)));
                figure(103);clf; imshow(ContAdjDSframe);     title('Contrast Lesser than Thresh = 24G');
                %ContAdjDSframe(ContAdjDSframe>=29000) = max(max(ContAdjDSframe(:,1)));
                %figure(104);clf; imshow(ContAdjDSframe);     title('Contrast Greater than Thresh = 29G');
                                
                
                if tracker > 1
                    previmage = imcrop(imread(deblank([FolderDir, list(f-1,:)])), scanBox);
%                     figure(105); clf; imshow(previmage);        title('Previous Image');
                    Idouble = im2double(previmage);        avg_bg = mean2(Idouble);              sigma_bg = std2(Idouble);                
                    low_bg = avg_bg-n*sigma_bg; if (low_bg<0); low_bg=0; end
                    high_bg = avg_bg+n*sigma_bg; if (high_bg>1); high_bg=1; end
                    
                        
                    ContAdjBGframe = imadjust(imguidedfilter(previmage),[low_bg high_bg],[]);  
%                     figure(106);clf; imshow(ContAdjBGframe);    title('Previous Image: Fore Image (Guided Filtered)');
                    %Adjust the contrast based on the standard deviation. %use medfilt2 or imguidedfilter for filtering
                    ContAdjBGframe(ContAdjBGframe<=25000)=max(max(ContAdjBGframe));
%                     figure(107);clf; imshow(ContAdjBGframe);     title('Contrast Lesser than Thresh = 24G');
                    ContAdjBGframe(ContAdjBGframe>=29000)=max(max(ContAdjBGframe));
%                     figure(108);clf; imshow(ContAdjBGframe);     title('Contrast Greater than Thresh = 29G');

                    TrackBox = [150 625 100 250];
                    movingObject = imcrop(ContAdjDSframe-ContAdjBGframe, TrackBox);
                    figure(105);clf; imshow(movingObject);     title('Guidewire Tip');
                    
                    cannImageLines = edge(movingObject,'canny');
%                     figure(110);clf; imshow(cannImageLines);    title('Prewitt ImageLines in Image');
                    [ImgLines,theta,rho] = hough(cannImageLines);
                    Peaks = houghpeaks(ImgLines, n, 'threshold', ceil(0.3*max(ImgLines(:))));
                    lines = houghlines(cannImageLines, theta, rho, Peaks,'FillGap',5,'MinLength',7); 
                    
                    [row,col] = find(L == 2)
                
                    figure(102); hold on;
                    max_len = 0;
                    for k = 1:length(lines)
                        xy = [lines(k).point1; lines(k).point2] + [TrackBox(1,1:2);TrackBox(1,1:2)];
%                        xy = [lines(k).point1; lines(k).point2] + reshape(TrackBox, 2,2)';
                        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

                        % Plot beginnings and ends of lines
                        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
                        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','green');

                        % Determine the endpoints of the longest line segment
                        len = norm(lines(k).point1 - lines(k).point2);
                        if ( len > max_len)
                            max_len = len;
                            xy_long = xy;
                        end
                    end
                    % highlight the longest line segment
                    plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');                    
                    dim = size(movingObject); counter=1;
%                     tic
%                     for i = 1: dim(1)
%                         %for ij = 1: dim(1)
%                            map = (kolormap(movingObject(i, :), max(max(movingObject))));
%                          %  counter=counter+1;
%                         %end
%                     end
%                     toc
% 
                 end
                 tracker = tracker + 1;
            end
        end
    end
%end