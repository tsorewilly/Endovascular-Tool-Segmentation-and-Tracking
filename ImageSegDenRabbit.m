%function ImageSegDenRabbit(DayDir)
    warning off;
    DayDir = 'E:\2020.05.22 (Complete Experiment)\Omisore-01\Omisore-14.51-Push-_02';
    TrialLists=dir(fullfile(DayDir));    
    FldList = {TrialLists.name};
    AnInfo = readtable('G:\Documents\My Postdoc Work\Intelligent Robot Navigation\Segmentation\imageAnnotation.xlsb');
    for fld=3:length(FldList)
        %Load Image Frames        
        SubjectDir = dir(fullfile([DayDir, '\', FldList{fld}]));
        SubjectDir = {SubjectDir.name};
        FolderDir = [DayDir, '\', FldList{fld}, '\', SubjectDir{find(not(cellfun('isempty',strfind(SubjectDir, 'Pull'))))}];
        Lists=dir(fullfile(FolderDir, '*.raw'));
        imageLists = {Lists.name};
        img_NUM = length(imageLists);
        LagLead=1; incr = 3; %While incr of 11 is 1 Sec Image Data Frames, LagLead can be optimized using rawViewer.exe

        %Query annotationTabe for trial bounds
        FldDtls = split(FolderDir, '\');
        AnnId = find(strcmp(AnInfo.FolderTag, FldDtls{end-2}) .* strcmp(AnInfo.Subject, FldDtls{end-1}) .* strcmp(AnInfo.ActionTag, FldDtls{end})==1);

        if strfind(lower(FolderDir), 'push')
            acType = 1; BgFrame = LagLead; %endId = img_NUM-LagLead; 
        elseif strfind(lower(FolderDir), 'pull')    
            acType = 2; BgFrame = img_NUM-LagLead; %endId = LagLead; incr = -20;
        else
            warndlg('Selected Folder May not Containt Push/Pull Actions');
            actionType = questdlg('Which Catheterization Action Have You Chosen', 'Action Type', 'Push', 'Pull','');
            switch actionType
                case 'Push'
                    acType = 1; BgFrame = LagLead; %endId = img_NUM-LagLead; 
                case 'Pull'
                    acType = 2; BgFrame = img_NUM-LagLead; %endId = LagLead; incr = -20;
            end
        end
         %%
        

        se = strel('diamond',7); %size(se.Neighborhood)
        se_wire=strel('line',17,100); 

        n = 1; %ceil(6.*rand(1,1)) % n = neighborhood

        Options.FrangiScaleRange=[7 7];

        %Define the background   
        bgFrame = readrawfile(fullfile(FolderDir, imageLists{BgFrame}));
        %figure(101); clf; imshow(bgFrame,[]); title('Background Frame');

        % hWaitBar = waitbar(0, 'Reading the CT Images');
        tracker=1;
        for f = LagLead : incr : img_NUM-LagLead 

            startTime = tic;
            dsFrame = readrawfile(fullfile(FolderDir, imageLists{f}));
            %figure(102); clf; imshow(dsFrame,[]); title('Current Frame Image');

            movingObject = dsFrame - bgFrame;
            movingObject=imtophat(movingObject, se);%[0 1 0; 0 1 0; 1 1 1; 0 1 0; 0 1 0]);        
            movingObject(movingObject<=-40000)=1;   movingObject(movingObject>=40000)=1;  
            %figure(103);clf; imshow((movingObject), []);     title('Moving Object');         

            %Calculate mean and std dev of the image for intensity adjustment.                                  
            Idouble = im2double(movingObject);        avg_fg = mean2(Idouble);              sigma_fg = std2(Idouble);
            low_fg = abs(avg_fg-n*sigma_fg);          high_fg = abs(avg_fg+n*sigma_fg);         
            %if (low_fg<0); low_fg=0; end    if (high_fg>1); high_fg=1; end
            lh=MinMax([low_fg high_fg], max(max(movingObject)), min(min(movingObject)));%Normalize with minmax method.

            %Adjust the contrast based on the standard deviation. %use medfilt2 or imguidedfilter for filtering
            ContAdj_MedFilter = imadjust(medfilt2(movingObject),[lh(1) lh(2)]);
            %figure(103);clf; imshow(ContAdj_MedFilter);    title('Current Adjusted Fore Image (Median Filter)'); 

            FrameVesselness = imadjust(FrangiFilter2D(ContAdj_MedFilter, Options), [0 4.3E-05]);
            %FrameVesselness = imadjust(FrameVesselness, [0.5 1]);
            %figure(104);clf; imagesc(FrameVesselness);  title('Vesselness Mapping after Median Filtering'); 
            %figure(101);clf; imshowpair(dsFrame,FrameVesselness,'montage');


            %FrameVesselness(FrameVesselness>0.5)=1; FrameVesselness(FrameVesselness<1)=0;     
            BBox = AnInfo(AnnId,5:end);
            Width = 100;    uRange=20;
            MSeq = BBox.StartY-BBox.MidY;
            FSeq = BBox.MidY-BBox.EndY;
            IPos = [BBox.StartX BBox.StartY];
            MPos = [BBox.MidX BBox.MidY];
            FPos = [BBox.EndX BBox.EndY];
            %BLPos = IPos-[Width 0]; BRPos = IPos+[Width 0];
            %TLPos = FPos-[Width 0]; TRPos = FPos+[Width 0];

            IM_Path = IPos + [ceil(([1:MSeq]/MSeq)*(MPos(1)-IPos(1))); ceil(([1:MSeq]/MSeq)*(MPos(2)-IPos(2)))]';
            MF_Path = MPos + [ceil(([1:FSeq]/FSeq)*(FPos(1)-MPos(1))); ceil(([1:FSeq]/FSeq)*(FPos(2)-MPos(2)))]';            
            %BLR_Path = BLPos + [ceil(([1:Width]/Width)*(BRPos(1)-BLPos(1))); ceil(([1:Width]/Width)*(BRPos(2)-BLPos(2)))]';
            %TRL_Path = TLPos + [ceil(([1:Width]/Width)*(TRPos(1)-TLPos(1))); ceil(([1:Width]/Width)*(TRPos(2)-TLPos(2)))]';
            %TrackLine = [BLR_Path; IM_Path+[Width 0]; MF_Path+[Width 0]; TRL_Path; flip(MF_Path)-[Width 0]; flip(IM_Path)-[Width 0]; BLR_Path(1,:)];

            CenterLine = [IM_Path; MF_Path];    
            BCells = [reshape(CenterLine(:,1) - [Width:-1:-Width], 1, []); reshape(CenterLine(:,2) - zeros(1,(2*Width)+1),1, [])]'; 
            DenoiseVessel = zeros(size(FrameVesselness));
            for i = 1:size(BCells,1)                
                DenoiseVessel(BCells(i,2), BCells(i,1)) = FrameVesselness(BCells(i,2), BCells(i,1));
            end
            %figure(105); imagesc(DenoiseVessel); title('Newly Denoised Frames');

            %===================================================================
%           tic
%           for Ci=1:length(CenterLine)
%               B_C_ells{Ci,1:(2*Width)+1} = CenterLine(:,1) - [Width:-1:-Width]; 
% %               for BiD = 1: Width
% %                   BCellsReshpId {Ci,BiD} = CenterLine(Ci,:) - [Width 0];
% %               end
%           end
%           toc;

            %This part needs further optimization
            innerLoopTime = tic;
            BCells1 = CenterLine(:,1) - [Width:-1:-Width];         
            BCells2 = CenterLine(:,2) - zeros(1,(2*Width)+1);
            for Ci=1:length(CenterLine)
                for BiD = 1: (2*Width)+1
                    BCellsReshpId{Ci,BiD} = [BCells2(Ci,BiD) BCells1(Ci,BiD)];
                    BCellsReshpData(Ci,BiD) = DenoiseVessel(BCells2(Ci,BiD), BCells1(Ci,BiD));
                end
            end            
            %===================================================================
            B_Size = size(BCellsReshpId);
            HorY = [floor(1:B_Size(1)/uRange:B_Size(1)) B_Size(1)]'; %Note that Index of y comes first
            HorX = [floor(1:B_Size(2)/uRange:B_Size(2)) B_Size(2)]';

            DenoiseVessel = zeros(size(FrameVesselness)); %cleanse all values in DenoisedVessel

            for Ry=1:uRange
                HorY_Range = HorY(Ry):HorY(Ry+1);
                for Rx=1:uRange                    
                    HorX_Range = HorX(Rx):HorX(Rx+1);
                    HorX_lent = length(HorX_Range);
%                     ChunkIndex{(3*Ry)-2, Rx} = BCellsReshpId(HorY_Select_Range, HorX_Select_Range);
%                     dsChunk = BCellsReshpData(HorY_Select_Range, HorX_Select_Range);
%                     ChunkIndex{(3*Ry)-1, Rx} = dsChunk;
%                     %Filter All-Zero Chunks and record the percentage as below                   
%                     %ChunkIndex{(Ry), Rx} = [num2str((nnz(~dsChunk)/(length(HorY_Select_Range)*length(HorX_Select_Range)))*100), '%']; 
%                     ChunkIndex{(Ry), Rx} = nnz(~dsChunk);     
                    ChunkIndex{Ry, Rx} = BCellsReshpId(HorY_Range, HorX_Range);
                    dsChunk = BCellsReshpData(HorY_Range, HorX_Range);
                    ChunkData{Ry, Rx} = dsChunk;
                    %Filter All-Zero Chunks and record the percentage as below                   
                    ChunkNonBGCent(Ry, Rx) = nnz(~dsChunk)/(length(HorY_Range)*length(HorX_Range));                    
                end
                pxTol = 2;  %2 is the optimal value, but filters out some tool pixels
                [BgIndex, ToolIndex] = mink(ChunkNonBGCent,pxTol,2); 

                for px=1:pxTol
                    cInd = cell2mat(ChunkIndex{Ry, ToolIndex(Ry, px)});
                    %(HorX_lent*ToolIndex(Ry, px))-HorX_lent+1 : HorX_lent*ToolIndex(Ry, px)
                    cselChunk = ChunkData{Ry, ToolIndex(Ry, px)};
                    for cInd_y = 1:size(cInd,1)
                        DenoiseVessel(cInd(cInd_y), cInd(cInd_y,[2:2:end])) = cselChunk(cInd_y,:);
                    end
                end
            end
            %figure(103); imagesc(DenoiseVessel); title('Moving Guidewire');
            %figure(104); histogram(DenoiseVessel); title('Tool Noise Distritubion');
            %toc(innerLoopTime);

            toc(startTime);  
            fEquiv = floor(MinMax(f, img_NUM, LagLead, 1, size(DenoiseVessel, 1)));
           % dYcoord = size(DenoiseVessel, 1)-fEquiv;%-100;
            dYcoord = fEquiv;%-100;
            dYc=[];
            %Ycoord = BBox.TrackY-fEquiv;
            while nnz(find(DenoiseVessel(abs(dYcoord), :)>0))>0
                if f<size(DenoiseVessel, 1)-BBox.TrackY
                    break;
                else
                    dYcoord = dYcoord-10;
                    dYc = [dYc dYcoord];
                end
            end
            if(dYcoord<1400);  dYcoord = dYcoord +10; end
            dXcoord = find(DenoiseVessel(dYcoord, :)>0);
            %checkUp(DenoiseVessel, dYcoord, dXcoord);
            dXcoord = mean(median(find(DenoiseVessel(dYcoord, :)>0)));

            if isnan(dXcoord)
                dXcoord = mean(BBox.StartX, BBox.EndX);%nnz(DenoiseVessel(dYcoord, :));
            end
            gWPos(tracker,:) = [dXcoord, dYcoord];
            hold on;
            text(dXcoord,dYcoord,['X: ',num2str(dXcoord), ', Y: ', num2str(dYcoord)]);
            disp(['X: ',num2str(dXcoord), ', Y: ', num2str(dYcoord)]);
            hold off;
            try
                connectGuideWire = imdilate(DenoiseVessel, se_wire);
                [Rdg, Riv, Edg] = detectRidges(connectGuideWire);
                figure(106);clf; imshow(connectGuideWire);

                % ---------    Plotting    --------
                [RwRdg, ClRdg]   = find(Rdg);      [RwRiv, ClRiv]   = find(Riv);      [RwEdg, ClEdg]   = find(Edg);
                    %figure(101); clf; imshow(bgFrame,[]); title('Background Frame'); hold on;
                    plot(ClRdg,RwRdg,'r.'); 
                    plot(ClRiv,RwRiv,'b.');
                    plot(ClEdg,RwEdg,'c.');
                connectGuideWire = bwareaopen(connectGuideWire,350);
                connectedGuideWire = bwmorph(connectGuideWire, 'skel', inf);
            catch err
            end
            %figure(107);clf; imshow(connectedGuideWire);
            for ddd=1

%             cannImageLines = edge(Edg,'canny');
%             figure(106);clf; imshow(cannImageLines);    title('Canny ImageLines in Image');
%             [ImgLines,theta,rho] = hough(cannImageLines);
%             Peaks = houghpeaks(ImgLines, n, 'threshold', ceil(0.3*max(ImgLines(:))));
%             lines = houghlines(cannImageLines, theta, rho, Peaks,'FillGap',5,'MinLength',7);  
%             
% 
% 
%             figure(102); hold on;
%             max_len = 0;
%             for k = 1:length(lines)
%                 xy = [lines(k).point1; lines(k).point2];
%                 plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%                 % Plot beginnings and ends of lines
%                 plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%                 plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','green');
% 
%                 % Determine the endpoints of the longest line segment
%                 len = norm(lines(k).point1 - lines(k).point2);
%                 if (len > max_len)
%                     max_len = len;
%                     xy_long = xy;
%                 end
%             end
%             % highlight the longest line segment
%             plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');                    
%             dim = size(movingObject); counter=1;
%     %                     tic
%     %                     for i = 1: dim(1)
%     %                         %for ij = 1: dim(1)
%     %                            map = (kolormap(movingObject(i, :), max(max(movingObject))));
%     %                          %  counter=counter+1;
%     %                         %end
%     %                     end
%     %                     toc
            end 
        
            if tracker > 1
                CosTheta = max(min(dot(gWPos(tracker-1,:), gWPos(tracker,:))/(norm(gWPos(tracker-1,:))*norm(gWPos(tracker,:))),1),-1);
                gWAngl(tracker) = acosd(CosTheta);                 
            end
            tracker = tracker + 1;
        end
        trialsDetails{fld,1} = FolderDir;
        trialsDetails{fld,2} = [gWPos gWAngl'];
        gWPos=[];
        vars2keep = {'fld','FldList','AnInfo','TrialLists','DayDir','trialsDetails'}; 
        clearvars('-except',vars2keep{:});
        fileSaveDir  = split(DayDir, '\');
        save(strcat(fileSaveDir(end), "-2.mat"), 'trialsDetails')
    end
    warning on;    
    
%end