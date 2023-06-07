clear; clc;
warning off;
DataDir = 'G:\Documents\My Postdoc Work\Publications\IJPRAI\';
AnInfo = readtable([DataDir, 'DataLinks.xlsx']);

for AnnId = 1:size(AnInfo,1)
    %Load Image Frames
    FolderDir ="E:";
    dsTrial = AnInfo(AnnId, :);
    FolderDir = join([FolderDir, string(dsTrial.FolderTag), string(dsTrial.Subject), string(dsTrial.ActionTag)],'\');
    
    Lists=dir(fullfile(FolderDir, '*.raw'));    
    imageLists = {Lists.name};
    img_NUM = length(imageLists);
    LagLead=5; incr = 50; %While incr of 11 is 1 Sec Image Data Frames, LagLead can be optimized using rawViewer.exe

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
    tracker=1; oldXcoord=1; oldYcoord=1;

    se = strel('diamond',7); %size(se.Neighborhood)
    se_wire=strel('line',17,100); 

    n = 1; %ceil(6.*rand(1,1)) % n = neighborhood

    Options.FrangiScaleRange=[7 7];

    %Define the background   
    bgFrame = readrawfile(fullfile(FolderDir, imageLists{BgFrame}));
    %figure(101); clf; imshow(bgFrame,[]); title('Guidewire Tracking Frames');
    Temp=zeros(size(bgFrame));
    % hWaitBar = waitbar(0, 'Reading the CT Images');

        for f = LagLead : incr : img_NUM-LagLead
            startTime = tic;
            dsFrame = readrawfile(fullfile(FolderDir, imageLists{f}));
            %figure(101); clf; imshow(dsFrame,[]); title('Guidewire Navigation Frames');
            %imwrite(dsFrame,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col1.png'])

            movingObject = dsFrame - bgFrame;
            movingObject=imtophat(movingObject, se);%[0 1 0; 0 1 0; 1 1 1; 0 1 0; 0 1 0]);        
            movingObject(movingObject<=-40000)=1;   movingObject(movingObject>=40000)=1;  

            %Calculate mean and std dev of the image for intensity adjustment.                                  
            Idouble = im2double(movingObject);        avg_fg = mean2(Idouble);              sigma_fg = std2(Idouble);
            low_fg = abs(avg_fg-n*sigma_fg);          high_fg = abs(avg_fg+n*sigma_fg);         
            %if (low_fg<0); low_fg=0; end    if (high_fg>1); high_fg=1; end
            lh=MinMax([low_fg high_fg], max(max(movingObject)), min(min(movingObject)));%Normalize with minmax method.

            %Adjust the contrast based on the standard deviation. %use medfilt2 or imguidedfilter for filtering
            ContAdj_MedFilter = imadjust(medfilt2(movingObject),[lh(1) lh(2)]);
            %figure(103);clf; imshow(ContAdj_MedFilter);    title('Median Filtering Result'); 
            %imwrite(ContAdj_MedFilter,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col2.png'])
            
            FrameVesselness = imadjust(FrangiFilter2D(ContAdj_MedFilter, Options), [0 4.3E-05]);
            %FrameVesselness = imadjust(FrameVesselness, [0.5 1]);
            %figure(101);clf; imagesc(FrameVesselness);  %title('Vesselness Measure Result');
            figure(101);clf; imshow(FrameVesselness);  %title('Vesselness Measure Result');
            %imwrite(FrameVesselness,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col3.png'])
            saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col2.png'])
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

            CenterLine = [IM_Path; MF_Path];    
            BCells = [reshape(CenterLine(:,1) - [Width:-1:-Width], 1, []); reshape(CenterLine(:,2) - zeros(1,(2*Width)+1),1, [])]'; 
            DenoiseVessel = zeros(size(FrameVesselness));
            for i = 1:size(BCells,1)                
                DenoiseVessel(BCells(i,2), BCells(i,1)) = FrameVesselness(BCells(i,2), BCells(i,1));
            end
            %figure(105); imagesc(DenoiseVessel); title('Guidewire Extraction');
            %imwrite(DenoiseVessel,[DataDir, 'Evaluation\Validation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col3.png'])
            %===================================================================
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
            %figure(104); histogram(DenoiseVessel); title('Tool Noise Distritubion');
            %toc(innerLoopTime);
            SegTime=toc(startTime);
            TrackStart=tic;
            %toc(startTime);  
            fEquiv = floor(MinMax(f, img_NUM, LagLead, 1, size(DenoiseVessel, 1)));
            dYcoord = size(DenoiseVessel, 1)-fEquiv;%-100;
            dYc=[];
            %Ycoord = BBox.TrackY-fEquiv;
            while nnz(find(DenoiseVessel(dYcoord, :)>0))>0
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
            %hold on;
            %text(dXcoord,dYcoord,['X: ',num2str(dXcoord), ', Y: ', num2str(dYcoord)]);
            %disp(['X: ',num2str(dXcoord), ', Y: ', num2str(dYcoord)]);
            %hold off;                
    %             figure(106);clf; imshow(connectGuideWire);
            figure(101); imshow(DenoiseVessel); %title('Guidewire Tracking');
            %imwrite(DenoiseVessel,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col4.png'])
            
            % ---------    Guidewire Tracking, (1) Whole Body, (2) Tip    --------
            saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col3.png'])
            connectGuideWire = imdilate(DenoiseVessel, se_wire);
            %connectGuideWire = imdilate(FrameVesselness, se_wire);
            [Rdg, Riv, Edg] = detectRidges(connectGuideWire);
            % ---------    Plotting    --------
            [RwRdg, ClRdg]   = find(Rdg);      [RwRiv, ClRiv]   = find(Riv);      [RwEdg, ClEdg]   = find(Edg);
%             figure(101); clf; imshow(dsFrame,[]); 
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col1.png'])
%             
            %figure(101); clf; imshow(dsFrame,[]); hold on;
                %plot(ClRdg,RwRdg,'r.');                 plot(ClRiv,RwRiv,'b.');
                %plot(ClEdg,RwEdg,'c.');
%                 text(10,50,['Seg Time: ', num2str(SegTime), ' Secs'], 'Color','white','FontSize',14);
%                 if tracker > 1  %Change to if guidewire detected
%                     CosTheta = max(min(dot(gWPos(tracker-1,:), gWPos(tracker,:))/(norm(gWPos(tracker-1,:))*norm(gWPos(tracker,:))),1),-1);
%                     gWAngl(tracker) = acosd(CosTheta);             
% 
%                     text(10,150,['Distance:  ', num2str(distBTW(gWPos(tracker-1,:), gWPos(tracker,:))), ' px'], 'Color','white','FontSize',14);
%                     text(10,250,['Direction: ', num2str(gWAngl(tracker)), '^o'], 'Color','white','FontSize',14);
%                 else
%                     text(10,150,'Distance:  0 px', 'Color','white','FontSize',14);
%                     text(10,250,'Direction: 0^o', 'Color','white','FontSize',14);
%                 end
                tracker = tracker + 1;
                
                %Frame = getframe(gca);            FrameData = Frame.cdata;
                %figure; imshow(FrameData)
                %imwrite(FrameData,[DataDir, 'Evaluation\Seq', '.png'],'PNG');
                %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col5.png'])
                %imwrite(FrameData,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col4.png'],'PNG')
                %trackedimage = imread([DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col1.png']);
                %labeledImage = bwlabel(DenoiseVessel,8);

                %Copy aniogram for validation
                %mkdir(join([DataDir, extractBefore(string(dsTrial.FolderTag),11), string(dsTrial.Subject), string(dsTrial.ActionTag)],'\'));
                %copyfile(fullfile(FolderDir, imageLists{BgFrame}),join([DataDir, extractBefore(string(dsTrial.FolderTag),11), string(dsTrial.Subject), string(dsTrial.ActionTag), imageLists{BgFrame}],'\'));
        end

    %repDetails = [gWPos gWAngl'];
    msgbox(['Trial ', num2str(AnnId), ' Operation Completed'],'Success', 'modal');%,'custom',icondata,summer);
end