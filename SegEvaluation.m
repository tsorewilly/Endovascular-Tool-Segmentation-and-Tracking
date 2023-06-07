clear; clc;
warning off;
DataDir = 'G:\Documents\My Postdoc Work\Publications\QMIS\';
AnInfo = readtable([DataDir, 'DataLinks.xlsx']);

for AnnId = 13%size(AnInfo,1)
    %Load Image Frames
    
    FolderDir ="E:";
    dsTrial = AnInfo(AnnId, :);
    FolderDir = join([FolderDir, string(dsTrial.FolderTag), string(dsTrial.Subject), string(dsTrial.ActionTag)],'\');
    
    Lists=dir(fullfile(FolderDir, '*.raw'));    
    imageLists = {Lists.name};
    img_NUM = length(imageLists);
    LagLead=5; incr = 10; %While incr of 11 is 1 Sec Image Data Frames, LagLead can be optimized using rawViewer.exe

    %Query annotationTabe for trial bounds
    FldDtls = split(FolderDir, '\');
    AnnId = find(strcmp(AnInfo.FolderTag, FldDtls{end-2}) .* strcmp(AnInfo.Subject, FldDtls{end-1}) .* strcmp(AnInfo.ActionTag, FldDtls{end})==1);

    if strfind(lower(FolderDir), 'push')
        acType = 1; BgFrame = LagLead; %endId = img_NUM-LagLead; 
    elseif strfind(lower(FolderDir), 'pull')    
        acType = 2; BgFrame = img_NUM-LagLead; %endId = LagLead; incr = -20;
    else
        warndlg('Selected Folder May not Contain Push/Pull Actions');
        actionType = questdlg('Which Catheterization Action Have You Chosen', 'Action Type', 'Push', 'Pull','');
        switch actionType
            case 'Push'
                acType = 1; BgFrame = LagLead; %endId = img_NUM-LagLead; 
            case 'Pull'
                acType = 2; BgFrame = img_NUM-LagLead; %endId = LagLead; incr = -20;
        end
    end

    %
    tracker=1; oldXcoord=1; oldYcoord=1;

    se = strel('diamond',7); %size(se.Neighborhood)
    se_wire=strel('line',17,100); 

    n = 1; %ceil(6.*rand(1,1)) % n = neighborhood

    Options.FrangiScaleRange=[7 7];

    %Define the background   
    bgFrame = readrawfile(fullfile(FolderDir, imageLists{BgFrame}));
    figure(101); clf; imshow(bgFrame,[]); title('Guidewire Tracking Frames');
    Temp=zeros(size(bgFrame));
    hWaitBar = waitbar(0, 'Reading the CT Images');

        for f = LagLead : incr : img_NUM-LagLead
%             if(f==105 | f==205 | f==305 | f==405 | f==505 | f==605 | f==705)
%                 continue
%             end
            startTime = tic;
             dsFrame = readrawfile(fullfile(FolderDir, imageLists{f}));
             imwrite(dsFrame,[DataDir, 'Evaluation\temp.tiff'])
            figure(101); clf; imshow(dsFrame,[]); %title('Guidewire Navigation Frames');
            saveas(gcf,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col1.png'])
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
            %figure(1);clf; imshow(FrameVesselness);  title('Vesselness Measure Result');
            %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col2.png'])

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
            figure(101); imshow(DenoiseVessel); %title('Proposed Segmentation'); colormap(gray); axis off; axis image;
            saveas(gcf,[DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col2.png'])
            
            % ---------    Guidewire Tracking, (1) Whole Body, (2) Tip    --------
            %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col3.png'])
            connectGuideWire = imdilate(DenoiseVessel, se_wire);
            %connectGuideWire = imdilate(FrameVesselness, se_wire);
            [Rdg, Riv, Edg] = detectRidges(connectGuideWire);
            % ---------    Plotting    --------
            [RwRdg, ClRdg]   = find(Rdg);      [RwRiv, ClRiv]   = find(Riv);      [RwEdg, ClEdg]   = find(Edg);
%             figure(101); clf; imshow(dsFrame,[]); 
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col1.png'])
%             
            figure(101); clf; imshow(dsFrame,[]); hold on;
            %plot(ClRdg,RwRdg,'r.');                 plot(ClRiv,RwRiv,'b.');
            plot(ClEdg,RwEdg,'c.');
            saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col3.png'])
            
            tracker = tracker + 1;
          
            %dsFrame = imresize(dsFrame, [512 (512*(size(dsFrame,2)/size(dsFrame,1)))]);
%            Frame_res=(MinMax(dsFrame, max(max(dsFrame)), min(min(dsFrame)), 0, 255)); %Convert image into 0-255 and of double type
            %figure(101); imshow(Frame_res, []);  
            %saveas(gcf, [DataDir, 'Evaluation\Seq-Temp.tiff'])

%             %% -U01- Cervantes_SSG_Otsu
%             OutCervantSSG = Cervantes_SSG_Otsu(Frame_res);
%             figure(101); imshow(OutCervantSSG); %title('Cervantes Single-scale Garbon Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col3.png'])
% 
%             %% -U02- Implement Cosfire Method for Evaluation              
%             symft = struct(); % Symmetric filter params
%             symft.sigma = 4.4;  symft.len = 8;	symft.sigma0 = 30;  symft.alpha= 0.7;
%             
%             asymft = struct(); % Asymmetric filter params
%             asymft.sigma = 4.8; asymft.len = 22; asymft.sigma0 = 10;   asymft.alpha = 0.1;
%             
%             [output.respimage] = BCOSFIRE_media15(double(gray2rgb(Frame_res)), symft, asymft, 0.9);            
%             output.segmented = (output.respimage > 37);
%             OutCosfire = output.segmented;
%             figure(101); imshow(OutCosfire); %title('B-COSFIRE Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col4.png'])
%             
%             %% -U03- Implement Fast Deterministic Multi-scale Method for Evaluation  
%             addpath('EyeRetina');
%             figure(101); imshow(Frame_res, []);   saveas(gcf, [DataDir, 'Seq-Temp.tiff'])
%             %scanBox = [137 46 1217 1123];
%             scanBox = [221 753 1308 1206];
%             dsFrame_line = imresize(imcrop(imread([DataDir, 'Seq-Temp.tiff']), scanBox), [1440 1560]);
%             [GC,ATW,ATG,Vs,ATW2,VsM,dilateEdge] = FnTrackInit8(dsFrame_line);
% 
%             OutFastMScale = FnTrack21(GC,VsM,dilateEdge);
%             figure(101); imshow(OutFastMScale, []); %title('Det. Multi-scale Segmentation'); %colormap(gray); axis off; axis image; 
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col5.png'])
%             
%             %% -U04- Heneghan
%             addpath('9.VesselCharacterization');
%             OutHeneghan = Heneghan(gray2rgb(Frame_res));
%             figure(101); imshow(OutHeneghan); %title('Heneghan Segmentation');  %colormap(gray); axis off; axis image;           
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col6.png'])
%             
%             %% -N01- Heneghan-WHT1
% %             addpath('5.VesselSegmenter');
% %             %See code in 5.VesselSegmenter folder
% %             OutHeneghWHT1=wth1(Frame_res, Frame_res, 11, 11, 0.15, 180, 10);
% %             figure(6); imshow(OutHeneghWHT1); %title('Heneghan-WHT1 Segmentation');   %colormap(gray); axis off; axis image;           
% %             %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col10.png'])
%           
%             %% -U07- Coye Filter    
%             addpath('6.CoyeFilter');
%             OutCoyeFilter = coyeFilterMod(dsFrame_line);
%             figure(101); imshow(OutCoyeFilter, []); %title('Coye Filter Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col7.png'])
% 
%              %% -N02- using Principal Curvature for Vessel Extraction in retinal fundus Image             
% %             addpath('8.PrinCurv');
% %             OutPrinCurv = vesselSegPC(gray2rgb(Frame_res));
% %             figure(8); imshow(OutPrinCurv); %title('Principal Curvature Segmentation'); %colormap(gray); axis off; axis image;
% %             %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col9.png'])
% %             
%              %% -N03- Cervantes_MGMFMSGMLP_Threshold
% %             OutCervantThresh = Cervantes_MGMFMSGMLP_Threshold(Frame_res);
% %             figure(9); imshow(OutCervantThresh); %title('Cervantes Multi-scale ANN Segmentation'); %colormap(gray); axis off; axis image;
% %             %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col9.png'])
% 
%             %% -U08- Changwimaluang_GMF_LocalEntropy
%             OutChangGMF = Changwimaluang_GMF_LocalEntropy(Frame_res);
%             figure(101); imshow(OutChangGMF); %title('Chang GMF Local Entropy Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col8.png'])
% 
%             %% -U09- Kang_TH_GMF_Degree
%             OutKangTHGMF = Kang_TH_GMF_Degree(Frame_res);
%             figure(101); imshow(OutKangTHGMF); %title('Kang TH GMF Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col9.png'])
% 
%             %% -U10- Nguyen_MLD_Threshold
%             OutNguyenMLDThresh = Nguyen_MLD_Threshold(Frame_res);
%             figure(101); imshow(OutNguyenMLDThresh); %title('Nguyen MLD Threshold Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col10.png'])
% 
%             %% -U11- Qian_MTH_Background
%             OutQianMTHBkgd = Qian_MTH_Background(Frame_res);
%             figure(101); imshow(OutQianMTHBkgd); %title('Qian MTH Background Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col11.png'])
% 
%             %% -N04- Cruz_MGMF_Otsu
% %             OutCruzMGMF = Cruz_MGMF_Otsu(Frame_res);
% %             figure(12); imshow(OutCruzMGMF); %title('Cruz MGMF Otsu Segmentation'); %colormap(gray); axis off; axis image;
% %             %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col9.png'])
% 
%             %% -N05- Eiho_TH_Background
% %             OutEihoTH = Eiho_TH_Background(Frame_res);
% %             figure(101); imshow(OutEihoTH); %title('Eiho TH Bkgd Segmentation'); %colormap(gray); axis off; axis image;
% %             saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col12.png'])
% 
%             %% -N06- Kang_GMF_Degree
% %             OutKangGMF = Kang_GMF_Degree(Frame_res);
% %             figure(14); imshow(OutKangGMF,[]); %title('Kang GMF Segmentation'); %colormap(gray); axis off; axis image;
%             %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col9.png'])
% 
%             %% -N07- Kang_TH_GMF_MaxEntropy
% %             OutKangTHGMFMax = Kang_TH_GMF_MaxEntropy(Frame_res);
% %             figure(16); imshow(OutKangTHGMFMax); %title('Kang TH GMF MaxEntropy Segmentation'); %colormap(gray); axis off; axis image;
%             %saveas(gcf, [DataDir, 'Evaluation\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col9.png'])

             %% -U12- Li_Hess_VesselRepair
%             OutLiHessVesRep = Li_Hess_VesselRepair(Frame_res);
%             figure(101); imshow(OutLiHessVesRep); %title('Hessian Vessel Repair Segmentation'); %colormap(gray); axis off; axis image;
%             saveas(gcf, [DataDir, 'Li results\Seq' num2str(AnnId),'-Frm', num2str(f), '-Col12.png'])

        end

    %repDetails = [gWPos gWAngl'];
    msgbox(['Trial ', num2str(AnnId), ' Operation Completed'],'Success', 'modal');%,'custom',icondata,summer);
end