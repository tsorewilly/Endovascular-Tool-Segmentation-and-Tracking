clear; clc;
warning off;
FolderDir = 'G:\Documents\My Postdoc Work\Publications\IJPRAI\Evaluation\';
topN = 20;
NosieThresh = 20
scanBox = [221 753 1306 1206];

    %Load Ground Truth Image (GTI) 
    List_GTI=dir(fullfile(FolderDir, '*-Col1.png'));    
    imageGTI = {List_GTI.name};
    img_NUM = length(imageGTI);
    
    %Load Segmented Image Frames (SIF)
    List_SIF=dir(fullfile(FolderDir, '*-Col2.png'));    
    imageSIF = {List_SIF.name};
    img_NUM = length(imageSIF);
    fprintf('Id \t\t\t\tTracAcc\t\t\t\t\tTrackConn\t\t\tTrackArea \n')%

    curSeq = 'Seq1-'; cS_id=1; startId=1;
    [d_deta_p, d_deta_r, d_average, d_std, rIndex, tIndex] = deal(0);
    
    for f = 1:img_NUM        
        MthdEval{f,1} = imcrop(imread(fullfile(FolderDir, imageGTI{f})), scanBox);
        scanObject=MthdEval{f,1}; scanObject(scanObject>0)=1; MthdEval{f,1} = scanObject;

        MthdEval{f,2} = imcrop(rgb2gray(imread(fullfile(FolderDir, imageSIF{f}))), scanBox);
        scanObject=MthdEval{f,2}; scanObject(scanObject>0)=1; MthdEval{f,2} = scanObject;

        fromIndex = size(MthdEval{1},1); toIndex = 1;
        for row = fromIndex:-1:toIndex
            GW_PX_IN_GTI = find((MthdEval{f,1}(row,:)));
            GTI_r = nnz(MthdEval{f,1}(row, GW_PX_IN_GTI));
            SIF_r = nnz(MthdEval{f,2}(row, GW_PX_IN_GTI));
            if GTI_r == 0
                GW_Conn(row,1) = 0;
            else
                GTI_y2 = row;
                if startId > topN
                    GTI_y1 = row+topN;
                else
                    GTI_y1 = fromIndex;
                end
                GW_Conn(row,1) = SIF_r/GTI_r;
                GW_Area(row,1) = SIF_r;
                GT_Area(row,1) = GTI_r;
            end

            if SIF_r ~= 0
                SIF_y2 = row;
                if startId > topN
                    SIF_y1 = row+topN;
                else
                    SIF_y1 = fromIndex;
                end                
            end
            %disp([row GTI_y1 GTI_y2 SIF_y1 SIF_y2]);
            startId = startId + 1;
        end
        startId = 1;
        MthdEval{f,3} = GW_Conn;
        
        %Determine the tracking displacement and orientation changes
        GTI_x1 = ceil(median(find((MthdEval{f,1}(GTI_y1,:)))));
        GTI_x2 = ceil(median(find((MthdEval{f,1}(GTI_y2,:)))));
        
        SIF_x1 = ceil(median(find((MthdEval{f,2}(SIF_y1,:)))));
        SIF_x2 = ceil(median(find((MthdEval{f,2}(SIF_y2,:)))));
        
        %disp([GTI_x1 GTI_x2 GTI_y1 GTI_y2; SIF_x1 SIF_x2 SIF_y1 SIF_y2]);
        if(~isnan(GTI_x1) & ~isnan(GTI_x2) & ~isnan(SIF_x1) & ~isnan(SIF_x2))
            GTI_P = distBTW([GTI_x1 GTI_y1],[GTI_x2 GTI_y2]);
            GTI_R = acosd(max(min(dot([GTI_x1 GTI_y1], [GTI_x2 GTI_y2])/(norm([GTI_x1 GTI_y1])*norm([GTI_x2 GTI_y2])),1),-1));

            SIF_P = distBTW([SIF_x1 SIF_y1],[SIF_x2 SIF_y2]);
            SIF_R = acosd(max(min(dot([SIF_x1 SIF_y1], [SIF_x2 SIF_y2])/(norm([SIF_x1 SIF_y1])*norm([SIF_x2 SIF_y2])),1),-1));

            delta_P = abs(GTI_P-SIF_P);
            delta_R = abs(GTI_R-SIF_R);
                        
            tIndex = tIndex + 1;
            %disp([GTI_x1 GTI_x2 GTI_y1 GTI_y2 SIF_x1 SIF_x2 SIF_y1 SIF_y2]); disp([delta_P delta_R]);
        else
            [delta_P, delta_R] = deal(0);
        end

        % Calculate and store the mean ± std of the non-zero elements. 
        indexToNonZero = GW_Conn~=0;
        aver = mean(GW_Conn(indexToNonZero));
        stdev = std(GW_Conn(indexToNonZero));
        area = sum(GW_Area(:,1))/sum(GT_Area(:,1));
        
        
        
        if(isnan(aver))
            [aver, stdev] = deal(0);
        else
            rIndex=rIndex+1;
        end
        
        %Filter out cases affected by noise
        if(delta_P>NosieThresh)
            delta_P = NosieThresh;
            delta_R = 0;
        end

        MthdEval{f,4} = [delta_P, delta_R];
        MthdEval{f,5} = [aver, stdev];
        MthdEval{f,6} = area;
        
        d_deta_p = d_deta_p + delta_P;
        d_deta_r = d_deta_r + delta_R;
        d_average = d_average + aver;
        d_std = d_std + stdev;        
        
        if(~strcmp(extractBefore(imageGTI{f}, 6), curSeq))
            fprintf('%d\t\t\t%.4f ± %.4f\t\t\t%.4f ± %.4f\t\t\t%.4f\n', cS_id, d_deta_p/tIndex, d_deta_r/tIndex, (d_average/rIndex)*100, (stdev/rIndex)*100, area*100)
            [d_deta_p, d_deta_r, d_average, d_std, rIndex] = deal(0);
            curSeq = extractBefore(imageGTI{f}, 6);
            cS_id=cS_id+1;
        end
    end    
fprintf('%d\t\t\t%.4f ± %.4f\t\t\t%.4f ± %.4f\t\t\t%.4f\n', cS_id, d_deta_p/tIndex, d_deta_r/tIndex, (d_average/rIndex)*100, (stdev/rIndex)*100, area*100)