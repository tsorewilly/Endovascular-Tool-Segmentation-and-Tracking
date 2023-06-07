clear; clc;
warning off;
FolderDir = 'G:\Documents\My Postdoc Work\Publications\QMIS\Evaluation\';

MthdEvalName = {'B-GroundTruth','J-Kang','L-Qian','C-Proposed','D-Cervantes_SSG', 'E-Azzopardi','K-Nguyen', ...
    'F-Vlachos','G-Heneghan','H-CoyeFilter ','I-Chanwimaluang'};%, 'M-Li'};
  
scanBox = [221 753 1308 1206];
mLent = length(MthdEvalName);

for Seq = 1:12
    %Load Image Frames    
    Lists=dir(fullfile(FolderDir, ['Seq', num2str(Seq),'-*']));    
    imageLists = {Lists.name};
    img_NUM = length(imageLists);
    ind=1; fdex=1; MthdEval={};
    
    fprintf('Model     TPos	  FPos	TNeg      FNeg	  ACC	mACC	ICA_1	ICC_1	  SEN	  SPE	 FDR	 FOR	 DSC	 PPR	 MCC	 F1S	 AROC \n')%
    fprintf('**********************************************************************************************************************************************\n')%

    CumSum = {zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17), zeros(1,17)};        
    for f = 1:mLent:img_NUM
        for m=1:mLent
            n=m;%-1;
            if m ==1
                MthdEval{fdex,m} = imcrop(imread(fullfile(FolderDir, imageLists{ind})), scanBox);
            else
                MthdEval{fdex,m} = imcrop(rgb2gray(imread(fullfile(FolderDir, imageLists{ind}))), scanBox);
            end
           scanObject=MthdEval{fdex,m}; scanObject(scanObject>0)=1; MthdEval{fdex,m} = scanObject;
%            figure (m); imshow(MthdEval{fdex,m}, []); %title(MthdEvalName{m})  
           ind=ind+1;

           %Now compare against ground truth
           if m > 0           
               GT_GW  = MthdEval{fdex,1};
               MTD_GW = MthdEval{fdex,m};
              [TP, TN, FP, FN] = deal(0);
               for row = 1:size(GT_GW,1)
                   for col = 1:size(GT_GW,2)
                       %Analyze just guidewire pixels to compute evaluation metrics
                       if GT_GW(row,col)== 1 & MTD_GW(row,col)== 1      %Check True Positive (TP)               
                           TP = TP + 1;
                       elseif (GT_GW(row,col)==0 & MTD_GW(row,col)==1)  %Check False Positive (FP)
                           FP = FP + 1;
                       elseif GT_GW(row,col)== 0 & MTD_GW(row,col)== 0 %Check True Negative (TN)
                           TN = TN + 1;
                       elseif (GT_GW(row,col)==1 & MTD_GW(row,col)==0) %Check False Negative (TN)
                           FN = FN + 1;
                       end
%                        if GT_GW(row,col)== 1 || MTD_GW(row,col)== 1                    
%                            if (GT_GW(row,col)==1 & MTD_GW(row,col)==1) %Check True Positive (TP)
%                                TP = TP + 1;
%                            elseif (GT_GW(row,col)==0 & MTD_GW(row,col)==1) %Check False Positive (FP)
%                                FP = FP + 1;
%                            end
%                        elseif GT_GW(row,col)== 0 || MTD_GW(row,col)== 0                    
%                            %if GT_GW(row,col)==MTD_GW(row,col) %Check True Negative (TN)
%                            if (GT_GW(row,col)==0 & MTD_GW(row,col)==0)
%                                TN = TN + 1;
%                            %if GT_GW(row,col)~=MTD_GW(row,col) %Check False Negative (FN)
%                            elseif (GT_GW(row,col)==1 & MTD_GW(row,col)==0)
%                                FN = FN + 1;
%                            end
%                        end
                   end
               end
                
               %Convert both ground_truth and model-result into single rows            
               %combine them for intraclass correlation coefficient calculation            
               CombinedData = [GT_GW(:) MTD_GW(:)];              
               ICA_1 = ICC_A_1(CombinedData); 
               ICC_1 = ICC_C_1(CombinedData);
              
               %FN = (size(GT_GW,1)*size(GT_GW,2))-(TP+FP+TN);
               rndVal = 4;
                TPos = TP;
                FPos = FP;
                TNeg = TN;
                FNeg = FN;
                ACC = iifNan(round((TP+TN)/(TP+FP+TN+FN), rndVal));
                SEN = iifNan(round((TP)/(TP+FN), rndVal));
                SPE = iifNan(round((TN)/(TN+FP), rndVal));
                mACC = iifNan(round((SEN+SPE)/2, rndVal));
                FDR = iifNan(round(FP/(FP+TP), rndVal));
                FOR = iifNan(round(FN/(FN+TN), rndVal));
                AROC = iifNan(round(trapz([SPE, SEN])));
                DSC = iifNan(round((TP+TP)/(TP+TP+FP+FN), rndVal));
                PPR = iifNan(round((TP)/(TP+FP), rndVal));
                MCC = iifNan(round(((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)), rndVal));
                F1S = TP/(TP+((FP+FN)*0.5));

                MthdEval{fdex, (m + mLent - 1)} = [TPos, FPos, TNeg, FNeg, ACC, mACC,ICA_1, ICC_1, SEN, SPE, FDR, FOR, DSC, PPR, MCC, F1S, AROC];
                CumSum{1,n} = CumSum{1, n} + MthdEval{fdex, (m + mLent - 1)};
                MthdEval{(img_NUM/mLent)+1, (m + mLent - 1)} = (CumSum{1,n});
                MthdEval{(img_NUM/mLent)+2, (m + mLent - 1)} = (CumSum{1,n}/fdex);


                %disp([num2str(n), '|\t\t', num2str(TPos), '|', num2str(FPos), '|', num2str(TNeg), '|', num2str(FNeg), '|', num2str(ACC), '|', num2str(SEN),'|', num2str(SPE),'|', num2str(DSC),'|', num2str(NPR)]);
                fprintf('%s\t%6d\t%6d\t%6d\t%6d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',...
                    MthdEvalName{m}(1:4), TPos, FPos, TNeg, FNeg, ACC, mACC,ICA_1, ICC_1, SEN, SPE, FDR, FOR, DSC, PPR, MCC, F1S, AROC);
           end
        end
        fdex=fdex+1;
    %    pause(1)
    end
    SeqMthdEval{Seq, 1}=MthdEval;
    
%   SeqPlot1(Seq,:)=cell2mat(MthdEval(end, 4));
%   SeqPlot2(Seq,:)=cell2mat(MthdEval(end, 5));
fprintf('\n\n');
end

%%Generate Plots
% figure (m+1); 
% plot(SeqPlot1(:,7), SeqPlot1(:,8))
% [Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,'true');
% plot(X,Y)
% xlabel('False positive rate') 
% ylabel('True positive rate')
% title('ROC for Classification by Logistic Regression')
