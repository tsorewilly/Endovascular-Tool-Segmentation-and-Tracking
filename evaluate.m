clear; clc;
warning off;
FolderDir = 'G:\Documents\My Postdoc Work\Publications\IJPRAI\Evaluation\';

MthdEvalName = {'GroundTruth','01-Proposed','03-Cervantes_SSG_Otsu', '04-Cosfire',...
  '04-FastMScale','05-FstFramLinStruc','06-Heneghan','07-Heneghan-WTH1',...
  '08-CoyeFilter', '09-PrincipalCurv'};

scanBox = [221 753 1308 1206];
mLent = length(MthdEvalName);

for Seq = 1:12
    %Load Image Frames    
    Lists=dir(fullfile(FolderDir, ['Seq', num2str(Seq),'-*']));    
    imageLists = {Lists.name};
    img_NUM = length(imageLists);
    ind=1; fdex=1; MthdEval={};
    fprintf('Id   TPos      FPos          TNeg      FNeg        ACC     mACC        SEN        SPE           DSC         PPR            MCC         AROC\n')%

    CumSum = {zeros(1,12), zeros(1,12)};        
    for f = 1:mLent:img_NUM
        for m=1:mLent
            if m ==1
                MthdEval{fdex,m} = imcrop(imread(fullfile(FolderDir, imageLists{ind})), scanBox);
            else
                MthdEval{fdex,m} = imcrop(rgb2gray(imread(fullfile(FolderDir, imageLists{ind}))), scanBox);
            end
           scanObject=MthdEval{fdex,m}; scanObject(scanObject>0)=1; MthdEval{fdex,m} = scanObject;
           z1(m)=nnz(scanObject);
           %figure (m); imshow(MthdEval{fdex,m}, []); title(imageLists{f+m-1})  
           ind=ind+1;

           %Now compare against ground truth
           if m > 1           
               GT_GW  = MthdEval{fdex,1};
               MTD_GW = MthdEval{fdex,m};
              [TP, TN, FP, FN] = deal(0);
               for row = 1:size(GT_GW,1)
                   for col = 1:size(GT_GW,1)
                       %Analyze just guidewire pixels to compute evaluation metrics
                       if GT_GW(row,col)== 1 || MTD_GW(row,col)== 1                    
                           if (GT_GW(row,col)==1 & MTD_GW(row,col)==1) %Check True Positive (TP)
                               TP = TP + 1;
                           elseif (GT_GW(row,col)==0 & MTD_GW(row,col)==1) %Check False Positive (FP)
                               FP = FP + 1;
                           end
                       elseif GT_GW(row,col)== 0 || MTD_GW(row,col)== 0                    
                           %if GT_GW(row,col)==MTD_GW(row,col) %Check True Negative (TN)
                           if (GT_GW(row,col)==0 & MTD_GW(row,col)==0)
                               TN = TN + 1;
                           %if GT_GW(row,col)~=MTD_GW(row,col) %Check False Negative (FN)
                           elseif (GT_GW(row,col)==1 & MTD_GW(row,col)==0)
                               FN = FN + 1;
                           end
                       end
                   end
               end
                
               FN = (size(GT_GW,1)*size(GT_GW,2))-(TP+FP+TN);
               rndVal = 4;
                TPos = TP;
                FPos = FP;
                TNeg = TN;
                FNeg = FN;
                ACC = iifNan(round((TP+TN)/(TP+FP+TN+FN), rndVal));
                SEN = iifNan(round((TP)/(TP+FN), rndVal));
                SPE = iifNan(round((TN)/(TN+FP), rndVal));
                mACC = iifNan(round((SEN+SPE)/2, rndVal));
                FOR = iifNan(round(FN/(FN+TN), rndVal));
                AROC = iifNan(round(trapz([SEN, SPE])));
                DSC = iifNan(round((TP+TP)/(TP+TP+FP+FN), rndVal));
                PPR = iifNan(round((TP)/(TP+FP), rndVal));
                MCC = iifNan(round(((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)), rndVal));

                MthdEval{fdex, (m + mLent - 1)} = [TPos, FPos, TNeg, FNeg, ACC, mACC, SEN, SPE, DSC, PPR, MCC, AROC];
                CumSum{1,m-1} = CumSum{1, m-1} + MthdEval{fdex, (m + mLent - 1)};
                MthdEval{(img_NUM/mLent)+1, (m + mLent - 1)} = (CumSum{1,m-1});
                MthdEval{(img_NUM/mLent)+2, (m + mLent - 1)} = (CumSum{1,m-1}/fdex);


                %disp([num2str(m-1), '|\t\t', num2str(TPos), '|', num2str(FPos), '|', num2str(TNeg), '|', num2str(FNeg), '|', num2str(ACC), '|', num2str(SEN),'|', num2str(SPE),'|', num2str(DSC),'|', num2str(NPR)]);
                fprintf('%d\t%d\t\t%d\t\t%d\t\t%d\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n',...
                    m-1, TPos, FPos, TNeg, FNeg, ACC, mACC, SEN, SPE, DSC, PPR, MCC, AROC);
           end
        end
        fdex=fdex+1;
    %    pause(1)
    end 
    SeqMthdEval{Seq, 1}=MthdEval;
    SeqPlot1(Seq,:)=cell2mat(MthdEval(end, 4));
    SeqPlot2(Seq,:)=cell2mat(MthdEval(end, 5));
end

%%Generate Plots
figure (m+1); 
plot(SeqPlot1(:,7), SeqPlot1(:,8))
[Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,'true');
plot(X,Y)
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression')
