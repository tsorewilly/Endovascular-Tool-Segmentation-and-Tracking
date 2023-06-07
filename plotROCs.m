%clear; clc; load('EvaluationResults.mat')
figure(101); clf;
[dXcum, dYcum, dAUCcum, Count] = deal(0); %[Xc, Yc, Tc, AUCc] = deal(0);%zeros(3,1)

for Seq = 1:12
    for mth = 2:size(MthdEvalName,2)
        [Xc, Yc, Tc, AUCc] = deal(0);
        for Frm = 1 : size(SeqMthdEval{Seq,1},1)-2
            MethodResult_res = (SeqMthdEval{Seq,1}{Frm,mth});
            GroundTruth_res = (SeqMthdEval{Seq,1}{Frm,1});
            if(nnz(SeqMthdEval{Seq,1}{Frm,1}))==0
                continue
            end

            GroundTruth_res(GroundTruth_res>0)=1; %Change to two-class data
            MethodResult_res(MethodResult_res>0)=1; %Change to two-class data
            
            actual=(double(GroundTruth_res(:)))>0;
            predic=double(MethodResult_res(:));        
            
            [X,Y,T,AUC] = perfcurve(categorical(actual), predic, 'true');
            Xc = Xc + X(2);
            Yc = Yc + Y(2);
            Tc = Tc + T(1);
            AUCc = AUCc + AUC;
        end
        Count = Count + Frm;
        dXcum(Seq, mth) = Xc/Frm;
        dYcum(Seq, mth) = Yc/Frm;
        dAUCcum(Seq, mth) = AUCc/Frm;
        dTcum(Seq, mth) = Tc/Frm;
    end
end

figure(101); hold on; plot([0 (dXcum) 1], [0 mean(dYcum) 1], 'LineWidth', 2); %Plot ROC Curve

xlabel('False positive rate'); ylabel('True positive rate')
lgd = legend({'Kang et al.[40]','Qian et al.[42]','Proposed Method','Cervantes-SSG  et al.[33]',...
    'Azzopardi et al.[34]','Nguyen et al.[41]','Vlachos & Dermatas [35]','Heneghan et al.[36]', ...
    'CoyeFilter [38]','Chanwimaluang & Fan [39]','Location','southeast'}, 'FontSize',15, 'FontWeight', 'Bold');
title (lgd, 'Comparison of ROC curves');    grid on;    hold off;