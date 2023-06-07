%% Create datasets

% [X,Y] = ndgrid(1:100,1:100);

% Ground Truth image
% Circle= zeros(100); Circle((X-50).^2 + (Y-50).^2 - 500 <= 0) = 1;

% Test image (parameterised by threshold)
% pkg load statistics % if using octave - needed for 'mvnpdf'
% pkg load image      % if using octave - needed for 'mat2gray'
% Gaussian = mvnpdf([X(:), Y(:)], [45, 45], [500,0;0,500]);
% Gaussian = reshape(Gaussian, size(X));
% Gaussian = mat2gray(Gaussian);
for Seq = 1:12
    for mth = 2:size(MthdEvalName,2)
        %for Frm = 1 : size(SeqMthdEval{Seq,1},1)-2
            Gaussian = (SeqMthdEval{Seq,1}{1, mth});
            Circle = (SeqMthdEval{Seq,1}{1, 1});
            %actual=(double(GroundTruth_res(:)))>0;
            %predic=double(MethodResult_res(:)); 
            %% Generate ROC curve for a range of thresholds
            ThresholdRange = 0 : 0.025 : 1;
            TPs = zeros(size(ThresholdRange));
            FPs = zeros(size(ThresholdRange));
            Ind = 0;
            for Threshold = ThresholdRange 
                Ind = Ind + 1;
                TP = Circle(:) .* (Gaussian(:) > Threshold);
                T  = Circle(:);
                TPR = sum(TP(:)) / sum(T(:));
                TPs(Ind) = TPR;

                FP = (1 - Circle(:)) .* (Gaussian(:) > Threshold);
                N  = (1 - Circle(:));
                FPR = sum(FP(:)) / sum(N(:));
                FPs(Ind) = FPR;
            end

            %% Plot curve
            plot(FPs, TPs, 'linewidth', 3, 'marker', 'o', 'markersize', 10, 'markeredgecolor', 'k', 'markerfacecolor', 'g');
            hold on; 
            plot(ThresholdRange, ThresholdRange, 'r-.', 'linewidth', 3);
        %end
    end
end
axis([0,1,0,1]);

xlabel('False positive rate'); ylabel('True positive rate')
lgd = legend({'Kang et al.[40]','Qian et al.[42]','Proposed Method','Cervantes-SSG  et al.[33]',...
    'Azzopardi et al.[34]','Nguyen et al.[41]','Vlachos & Dermatas [35]','Heneghan et al.[36]', ...
    'CoyeFilter [38]','Chanwimaluang & Fan [39]','Location','southeast'}, 'FontSize',15, 'FontWeight', 'Bold');
title (lgd, 'Comparison of ROC curves');    grid on;    hold off;