%% Function "eval_metr" calculates the metrics for model validatio and evaluation 
%   Input parameters                                                      
%       * class_1:  Column vector that represents the data of the first class.
%       * class_2:  Column vector that represents the data of the second class.
%       * dispp:    (Optional) If dispp is 1, the ROC Curve will be dispayed inside the active figure.
%       * dispt:    (Optional) If dispt is 1, the optimum threshold parameters obtained will be displayed 
%                                                                         
%   Output variables                                                      
%       * ROC_data: Struct that contains all the curve parameters.        
%           - param:    Struct that contains the cuantitative parameters of the obtained curve, which are:                 %
%               + Threshold:Optimum threshold calculated in order to maximice the sensitivity and specificity
%               values which is colocated in the nearest point to (0,1).
%               + AROC:     Area under ROC curve.                         
%               + Accuracy: Accuracy.                                     
%               + Sensi:    Sensitivity (i.e., recall, hit rate, or true positive rate).                               %
%               + Speci:    Specificity (i.e., selectivity, or true negative rate).                               %
%               + PPV:      Positive predicted value (i.e., precision).   
%               + NPV:      Negative predicted value.                     
%               + FNR:      False negative rate (i.e., miss rate).        
%               + FPR:      False positive rate (i.e., fall-out).         
%               + F1_score: F1 score (harmonic mean of precision and sensitivity).                                 %
%               + MCC:      Matthews correlation coefficient.             
%           - curve:    Matrix that contains the specificity and specificity of each threshold point in columns.          %
% 
%   Example of use:                                                       
%       roc_curve(class_1, class_2);                                      
% ----------------------------------------------------------------------- 
function ROC_data = roc_curve(class_1, class_2, dispp, dispt)

    % Setting default parameters and detecting errors
    if(nargin<4), dispt = 1;    end
    if(nargin<3), dispp = 1;    end
    if(nargin<2), error('Params "class_1" or "class_2" are not indicated.'); end
    class_1 = class_1(:);
    class_2 = class_2(:);
    
    % Calculating the threshold values between the data points
    s_data = unique(sort([class_1; class_2]));          % Sorted data points
    s_data(isnan(s_data)) = [];                 % Delete NaN values
    d_data = diff(s_data);                      % Difference between consecutive points
    if(isempty(d_data)), error('Both class data are the same!'); end
    d_data(length(d_data)+1,1) = d_data(length(d_data));% Last point
    thres(1,1) = s_data(1) - d_data(1);                 % First point
    thres(2:length(s_data)+1,1) = s_data + d_data./2;   % Threshold values
        
    % Calculating the sensitivity and specificity of each threshold
    curve = zeros(size(thres,1),2);
    distance = zeros(size(thres,1),1);
    for id_t = 1:1:length(thres)
        TP = length(find(class_2 >= thres(id_t)));    % True positives
        FP = length(find(class_1 >= thres(id_t)));    % False positives
        FN = length(find(class_2 < thres(id_t)));     % False negatives
        TN = length(find(class_1 < thres(id_t)));     % True negatives
        
        curve(id_t,1) = TP/(TP + FN);   % Sensitivity
        curve(id_t,2) = TN/(TN + FP);	% Specificity
        
        % Distance between each point and the optimum point (0,1)
        distance(id_t)= sqrt((1-curve(id_t,1))^2+(curve(id_t,2)-1)^2);
    end
    
    % Optimum threshold and parameters
    [~, opt] = min(distance);
    TP = length(find(class_2 >= thres(opt)));    % No. true positives
    FP = length(find(class_1 >= thres(opt)));    % No. false positives 
    FN = length(find(class_2 < thres(opt)));     % No. false negatives                                 
    TN = length(find(class_1 < thres(opt)));     % No. true negatives       
    
    % Output parameters
    param.Sensi = curve(opt,1);                 % Sensitivity
    param.Speci = curve(opt,2);                 % Specificity
    param.AROC  = abs(trapz(1-curve(:,2),curve(:,1))); % Area under curve
    param.Accuracy = (TP+TN)/(TP+TN+FP+FN);     % Accuracy
    param.PPV   = TP/(TP+FP);                   % Positive predictive value
    param.NPV   = TN/(TN+FN);                   % Negative predictive value
    param.FNR   = FN/(FN+TP);                   % False negative rate
    param.FPR   = FP/(FP+TN);                   % False positive rate
    param.F1_score = 2*TP/(2*TP+FP+FN);         % F1 score
    param.MCC   = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));  % Matthews correlation coefficient
    
    param.TP = TP;    % No. true positives
    param.FP = FP;    % No. false positives 
    param.FN = FN;    % No. false negatives                                 
    param.TN = TN;    % No. true negatives  
    
    % Plotting if required
    if(dispp == 1)
        figure(10001)
        fill_color = [11/255, 208/255, 217/255];
        fill([1-curve(:,2); 1], [curve(:,1); 0], fill_color,'FaceAlpha',0.5);
        hold on; plot(1-curve(:,2), curve(:,1), '-b', 'LineWidth', 2);
        hold on; plot(1-curve(opt,2), curve(opt,1), 'or', 'MarkerSize', 10);
        hold on; plot(1-curve(opt,2), curve(opt,1), 'xr', 'MarkerSize', 12);
        hold off; axis square; grid on; xlabel('1 - Specificity'); ylabel('Sensitivity');
        title(['AROC = ' num2str(param.AROC)]);
    end
    
    % AROC warning
    if param.AROC < 0.5
        warning('Since AROC is less than 0.5, you should swap the classes: roc_curve(class_2,class_1).');
    end
    
    % Log screen parameters if required
    if(dispt == 1)
        fprintf('\n ROC CURVE PARAMETERS\n');
        fprintf(' ------------------------------\n');
        fprintf('  - Sensitivity:  %.4f\n', param.Sensi);
        fprintf('  - Specificity:  %.4f\n', param.Speci);
        fprintf('  - AROC:         %.4f\n', param.AROC);
        fprintf('  - Accuracy:     %.4f\n', param.Accuracy);
        fprintf('  - PPV:          %.4f\n', param.PPV);
        fprintf('  - NPV:          %.4f\n', param.NPV);
        fprintf('  - FNR:          %.4f\n', param.FNR);
        fprintf('  - FPR:          %.4f\n', param.FPR);
        fprintf('  - F1 score:     %.4f\n', param.F1_score);
        fprintf('  - MCC:          %.4f\n', param.MCC);
        fprintf(' \n');
    elseif(dispt == 2)
%         fprintf('ACC\t\t\tSEN\t\t\tSPE\t\t\tAROC\t\tMCC\t\tF1\t\tFPR\t\tFNR\t\tNPV\t\tPPV\t\t\n');
%         fprintf('%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f',...
%             param.Accuracy, param.Sensi,param.Speci, param.AROC,param.PPV, param.NPV, param.FNR, param.FPR, ...
%             param.F1_score,param.MCC);
%         fprintf(' \n');
    end
    
    % Assinging parameters and curve data
    ROC_data.param = param;
    ROC_data.curve = curve;
end