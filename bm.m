% results = bm([contingency matrix])
% 
% contingency matrix  ->  contingency table of network prediction and training targets
%                             ie. classification x target'
%                                         
% results             <-  structure of results with the following fields
%                             contingencyMatrix:  As provided as argument 1 to function
%                             recall:             Ratio of times class correctly predicted
%                                                 to total times it was actually that class
%                             precision:          Ratio of times class correctly predicted 
%                                                 to total times it was predicted to be that class
%                             weightedAverage:    Ratio of correct predicted classes
%                                                 to total number of cases
%                             F:                  Harmonic mean
%                             Fall:               F for all classes
%                             G:                  Geometric mean
%                             Gall:               G for all classes
%                             bookmakerMatrix:    Matrix of bookmaker results per class
%                             bookmaker:          Bookmaker result for all classes
%
% matlab: sf 16/03/03
% octave: dp 11//03 - also extended to weight k>2 combinations of bm, F and G correctly


function results = bm(cm)

if (size(cm,1) ~= size(cm,2))
    error('Contingency matrix must be square'); 
else
    k = size(cm,1);
end

N = sum(sum(cm));
rprob = sum(cm) ./ N;
pprob = sum(cm') ./ N;

recall = diag(cm)' ./ sum(cm);
precision = diag(cm')' ./ sum(cm');
wav = sum(diag(cm)) ./ N;

F = (2.*precision.*recall) ./ (precision + recall);
G = sqrt(precision.*recall);

% Fall - use weighted harmonic mean 
Fall = 1 / sum(rprob./F); %(assume real distribution)
Fall = k / sum(ones(1,k)./F); %(assume equiprobable)

% Gall - use weighted geometric mean
Gall = prod(G.^rprob); %(assume real distribution)
Gall = (prod(G)^(1/k)); %(assume equiprobable)

mask = diag(ones(1,k));
maskc = reshape(mask,k*k,1);
ind = find(maskc==0);
maskc(ind) = -1;
mask = reshape(maskc,k,k);

prob = rprob;
prob = prob(ones(k,1),:)';
probc = reshape(prob,k*k,1);
probc(ind) = 1-probc(ind);
prob = reshape(probc,k,k);

bmcm = cm ./ prob;    
bmcm = bmcm .* mask;
bms = sum(bmcm') / N;
bm = bms * pprob';

%results.contingencyMatrix = cm;
results.N = N;
results.recall = recall;
results.precision = precision;
results.weightedAverage = wav;
results.F = F;
results.Fall = Fall;
results.G = G;
results.Gall = Gall;
%results.bookmakerMatrix = bmcm;
results.bookmakerSum = bms;
results.bookmaker = bm;
