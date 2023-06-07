clear;  close all
%% Sparse_main.m %%
% This application is an implementation of connectivity testing algorithm
% for sparse images published by Dana Ron and Gilad Tsur
%Add Mex file DLL to path
path(path,'../MatlabDll_BFS/Debug');
ImConnected = 0; %Determines if Image is connected or not
nTests = 100;
Draw = 1; %enable/disable drawing
Verbose = 0; %verbosity level
eps = 0.012; %epsilon of image
LargestCE = 1335; %size of largest connectivity element
fileName = 'Seq1-Frm155-Col4'; %image file name
showIm = 0;
SaveImage = 0; %enable/disable image saving
thresh = 0.1;
EdgeCutoff = 10;
doMed = 0; %perform median filtering before applying threshold
doMedSecond = 0; %perform median filtering after applying threshold
doFilterSmallParticles = 1; %enable/disable filtering of sample particles
PartFiltSize = 5;
PartFiltThresh = 12;
doThinning = 1; %enable/disable thinning
median_size = 9;
%Pre-process image
file=['G:\Documents\My Postdoc Work\Publications\IJPRAI\Results\', fileName, '.png'];
Im = im2single(imread(char(file)));
if(size(Im,3) > 1)
    Im = Im(:,:,1);
end

if showIm == 1
    subplot(221)
    imshow(Im);
end

%Perform median filtering
if doMed == 1
    if(median_size > 0)
        Im = medfilt2(Im, [median_size median_size]);
    end
end

if showIm == 1
    subplot(222)
    imshow(Im);
end

%Sobel filtering
ImGradX = imfilter(Im, [-1 0 1; -2 0 2; -1 0 1], 'symmetric');
ImGradY = imfilter(Im, [1 2 1; 0 0 0; -1 -2 -1], 'symmetric');
Im = sqrt(ImGradX.^2 + ImGradY.^2);
if showIm == 1
    subplot(223)
    imshow(Im);
end

Im(find(Im > thresh)) = 1;
Im(find(Im <= thresh)) = 0;

%Perform median filtering again
if doMedSecond == 1
    median_size = 3;
    if(median_size > 0)
        Im = medfilt2(Im, [median_size median_size]);
    end
end

if doFilterSmallParticles == 1
    halfFiltsize = (PartFiltSize-1)/2;
    for row= PartFiltSize:size(Im,2) - PartFiltSize
        for col= PartFiltSize:size(Im,1) - PartFiltSize
            box = Im(col-halfFiltsize:col+halfFiltsize, row-halfFiltsize:row+halfFiltsize);
            if length(find(box == 1)) < PartFiltThresh
                Im(col,row) = 0;
            end
        end
    end
end

if doThinning == 1
    ImCleaned = Im;
    for row= 2:size(Im,2) - 1
        for col= 2:size(Im,1) - 1
            if (Im(col+1,row) + Im(col-1,row) + Im(col,row+1) + Im(col,row-1) )== 4
                ImCleaned(col,row) = 0;
            end
        end
    end
    Im = ImCleaned;
end

if EdgeCutoff > 0
    Im = Im(EdgeCutoff:end-EdgeCutoff, :);
    Im = Im(:, EdgeCutoff:end - EdgeCutoff);
end
if showIm == 1
    subplot(224)
    imshow(Im);
end

Im = round(Im);
%Save Image
if SaveImage == 1
fileName=['../images/', fileName, '_preProc_sparse_thresh_',num2str(thresh)];
    imwrite(Im, [fileName, '.jpg'], 'jpeg', 'quality', 100);
end

%Test if image is sparse
SparsityBorder = ((size(Im,1) + size(Im,2)) / 2) ^ (4/3);
if length(find(Im == 1)) < SparsityBorder
    if(Verbose == 1)
        fprintf(['Image is sparse, total of ', num2str(length(find(Im== 1))),' 1-pixels\n']);
    else
        fprintf(['Image is sparse\n']);
    end
else
    if(Verbose == 1)
        fprintf(['Image is dense ', num2str(length(find(Im == 1))), 'pixels. For Sparse requires up to ', num2str(SparsityBorder),'. Exiting\n']);
        return;
    else
        fprintf('Image is dense. Exiting\n');
        return
    end
end

if showIm == 1
    return;
end
%Test if image is connected
%BBsizeNorm = [0.1:0.1:1 1.2:0.2:10];
BBsizeNorm = 1;
nBBsizeNorm = length(BBsizeNorm);
%BFSNorm = [1:10 20:10:100 200 500];
BFSNorm = [1];
nBFSNorm = length(BFSNorm);
NodeSampleNorm = [1:10 20:10:100];
%NodeSampleNorm = [1];
nNodeSampleNorm = length(NodeSampleNorm);
SimResult = zeros(6,30);
%for k=1:nBBsizeNorm
%for k=1:nBFSNorm
for k=1:nNodeSampleNorm
    SimResult(:,k) = Sparse_func(ImConnected, nTests, Draw, Verbose,Im, eps, BBsizeNorm(1), BFSNorm(1), NodeSampleNorm (k), LargestCE);
end

figure(2)
plot(NodeSampleNorm, SimResult(1,1:NodeSampleNorm)/nTests , '—rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',4);
title('Success rate as function of node sample norm');
figure(3)
semilogx(NodeSampleNorm, SimResult(3,1:nNodeSampleNorm) , '—rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',5);
hold on
semilogx(NodeSampleNorm, SimResult(6,1:nNodeSampleNorm), ':rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',5);
