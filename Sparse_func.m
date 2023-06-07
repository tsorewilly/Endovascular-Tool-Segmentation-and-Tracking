function [Stats] = Sparse_func(ImConnected, nTests, Draw, Verbose, Im,eps, BBNorm, BFSNorm, NodeSampleNorm, LargestCE)
%ImConnected: tells the program if the image is connected or not
%nTests: Number of tests to perform
%Draw: Enable/disable drawing
%Verbose: Verbosity level of function
%Im: Input image
%eps: Image epsilon
%BBNorm: Normalization factor of back-bone matrix size
%BFSNorm: Normalization factor of BFS size
%NodeSampleNorm:Normalization factor of number of sampled pixels
%LargestCE: Size of largest connectivity element
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stats description
%%%%%%%%%%%%%%%%%%%%%%%%%%
    Stats = zeros(6,1);
    %Stats(1,1) = Number of successful tests
    %Stats(2,1) = Number of failed tests
    %Stats(3,1) = average number of Nodes covered per test
    %Stats(4,1) = number of tests completed in BB phase
    %Stats(5,1) = number of tests completed in 2nd phase
    %Stats(6,1) = Algorithm complexity in worst-case
    if ImConnected == 1
        fprintf('Program does not support connected images\n');
        return;
    end
    size(Im);

    % Measure w(M) and calculate epsilon
    Hamm = length(find(Im == 1));
    if(Verbose == 1)
        fprintf('Hamming distance = %d\n', Hamm);
    end

    %build backbone graph
    BBsize = round(sqrt(Hamm)/BBNorm) %width of each backbone sub-matrix
    if BBsize > min(size(Im,1), size(Im,2))
        BBsize = min(size(Im,1), size(Im,2))
    end
    if BBsize == 0
    BBsize = 1;
    end
    %number of BB submatrices along the X dimension
    BBdimX = ceil((size(Im,2) / BBsize))
    %number of BB submatrices along the Y dimension
    BBdimY = ceil((size(Im,2) / BBsize))
    if(Verbose == 1)
        fprintf('backbone is %d x %d submatrices, each of size %d\n', BBdimY, BBdimX, BBsize);
    end

    BBsampleNum = round(sqrt(Hamm)*log(Hamm)*BBNorm*BBNorm); 
    %increase number of sampled pixels by BBNorm^2 since there are BBNorm^2 more submatrices
    if BBsampleNum > Hamm
        BBsampleNum = Hamm;
    end

    if BBsampleNum == 0
        BBsampleNum = 1;
    end

    if(Verbose == 1)
        fprintf('Sampling %d one-pixels of matrix\n', BBsampleNum);
    end
    Stats(6,1) = BBsampleNum + LargestCE;

    for Test = 1:nTests
        BBgraph = zeros(BBdimY, BBdimX);
        
        %Randomally select BBsampleNum 1-pixels from image
        NodesInTest = BBsampleNum;
        pVect = find(Im == 1);
        randidx = randperm(length(pVect));
        IdxSel = pVect(randidx(1:BBsampleNum));
        MarkNodes = Im;
        MarkNodes = double(MarkNodes);
        
        %array of coordinates of selected nodes
        SelNodes = zeros(ceil(BBsampleNum),2);
        [row col] = ind2sub(size(Im),IdxSel(1:BBsampleNum));
        SelNodes(:,:) = [row col];
        
        %fill backbone graph according to sampled 1-pixels
        %save the coordinates of every 1-pixel that is the first to fill a sub-matrix
        FirstPixinBB = [0 0];
        for k=1:BBsampleNum
            CoorY = ceil(SelNodes(k,1)/BBsize);
            CoorX = ceil(SelNodes(k,2)/BBsize);
            if BBgraph(CoorY, CoorX) == 0
                FirstPixinBB = [FirstPixinBB; SelNodes(k,1) SelNodes(k,2)];
                BBgraph(CoorY, CoorX) = 1;
            end
        end
        
        FirstPixinBB = FirstPixinBB(2:end, :);
        %check connectivity of BBgraph
        StartNd = [ceil(SelNodes(1,1)/BBsize) ceil(SelNodes(1,2)/BBsize)];
        [BBgraphRemain IsConn BBgraphNnodeNum] = MatlabDll(StartNd(1),StartNd(2),BBgraph, BBdimX*BBdimY, 1, Verbose);
        
        if BBgraphNnodeNum ~= length(find(BBgraph == 1))
            if Verbose == 1
                fprintf('Graph is not connected\n');
            end
            Stats(1,1) = Stats(1,1) + 1;
            Stats(3,1) = Stats(3,1)*((Test-1)/Test) + NodesInTest/Test;
            Stats(4,1) = Stats(4,1) + 1;
            continue;
        end
        
        %Perform pixel level BFS
        FineSampleNum = round((log(Hamm)/eps)*BBNorm*BBNorm/NodeSampleNorm);
        if FineSampleNum > length(find(BBgraph == 1))
            FineSampleNum = length(find(BBgraph == 1));
        end
        
        if FineSampleNum < 1
            FineSampleNum = 1;
        end
        %Randomally select BBSampleNum 1-pixels from FirstPixinBB
        randidx = randperm(size(FirstPixinBB, 1));
        randidx = randidx(1:FineSampleNum);
        SelNodes = FirstPixinBB(randidx,:);
        %Perfrom BFS from every sampled 1-pixel
        BFSStop = 8*sqrt(Hamm)/eps/BFSNorm;
        
        if BFSStop > Hamm
            BFSStop = Hamm;
        end
        
        if BFSStop < 1
            BFSStop = 1;
        end
        
        MarkNodes = Im;
        MarkNodes = double(MarkNodes);
        IsConn = 1;
        for Nd=1:length(SelNodes)
            [MarkNodes IsConn Nnodes] = MatlabDll(SelNodes(Nd,1), SelNodes(Nd,2),MarkNodes, BFSStop, Nd, Verbose);
            NodesInTest = NodesInTest + Nnodes;
            if IsConn == 0
                Stats(1,1) = Stats(1,1) + 1;
                Stats(3,1) = Stats(3,1)*((Test-1)/Test) + NodesInTest/Test;
                Stats(5,1) = Stats(5,1) + 1;
                if Verbose == 1
                    fprintf('IsConn is zero\n');
                end
                break;
            end
        end
        if IsConn == 1
        Stats(2,1) = Stats(2,1) + 1;
        Stats(3,1) = Stats(3,1)*((Test-1)/Test) + NodesInTest/Test;
        Stats(5,1) = Stats(5,1) + 1;
        end
    end
end