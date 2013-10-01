clear all;

useGpu = true;
tileSizeX = 8;
tileSizeY = 8;
fIn = 3;
fOut = 10;
numTilesX = 30;
numTilesY = 30;
numIters = 100;

%% init rand stream:
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

%% init tiledConv object
tiledConv = TiledConv(tileSizeX,tileSizeY,fIn,fOut,useGpu);
tiledConvInv = TiledConv(tileSizeX,tileSizeY,fOut,fIn, useGpu);

%% initialize Weights for tiledConv:
W = randn([tileSizeX+1,tileSizeY+1,fIn,tileSizeX,tileSizeY,fOut],'double');
Winv = tiledConv.invertW(W);

%% initialize input for tiledConv:
in = rand([numTilesX*tileSizeX,numTilesY*tileSizeY,fIn],'double');

if useGpu
    W = gpuArray(W);
    in = gpuArray(in);
end

%% iterate several times up and down:

WfullSplitted = tiledConv.W2WfullSplitted(W);
WfullSplittedInv = tiledConvInv.W2WfullSplitted(Winv);

tic;
for k=1:numIters
    
    in = (in>0.5);
    
    out5D = tiledConv.convWfullSplitted(WfullSplitted,in,'full');
    out = permute(out5D,[3 1 4 2 5]);
    out = reshape(out,[size(out,1)*size(out,2) size(out,3)*size(out,4) size(out,5)]);
    in5D = tiledConvInv.convWfullSplitted(WfullSplittedInv,out,'valid');
    in = permute(in5D,[3 1 4 2 5]);
    in = reshape(in,[size(in,1)*size(in,2) size(in,3)*size(in,4) size(in,5)]);
    Wbackprop = tiledConvInv.backpropErrorW(out,in5D);
    
    %print progress:
    fprintf(['k = ' num2str(k) ' of ' num2str(numIters) ' toc=' num2str(toc) '\n']);
end
disp(toc)

% in = gather(in)
% save('inCPU.mat','in')