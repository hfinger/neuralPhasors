clear all;

tileSizeX = 16;
tileSizeY = 16;
fIn = 10;
fOut = 10;
numTilesX = 20;
numTilesY = 15;
useGpu = false;
numIters = 20;

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
in = randn([numTilesX*tileSizeX,numTilesY*tileSizeY,fIn],'double');

if useGpu
    W = gpuArray(W);
    in = gpuArray(in);
end

%% iterate several times up and down:
tic;
for k=1:numIters
    out5D = tiledConv.convW(W,in,'full');
    out = permute(out5D,[3 1 4 2 5]);
    out = reshape(out,[size(out,1)*size(out,2) size(out,3)*size(out,4) size(out,5)]);
    in5D = tiledConvInv.convW(Winv,out,'valid');
    in = permute(in5D,[3 1 4 2 5]);
    in = reshape(in,[size(in,1)*size(in,2) size(in,3)*size(in,4) size(in,5)]);
    Wbackprop = tiledConvInv.backpropErrorW(out,in5D);
    
    %print progress:
    fprintf(['k = ' num2str(k) ' of ' num2str(numIters) ' toc=' num2str(toc) '\n']);
end
disp(toc)
