tileSizeX = 12;
tileSizeY = 12;
fIn = 10;
fOut = 14;
numTilesX = randi(3,1);
numTilesY = randi(3,1);
useGpu = true;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

%% a tiled Conv with just one weight element 1 and one input element 1 should output just one element 1
tiledConv = TiledConv(tileSizeX,tileSizeY,fIn,fOut,useGpu);

W = zeros(tileSizeX+1,tileSizeY+1,fIn,tileSizeX,tileSizeY,fOut);
Wdim1=randi(size(W,1),1);
Wdim2=randi(size(W,2),1);
Wdim3=randi(size(W,3),1);

in = zeros(numTilesX*tileSizeX,numTilesY*tileSizeY,fIn);
inPosInTileX = randi(tileSizeX,1);
inPosInTileY = randi(tileSizeY,1);
inDim1 = (randi(numTilesX,1)-1) * tileSizeX + inPosInTileX;
inDim2 = (randi(numTilesY,1)-1) * tileSizeY + inPosInTileY;
inDim3 = Wdim3;
in(inDim1,inDim2,inDim3) = 1;

Wdim4 = -Wdim1 + inPosInTileX + 1;
Wdim5 = -Wdim2 + inPosInTileY + 1;
Wdim6=randi(size(W,6),1);

Wdim4 = mod(Wdim4-1,size(W,4)) + 1;
Wdim5 = mod(Wdim5-1,size(W,5)) + 1;
W(Wdim1,Wdim2,Wdim3,Wdim4,Wdim5,Wdim6) = 1;

if useGpu
    W = gpuArray(W);
    in = gpuArray(in);
end

out = tiledConv.convW(W,in,'full');

%% check
assert(gather(sum(out(:)==1))==1)

%% check inversion:
tiledConvInv = TiledConv(tileSizeX,tileSizeY,fOut,fIn, useGpu);
Winv = tiledConv.invertW(W);
assert(numel(unique(tiledConv.idsW2Winv)) == numel(tiledConv.idsW2Winv));
out = permute(out,[3 1 4 2 5]);
out = reshape(out,[size(out,1)*size(out,2) size(out,3)*size(out,4) size(out,5)]);
reconstr = tiledConvInv.convW(Winv,out,'valid');
assert(gather(sum(reconstr(:)==1))==1)

%% check if reconstruction is the same as the input
reconstr = permute(reconstr,[3 1 4 2 5]);
reconstr = reshape(reconstr,[size(reconstr,1)*size(reconstr,2) size(reconstr,3)*size(reconstr,4) size(reconstr,5)]);
assert(isequal(gather(reconstr),gather(in)));

%% check backprop error / try to reconstruct W from in and out:
in2 = zeros([size(in,1)+tileSizeX*2 size(in,2)+tileSizeY*2 size(in,3)]);
if useGpu
    in2 = gpuArray(in2);
end
in2(tileSizeX+1:end-tileSizeX,tileSizeY+1:end-tileSizeY,:) = in;
out2 = tiledConv.convW(W,in2,'valid');
Wbackprop = tiledConv.backpropErrorW(in2,out2);
assert(isequal(gather(W),gather(Wbackprop)))

%%
tmp=gather(out2);
size(tmp)
sum(tmp(:))