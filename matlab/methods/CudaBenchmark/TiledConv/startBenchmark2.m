function executionTime = startBenchmark2( useGpu )
%% A simple benchmark optimizing a tiled convolutional autoencoder...

if nargin<1
  useGpu = false;
end

numIters = 100;

tileSizeX = 5;
tileSizeY = 5;
fIn = 3;
fOut = 10;
numTilesX = 80;
numTilesY = 60;
numSamples = 20;
lrate = 1e-4;

%% init rand stream:
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

trainData = rand([numTilesX*tileSizeX,numTilesY*tileSizeY,fIn,numSamples],'double');

tic;
[ Wforw, Wback ] = trainAE( trainData, useGpu, tileSizeX, tileSizeY, fIn, fOut, numIters, lrate );

executionTime = toc;
disp(['execution time: ' num2str(executionTime)])
