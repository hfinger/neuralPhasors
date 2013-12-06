function [ Wforw, Wback ] = trainAE( trainData, useGpu, tileSizeX, tileSizeY, fIn, fOut, numIters, lrate )

%% init tiledConvForw object
tiledConvForw = TiledConv(tileSizeX,tileSizeY,fIn,fOut,useGpu);
tiledConvBack = TiledConv(tileSizeX,tileSizeY,fOut,fIn, useGpu);

%% initialize Weights for tiledConvForw:
Wforw = randn([tileSizeX+1,tileSizeY+1,fIn,tileSizeX,tileSizeY,fOut],'double');
Wback = tiledConvForw.invertW(Wforw);
if useGpu
    Wforw = gpuArray(Wforw);
    Wback = gpuArray(Wback);
end


%% iterate several times up and down:
for k=1:numIters
    
    error = 0;
    
    if useGpu
        WForwGrad = parallel.gpu.GPUArray.zeros(size(Wforw),'double');
        WBackGrad = parallel.gpu.GPUArray.zeros(size(Wback),'double');
    else
        WForwGrad = zeros(size(Wforw),'double');
        WBackGrad = zeros(size(Wback),'double');
    end
    
%     WforwInv = tiledConvForw.invertW(Wforw);
    WbackInv = tiledConvBack.invertW(Wback);
    
    WfullSplittedForw = tiledConvForw.W2WfullSplitted(Wforw);
    WfullSplittedBack = tiledConvBack.W2WfullSplitted(Wback);
    
%     WfullSplittedForwInv = tiledConvBack.W2WfullSplitted(WforwInv);
    WfullSplittedBackInv = tiledConvForw.W2WfullSplitted(WbackInv);

    for s=1:size(trainData,4)
        if useGpu
            in = gpuArray(trainData(:,:,:,s));
        else
            in = trainData(:,:,:,s);
        end
        in5D = tiledConvForw.compact2expanded(in);

        hid5D = tiledConvForw.convWfullSplitted(WfullSplittedForw,in,'valid');
        hidAct5D = sigm(hid5D);
        hid = tiledConvForw.expanded2compact(hid5D);
        hidAct = tiledConvForw.expanded2compact(hidAct5D);

        recons5D = tiledConvBack.convWfullSplitted(WfullSplittedBack,hidAct,'valid');
        recons5DAct = sigm(recons5D);
        
        %% error backprop
        reconsDelta = in5D(2:end-1, 2:end-1, :, :, :) - recons5DAct;
        
        error = error + sum(reconsDelta(:).^2);
        
        reconsDelta = - reconsDelta * lrate;
        reconsDelta = reconsDelta .* (recons5DAct.*(1-recons5DAct));
        reconsDelta3D = tiledConvForw.expanded2compact(reconsDelta);
        hidDelta = tiledConvForw.convWfullSplitted(WfullSplittedBackInv,reconsDelta3D,'full');
        hidDelta = hidDelta .* (hidAct5D.*(1-hidAct5D));
        
        WForwGrad = WForwGrad - tiledConvForw.backpropErrorW(in,hidDelta);
        WBackGrad = WBackGrad - tiledConvBack.backpropErrorW(hid,reconsDelta);

%           disp(s)
    end
    
    Wforw = Wforw + WForwGrad;
    Wback = Wback + WBackGrad;
    
    %print progress:
    fprintf(['k = ' num2str(k) ' of ' num2str(numIters) ' toc=' num2str(toc) ' error=' num2str(error) '\n']);
end

Wforw = gather(Wforw);
Wback = gather(Wback);

end

function X = sigm(P)
    X = 1./(1+exp(-P));
end