function [ ys ] = runKuramoto( C,D,t_max,approx,useGPU )
% Kuramoto Model:
% C                    = Matrix of coupling weights (NxN) between pairs of regions (can be directed and weighted) (targetNeuronxSourceNeuron)
% D                    = Matrix of distances (in mm) (NxN) (targetNeuronxSourceNeuron)
% k                    = Coupling strength
% f                    = frequency of all oscillators (Hz)
% v                    = Velocity (in m/s)
% t_max                = Total time of simulated activity (in seconds)
% dt                   = integration step (smaller than smaller delays) (in seconds) (eg. 0.0001 s = 0.1 ms)
% sampling             = sampling for saving simulated activity (eg. 10) => 10*dt  = 1ms
% sig_n                = standard deviation of noise in the phases (can be zero, in radians)
% d                    = initial seconds to remove (transient dynamics) (in seconds)
% approx               = use approximation of weak coupling (true/false)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=400;
f=60;
v=1.7;
dt=0.0001;
sampling=10000;
sig_n=0;
d=0;
verbose=true;

sig_noise  = sig_n*sqrt(dt);    % Scale noise per time step !!! Should we scale with f?

timeDelays = round(D/(v*1e3));
stepsDelay = round(D/(v*dt*1e3));
C = k*C*dt;

N = size(C,1);
numIters = t_max/dt;
y = 2*pi*rand(N,1, 'double');
ys = zeros(N,numIters/sampling, 'double');

if useGPU
    y = gpuArray(y);
    C = gpuArray(C);
    stepsDelay = gpuArray(stepsDelay);
    timeDelays = gpuArray(timeDelays);
end

if ~approx
    ringBufferSize = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays
    if useGPU
        ringBuffer = parallel.gpu.GPUArray.rand([ringBufferSize, N],'double');
    else
        ringBuffer = 2*pi*rand([ringBufferSize, N],'double');
    end
end

tic;
t=0;
for i=0:numIters-1
    t = t + dt;
    
    %update state variable:
    if approx
      y = y + sum(C .* sin( bsxfun(@minus,y',y) - timeDelays ) ,2);
    else
      curPosInRingBuff = mod(i,ringBufferSize);
      yDelayed = ringBuffer(sub2ind(size(ringBuffer),mod(curPosInRingBuff-1-stepsDelay,ringBufferSize)+1,repmat(1:size(ringBuffer,2),[N 1])));
      y = y + sum(C .* sin( bsxfun(@minus,yDelayed,y)) ,2);
    end
%     if approx
% %         tmp0 = y*parallel.gpu.GPUArray.ones(1,size(y,1));
% %         tmp1 = tmp0' - tmp0 - timeDelays;
% tmp1 = bsxfun(@minus,y,y') - timeDelays;
%         tmp2 = sin( tmp1 );
%         tmp3 = C .* tmp2;
%         y = y + sum( tmp3,2);
%     else
%         curPosInRingBuff = mod(i,ringBufferSize);
%         yDelayed = ringBuffer(sub2ind(size(ringBuffer),mod(curPosInRingBuff-1-stepsDelay,ringBufferSize)+1,repmat(1:size(ringBuffer,2),[N 1])));
%         y = y + sum(C .* sin( yDelayed - y*ones(1,size(y,1))) ,2);
%     end
    y = y + 2*pi*f*dt;
%     if useGPU
%         y = y + sig_noise*parallel.gpu.GPUArray.randn([N,1],'double');
%     else
%         y = y + sig_noise*randn([N,1],'double');
%     end
    
    %update ring buffer
    if ~approx
        ringBuffer(curPosInRingBuff+1,:) = y;
    end
    
    %save every sampling steps:
    if mod(i,sampling)==0
        ys(:,i/sampling+1) = gather(y);
    end
    
    %print progress:
    if verbose && mod(i,round(numIters/100))==0
        fprintf(['t = ' num2str(t) ' of ' num2str(t_max) ' toc=' num2str(toc) '\n']);
    end
end

if d>0
    numSamplesDelete = round(d/(dt*sampling));
    ys = ys(numSamplesDelete+1:end);
end

end

