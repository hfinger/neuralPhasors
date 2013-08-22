function [ ys ] = runRatemodel(C,D,k,tau,v,t_max,dt,sampling,sig_n,d,verbose)
% Rate Model:
% C                    = Matrix of coupling weights (NxN) between pairs of regions (can be directed and weighted) 
% D                    = Matrix of distances (in mm) (NxN)
% k                    = Coupling strength
% tau                  = Time scale (ms)
% v                    = Velocity (in m/s)
% t_max                = Total time of simulated activity (in seconds)
% dt                   = integration step (smaller than smaller delays) (in seconds) (eg. 0.0001 s = 0.1 ms)
% sampling             = sampling for saving simulated activity (eg. 10) => 10*dt  = 1ms
% sig_n                = standard deviation of noise in the phases (can be zero, in radians) 
% d                    = initial seconds to remove (transient dynamics) (in seconds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau        = tau/1e3;

sig_noise  = sig_n*sqrt(dt/tau);    % Scale noise per time step !!! why scale with tau?

stepsDelay = round(D/(v*dt*1e3));
C = k*C*dt/tau;

N = size(C,1);

ringBufferSize = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays
initBuf = sig_n*randn(N,ringBufferSize);
ringBuffer = initBuf';
y = ones(size(C,1),1, 'double');
y = y.*initBuf(:,end);

numIters = t_max/dt;

ys = zeros(N,numIters/sampling, 'double');

t=0;
for i=0:numIters-1
    t = t + dt;
    curPosInRingBuff = mod(i,ringBufferSize);
    yDelayed = ringBuffer(sub2ind(size(ringBuffer),mod(curPosInRingBuff-1-stepsDelay,ringBufferSize)+1,repmat(1:size(ringBuffer,2),[N 1])));
    
    %update state variable:
    dy = sum(C.*yDelayed,2) + sig_noise*randn(N,1);
    y = (1-dt/tau)*y + dy;
    ringBuffer(curPosInRingBuff+1,:) = y;
    
    %save every sampling steps:
    if mod(i,sampling)==0
        ys(:,i/sampling+1) = y;
    end
    
    %print progress:
    if verbose && mod(i,round(numIters/100))==0
        fprintf(['t = ' num2str(t) ' of ' num2str(t_max) '\n']);
    end
end

end

