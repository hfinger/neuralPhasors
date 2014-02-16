function [ ys ] = runKuramoto( C,D,k,f,v,t_max,dt,sampling,sig_n,d,verbose,approx,invertSin,startState )
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

sig_noise  = sig_n*sqrt(dt);    % Scale noise per time step !!! Should we scale with f?

stepsDelay = round(D/(v*dt*1e3));
timeDelays = round(D/(v*1e3));
C = k*C*dt;
if invertSin
  C = -C;
end

N = size(C,1);
if isempty(startState)
  y = 2*pi*rand(N,1, 'double');
else
  y = startState(:,end);
end


if ~approx
  ringBufferSize = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays
  disp(['size of ring buffer: ' num2str(ringBufferSize)]);
  if isempty(startState)
    initBuf = 2*pi*rand(N,ringBufferSize);
  else
    initBuf = startState(:,end-ringBufferSize+1:end);
  end
  ringBuffer = initBuf';
end

numIters = t_max/dt;

ys = zeros(N,numIters/sampling, 'double');

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
    y = y + 2*pi*f*dt;
    y = y + sig_noise*randn(N,1);
    
    %update ring buffer
    if ~approx
      ringBuffer(curPosInRingBuff+1,:) = y;
    end
    
    %save every sampling steps:
    if mod(i,sampling)==0
        ys(:,i/sampling+1) = y;
    end
    
    %print progress:
    if verbose && mod(i,round(numIters/100))==0
        fprintf(['t = ' num2str(t) ' of ' num2str(t_max) '\n']);
    end
end

if d>0
  numSamplesDelete = round(d/(dt*sampling));
  ys = ys(:,numSamplesDelete+1:end);
end

end

