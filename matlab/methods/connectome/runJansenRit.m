function [ xColl, driverColl ] = runJansenRit( StartStates,Drivers,DrivFreqs,DrivPO,DrivStart,DrivDur,C,D,k,v,tMax,dt,d,snr,sampling,verbose,JRParams )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Scale Connectivity and Delay matrix
D = round(D/(v*dt*1e3));
C = k*C;

% initialize stuff
N = size(C, 1);
numIters = ceil(tMax/dt);
xColl = zeros(N, numIters/sampling, 'double');
driverColl = zeros(length(DrivPO), numIters/sampling, 'double');

% build buffer to realize delays in network
ringBufferSize = fix(max(D(:)))+1;
disp(['size of ring buffer: ' num2str(ringBufferSize)]);
ringBuffer = squeeze(StartStates(:,1,end-ringBufferSize+1:end))';

%% Run simulation

t = 0;
a = 0;
x = squeeze(StartStates(:,:,end));
inpAvg = 0;
for i = 0:numIters-1
    t = t + dt;
    
    % realize global model delays D
    curPosInRingBuff = mod(i,ringBufferSize);
    
    lookupInBuffer = mod(curPosInRingBuff-1-D,ringBufferSize)+1;
    linInd = sub2ind(size(ringBuffer),lookupInBuffer,repmat(1:size(ringBuffer,2),[N 1]));
    x1Delayed = ringBuffer(linInd);
    
    % extrinsic input updates
    if t > DrivStart && t < DrivStart + DrivDur
        drivPhase = sin(t*2*pi*DrivFreqs + DrivPO);
    else
        drivPhase = zeros(1,length(DrivPO));
    end
    driver = Drivers .* repmat(drivPhase,size(Drivers,1),1);
    driver = sum(driver,2);
    
    inp_tmp = sum(C .* WTP(x1Delayed', JRParams.e0, JRParams.u0, JRParams.r), 2);
    inp = 220 - inpAvg + snr * randn(N, 1) + inp_tmp;
    inpAvg = 0.9999 * inpAvg + 0.0001 * inp_tmp;
    
    % intrinsic state variable updates
    [ dx,da ] = JansenRitUpdate(x,a,driver,inp,JRParams.Ke,JRParams.Ki,JRParams.Ka,JRParams.Cpe,JRParams.Cpi,JRParams.Cep,JRParams.Cip,JRParams.Cii,JRParams.de1,JRParams.de2,JRParams.di1,JRParams.di2,JRParams.e0,JRParams.u0,JRParams.r,JRParams.S0);
    x = x + dt*dx;
    a = a + dt*da;
    
    %update ring buffer
    ringBuffer(curPosInRingBuff+1,:) = x(:,1);
    
    %save every sampling steps:
    if mod(i,sampling)==0 && t > d
        xColl(:,i/sampling + 1) = x(:,1);
        driverColl(:,i/sampling + 1) = drivPhase;
    end
    
    %print progress:
    if verbose && mod(i,round(numIters/100))==0
        fprintf(['t = ' num2str(t) ' of ' num2str(tMax) '\n']);
    end
    
end

xColl = xColl(:,(d/(dt*sampling))+1 : end);
driverColl = driverColl(:,(d/(dt*sampling))+1 : end);
