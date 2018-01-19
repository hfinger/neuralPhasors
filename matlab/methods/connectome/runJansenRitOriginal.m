function [ xColl, driverColl ] = runJansenRitOriginal( StartStates,Drivers,DrivFreqs,DrivPO,DrivStart,DrivDur,C,D,k,v,tMax,dt,initSampRem,nVar,nMu,rAvg,netInp,subInp,sampling,verbose,JRParams )
%UNTITLED Function that simulates a number of coupled Jansen-Rit neural
%masses over a certain time. A driver can be included to drive specific
%neural masses externally.
%
%   Input Parameters:
%
%       StartStates - 3-dimensional array with first dimension = number of
%                     neural masses & second dimension = populations, third = timesteps. Used as
%                     initial states of the simulation
%       Drivers     - 2-dimensional array with first dimension = number of
%                     neural masses & second dimension = number of drivers.
%                     Values indicate driving strength.
%       DrivFreqs   - Vector with length = number of drivers that indicates
%                     the frequency of each driver in Hz.
%       DrivPO      - Vector with length = number of drivers that indicates
%                     the phase offset of each driver
%       DrivStart   - Scalar, indicates when to start driving the network
%       DrivDur     - Scalar, indicates how long to drive the network
%       C           - N (destination) x N (origin) connectivity matrix with N being the number of
%                     neural masses
%       D           - N (destination) x N (origin) distance matrix with N being the number of
%                     neural masses
%       k           - Connectivity scaling constant
%       v           - Velocity of information transfer [m/s]
%       tMax        - Maximum simulation time [s]
%       dt          - Simulation stepsize [s]
%       initSampRem - Initial time interval to be cut-off [s]
%       nVar        - Variance of white noise used to drive neural masses
%       nMu         - Mean of white noise used to drive neural masses
%       rAvg        - if true, substract running average of input from
%                     other regions from input to neural masses
%       netInp      - scales input to [pyramidal, excitatory, inhibitory] 
%                     neurons from connectome
%       subInp      - scales input to [pyramidal, excitatory, inhibitory] 
%                     neurons from subcortical regions
%       sampling    - Inverse sampling frequency per sample
%       verbose     - If true, simulation progress will be displayed
%       JRParams    - Structure containing JansenRit parameters
%
%   Output:
% 
%       xColl       - 2-dimensional array with first dimension = number of neural
%                     masses & second dimension = timesteps. Values are
%                     post-synaptic potentials of pyramidal cells.
%       driverColl  - 2-dimensional array with first dimension = number of
%                     drivers & second dimension = timesteps. Values
%                     indicate voltage of drivers.

%% Initializations

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
%disp(['size of ring buffer: ' num2str(ringBufferSize)]);
ringBuffer = squeeze(StartStates(:,1,end-ringBufferSize+1:end))'; % timesteps x N

%ringBuffer = 2*JRParams.e0./(1 + exp(JRParams.r*(JRParams.u0 - ringBuffer))) - JRParams.S0;

%% Run simulation

t = 0;
x = squeeze(StartStates(:,:,end));
for i = 0:numIters-1
    t = t + dt;
    
    % realize global model delays D
    curPosInRingBuff = mod(i,ringBufferSize); 
    lookupInBuffer = mod(curPosInRingBuff-1-D,ringBufferSize)+1;
    linInd = sub2ind(size(ringBuffer),lookupInBuffer,repmat(1:size(ringBuffer,2),[N 1]));
    x1DelayedSpikeRate = ringBuffer(linInd);
    
    % extrinsic input updates
    if t > DrivStart && t < DrivStart + DrivDur
        drivPhase = sin(t*2*pi*DrivFreqs + DrivPO);
    else
        drivPhase = zeros(1,length(DrivPO));
    end
    driver = Drivers .* repmat(drivPhase,size(Drivers,1),1);
    driver = sum(driver,2);
    
    inp_net = sum(C .* x1DelayedSpikeRate, 2);
    inp_sub = nMu + nVar * randn(N, 1);
    
    % intrinsic state variable updates
    dx = odesys(x, driver, inp_sub, inp_net, JRParams);
    x = x + dt*dx;
    
    %update ring buffer
    
    if JRParams.use_out_psp
        x1spikeRateOut = x(:,7); % v_{o}
    elseif JRParams.use_sigm_as_out
        % use Sigm(v_{e} - v_{i})
        x1spikeRateOut = WTP(x(:,2) - x(:,3), JRParams.e0, JRParams.u0, JRParams.r);
    elseif JRParams.use_sigm_y0_as_out
        % use Sigm(v_{p})
        x1spikeRateOut = WTP(x(:,1), JRParams.e0, JRParams.u0, JRParams.r);
    else
        x1spikeRateOut = x(:,1); % v_{p}
    end
    
    ringBuffer(curPosInRingBuff+1,:) = x1spikeRateOut;
    
    %save every sampling steps:
    if mod(i,sampling)==0 && t > initSampRem
        xColl(:,i/sampling + 1) = x(:,1);
        driverColl(:,i/sampling + 1) = drivPhase;
    end
    
    %print progress:
    if verbose && mod(i,round(numIters/100))==0
        fprintf(['t = ' num2str(t) ' of ' num2str(tMax) '\n']);
    end
    
end

%% Store results

xColl = xColl(:,(initSampRem/(dt*sampling))+1 : end);
driverColl = driverColl(:,(initSampRem/(dt*sampling))+1 : end);

end

function dy = odesys(y, driver, inp_sub, inp_net, JRParams)

    % y(:,1) = v_{p}
    % y(:,2) = v_{e}
    % y(:,3) = v_{i}
    % y(:,4) = i_{p}
    % y(:,5) = i_{e}
    % y(:,6) = i_{i}
    
    % only if JRParams.use_out_psp then also:
    % y(:,7) = v_{o}
    % y(:,8) = i_{o}
    
    C1 = JRParams.Cpe;
    C2 = JRParams.Cep;
    C3 = JRParams.Cpi;
    C4 = JRParams.Cip;
    
    inpP =      WTP(y(:,2) - y(:,3) + driver, JRParams.e0, JRParams.u0, JRParams.r);
    inpE = C2 * WTP(C1 * y(:,1), JRParams.e0, JRParams.u0, JRParams.r) + inp_sub + inp_net;
    inpI = C4 * WTP(C3 * y(:,1), JRParams.e0, JRParams.u0, JRParams.r);
    
    dy = zeros(size(y));
    
    dy(:,1:3) = y(:,4:6);
    dy(:,4) = JRParams.Ke * inpP + JRParams.de1 * y(:,4) + JRParams.de2 * y(:,1);
    dy(:,5) = JRParams.Ke * inpE + JRParams.de1 * y(:,5) + JRParams.de2 * y(:,2);
    dy(:,6) = JRParams.Ki * inpI + JRParams.di1 * y(:,6) + JRParams.di2 * y(:,3);
    
    if JRParams.use_out_psp
        dy(:,7) = y(:,8);
        dy(:,8) = JRParams.Ko * inpP + JRParams.do1 * y(:,8) + JRParams.do2 * y(:,7);
    end
    
end
