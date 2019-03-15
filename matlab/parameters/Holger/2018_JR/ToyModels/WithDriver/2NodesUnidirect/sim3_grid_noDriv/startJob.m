clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000; % measured 3400 MB
params.Gridjob.wc_host = '!(*saturn*|*c01grid*|*ramsauer*|*kuma*|*taygete*|*helike*|*leda*|*orthosie*|*themisto*)';
params.Gridjob.jobname = 'Connectome';
params.Gridjob.continue = true;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.runOnlyJobIds = [];
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:29:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectomePaper.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectomePaper.k = 20; %global connection strength scaling
params.JansenRitConnectomePaper.v = 3; % velocity [m/s]
params.JansenRitConnectomePaper.tMax = 605; %max simulation time [seconds]
params.JansenRitConnectomePaper.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectomePaper.sampling = 2; % sampling every x steps
params.JansenRitConnectomePaper.noiseVar = 0; % variance of noise used to drive neural masses (22 was mentioned in msc thesis)
params.JansenRitConnectomePaper.noiseMu = [120, 320]; %num2cell(round(20:20:140)); %220; % mean of noise used to drive neural masses (220 was mentioned in msc thesis)
params.JansenRitConnectomePaper.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectomePaper.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.subInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.initSampRem = 124; %initial interval to remove [seconds]
params.JansenRitConnectomePaper.verbose = false; % if we want to print time steps to console

c_tmp = 135;
params.JansenRitConnectomePaper.cs = [c_tmp, c_tmp*0.8, c_tmp*0.25, c_tmp*0.25, 0]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)

params.JansenRitConnectomePaper.fTarget = 'drivFreq'; % [Hz]

params.JansenRitConnectomePaper.FCMeasure = {{'Coherence','phaseDiffHist'}}; % FC measure to use
params.JansenRitConnectomePaper.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectomePaper.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectomePaper.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectomePaper.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectomePaper.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectomePaper.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated

params.JansenRitConnectomePaper.u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V].. was 0 in master thesis
params.JansenRitConnectomePaper.He = 3.25e-3; % Average synaptic gain for excitatory synapses [V]
params.JansenRitConnectomePaper.Hi = 22e-3; % Average synaptic gain for inhibitory synapses [V]
params.JansenRitConnectomePaper.Te = 10e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
params.JansenRitConnectomePaper.Ti = 20e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
params.JansenRitConnectomePaper.e0 = 2.5; % determines maximum mean firing rate [1/s]
params.JansenRitConnectomePaper.r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]

params.JansenRitConnectomePaper.calcFC_nwin1 = false;
params.JansenRitConnectomePaper.subtract_S0 = false;
params.JansenRitConnectomePaper.use_moran = false;
params.JansenRitConnectomePaper.use_out_psp = false;
params.JansenRitConnectomePaper.use_sigm_as_out = false;
params.JansenRitConnectomePaper.use_inpP_as_out = true;
params.JansenRitConnectomePaper.use_sigm_y0_as_out = false;
params.JansenRitConnectomePaper.calcCohWithDriver = true;
params.JansenRitConnectomePaper.saveSpectrum = false;

% connectivity matrix
C = zeros(2);
C(1,2) = 0.1;
C(2,1) = 0;

% concat many to do multiple simulations simulatanously:
num_repeats = 16;
C_full = zeros(32, 32);
for k=1:num_repeats
    start_idx=2*(k-1)+1;
    C_full(start_idx:start_idx+1, start_idx:start_idx+1) = C;
end
C = C_full;

% delay matrix
d12 = 0:5:200; % distance between nodes at the ends of chain
D = zeros(2);

Delays = cell(1,length(d12));
for d=1:length(d12)
    D(1,2) = d12(d);
    D(2,1) = d12(d);
    
    D_full = zeros(32, 32);
    for k=1:num_repeats
        start_idx=2*(k-1)+1;
        D_full(start_idx:start_idx+1, start_idx:start_idx+1) = D;
    end
    Delays{d} = D_full;
end

params.JansenRitConnectomePaper.C = C; % connectivity matrix
params.JansenRitConnectomePaper.D = Delays; % distance matrix


%% set driver params
params.JansenRitConnectomePaper.drivPos = 1:32; % indices of network nodes to be driven
params.JansenRitConnectomePaper.drivScale = 0; % 0.25, 0.5 driver strength [mV]
phase_offsets = [0:0.0625:(1-0.0625)]*2*pi;
phase_offsets = cat(1, phase_offsets, zeros(size(phase_offsets)));
phase_offsets = phase_offsets(:);
params.JansenRitConnectomePaper.drivPO = phase_offsets'; % phase offset of drivers
params.JansenRitConnectomePaper.drivFreq = 11; % 10 and 11 frequency of driver [Hz]
params.JansenRitConnectomePaper.drivDur = 605; % duration of driver [s]

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);

%data = dataPaths();
%job = JansenRitConnectomePaper(fullfile(data.workdir, 'Holger/2018_JR/ToyModels/WithDriver/2NodesBidirect', my_foldername, 'temp_Connectome/jobDesc.mat'));
%job.finishJobs()