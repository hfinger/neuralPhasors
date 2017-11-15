classdef JansenRitConnectome < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = JansenRitConnectome(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.JansenRitConnectome.p = 2; %defines what kind of p norm to use for normalization of structural connectivity
      this.params.JansenRitConnectome.k = 8; %global connection strength scaling
      this.params.JansenRitConnectome.v = 4.5; % velocity [m/s]
      this.params.JansenRitConnectome.tMax = 305; %max simulation time [seconds]
      this.params.JansenRitConnectome.dt = 0.0001; % simulation step size [seconds]
      this.params.JansenRitConnectome.sampling = 10; % sampling every x steps
      this.params.JansenRitConnectome.noiseVar = 100; % variance of white noise used to drive neural masses
      this.params.JansenRitConnectome.noiseMu = 0; % mean of white noise used to drive neural masses
      this.params.JansenRitConnectome.runningAvg = false; % if true, substract running average from input to neural masses
      this.params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
      this.params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
      this.params.JansenRitConnectome.drivFreq = 29.5; %driving frequencies
      this.params.JansenRitConnectome.drivPos = [1 2]; % network node to drive
      this.params.JansenRitConnectome.drivRange = [1 1]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
      this.params.JansenRitConnectome.drivPO = [0]*pi; % phase offset of drivers
      this.params.JansenRitConnectome.drivScale = 0.; % amplitude/strength of drivers
      this.params.JansenRitConnectome.drivStart = 1; % timepoint at which to start driving [seconds]
      this.params.JansenRitConnectome.drivDur = 305; % driving duration [seconds]
      this.params.JansenRitConnectome.FCMeasure = 'Coherence';
      this.params.JansenRitConnectome.fullFC = true; %whether to store full functional connectivity matrix or only FC of driven region
      this.params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated FC to empirical FC
      this.params.JansenRitConnectome.nWindows = 1; %number of timewindows over which to calculate FC
      this.params.JansenRitConnectome.nBins = 80; % only important, if information theoretic FC measure is used. Number of bins to use to discretize phase signal
      this.params.JansenRitConnectome.storeY = false; % if true, raw signal will be stored on struct aswell
      this.params.JansenRitConnectome.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
      this.params.JansenRitConnectome.C = zeros(33); % Connectivity matrix
      this.params.JansenRitConnectome.D = ones(33); % Distance matrix
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});

    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%

      
      %%%% END EDIT HERE:                                        %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      % initialze simulation parameters
      sim.p = this.params.JansenRitConnectome.p;
      sim.k = this.params.JansenRitConnectome.k;
      sim.v = this.params.JansenRitConnectome.v;
      sim.tMax = this.params.JansenRitConnectome.tMax;
      sim.dt = this.params.JansenRitConnectome.dt;
      sim.sampling = this.params.JansenRitConnectome.sampling;
      sim.noiseVar = this.params.JansenRitConnectome.noiseVar;
      sim.noiseMu = this.params.JansenRitConnectome.noiseMu;
      sim.rAvg = this.params.JansenRitConnectome.runningAvg;
      sim.d = this.params.JansenRitConnectome.d;
      sim.verbose = this.params.JansenRitConnectome.verbose;
      sim.FCMeasure = this.params.JansenRitConnectome.FCMeasure;
      sim.fullFC = this.params.JansenRitConnectome.fullFC;
      sim.nWindows = this.params.JansenRitConnectome.nWindows;
      sim.nBins = this.params.JansenRitConnectome.nBins;
      sim.filterSig = this.params.JansenRitConnectome.filterSig;
      
      % Default Jansen Rit parameters
      He = 3.25e-3; % Average synaptic gain for excitatory synapses [V]
      Hi = 22e-3; % Average synaptic gain for inhibitory synapses [V]
      Te = 10e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
      Ti = 20e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
      cs = [128 128 64 64 0]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)
      JRParams.Cpe = cs(1); % Connection from pyramidal cells to excitatory interneurons
      JRParams.Cpi = cs(3); % Connection from pyramidal cells to inhibitory interneurons
      JRParams.Cep = cs(2); % Connection from excitatory interneurons to pyramidal cells
      JRParams.Cip = cs(4); % Connection from inhibitory interneurons to pyramidal cells
      JRParams.Cii = cs(5); % recurrent intrinsic connection from inhibitory interneurons to themselves (16)
      JRParams.u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V]
      JRParams.e0 = 2.5; % determines maximum mean firing rate [1/s]
      JRParams.r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]
      JRParams.Ke = He/Te; % excitatory synaptic constant
      JRParams.Ki = Hi/Ti; % inhibitory synaptic constant
      JRParams.de1 = -2/Te; % constant for first derivative of excitatory synapses
      JRParams.de2 = -1/Te^2; % constant for second derivative of excitatory synapses
      JRParams.di1 = -2/Ti; % constant for first derivative of inhibitory synapses
      JRParams.di2 = -1/Ti^2; % constant for second derivative of inhibitory synapses
      JRParams.Ka = 0; % spyke frequency adaptation rate constant [1/s] (1.9531)
      JRParams.S0 = WTP(0, JRParams.e0, JRParams.u0, JRParams.r); % mean firing rate constant offset
      
      % initialize driver parameters
      drivPos = this.params.JansenRitConnectome.drivPos;
      sim.drivRange = this.params.JansenRitConnectome.drivRange;
      sim.drivFreq = this.params.JansenRitConnectome.drivFreq;
      sim.drivPO = this.params.JansenRitConnectome.drivPO;
      sim.drivScale = this.params.JansenRitConnectome.drivScale;
      sim.drivStart = this.params.JansenRitConnectome.drivStart;
      sim.drivDur = this.params.JansenRitConnectome.drivDur;
      
      % bandpass filter parameters
      env.t_rm = 5; % time to start bandpassing the signal at [s]
      env.sigBandpass(1).Fst1 = 7.5; % end stop band [Hz] 
      env.sigBandpass(1).Fp1 = 8; % start pass band [Hz]
      env.sigBandpass(1).Fp2 = 11; % [Hz] end pass band
      env.sigBandpass(1).Fst2 = 11.5; % [Hz] start stop band
      env.sigBandpass(1).Ast1 = 9; % frequency attenuation in first stopband
      env.sigBandpass(1).Ap = 1; % passband ripples
      env.sigBandpass(1).Ast2 = 9; % frequency attenuation in second stopband
      
      % load C,D and anatomical indices
      if sum(sum(this.params.JansenRitConnectome.C)) == 0
          homotopeScaling = 0.1;
          [C,D] = getConnectome(1,sim.p,homotopeScaling);
      else
          C = this.params.JansenRitConnectome.C;
          D = this.params.JansenRitConnectome.D;
      end
      
      % get 3d coordinates for each node and for eeg electrodes
      %path_ROICoordinates = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e/fs_rois/ca03_fs_rois.mat';
      %path_electrodeCoordinates = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/ca_electrodeLocations/ca03_EEGLocations.txt';
      %roiData = load(path_ROICoordinates);
      %electrodeCoordinatesTable = readtable(path_electrodeCoordinates,'Delimiter',' ','ReadVariableNames',false);
      
      %roiCoordinates = roiData.fs_rois;
      %electrodeCoordinates = table2array(electrodeCoordinatesTable(:,2:end));
      
      % reorder roi coordinates after anatomical indices
      %roiCoordinates = roiCoordinates(resortIds,:);
      
      % create drivers:    
      %     find electrode position closest to node to be stimulated and extract
      %     stimulation strength for each node from gaussian centered around that
      %     electrode position
      drivers = zeros(size(C,1),size(drivPos',1));
      %ED = pdist2(electrodeCoordinates, roiCoordinates);
      for i = 1:size(drivers, 2)
      %    [~,pos] = min(ED(:,drivPos(i)));
      %    drivStrength = normpdf(ED(pos,:), 0, sim.drivRange(i));
      %    drivStrength = drivStrength - min(drivStrength);
      %    drivers(:,i) = (drivStrength/max(drivStrength)) * sim.drivScale;
           drivers(drivPos(i),i) = 1 * sim.drivScale; 
      end

      % run jansen rit simulation and store PSPs of pyramidal cells
      nStatesJR = 13;
      StartStates = zeros(size(C, 1), nStatesJR, 1/sim.dt);
      [ PSPs, Driver ] = runJansenRit( StartStates, drivers, sim.drivFreq, [0,sim.drivPO], sim.drivStart, sim.drivDur, C, D, sim.k, sim.v, sim.tMax, sim.dt, sim.d, sim.noiseVar, sim.noiseMu, sim.rAvg, sim.sampling, sim.verbose, JRParams);
      simResult.Y = vertcat(Driver,PSPs);
      sim.drivPos = drivPos + size(Driver, 1);
      simResult.sim = sim;
      clear PSPs
      clear Driver

      % extract peak frequencies of each node
      Fs = 1/(sim.dt*sim.sampling);
      Hs = spectrum.periodogram;
      freqs = zeros(size(simResult.Y,1),1);
      for n=1:size(simResult.Y,1)
        h = psd(Hs,simResult.Y(n,:),'Fs',Fs);
        f = h.Frequencies;
        d = h.Data;
        d(1:100) = 0;
        [~,target] = max(d);
        freqs(n) = f(target);
      end
      simResult.freqs = freqs;

      % apply bandpass filter
      [ simEval ] = calcEnvFC(env, sim.filterSig, simResult, 1, 0);

      % calculate FC for phase of signal
      if this.params.JansenRitConnectome.storeY
          Y_raw = simResult.Y;
      end
      winLength = floor((sim.tMax - sim.d) / sim.nWindows);
      simResult.Y = simEval.phaseBP{1,1};
      clear simEval
      WindowsStart = [1:winLength:sim.nWindows*winLength]; % starting points of time windows for which to evaluate coherence [seconds]
      WindowsEnd = [winLength:winLength:sim.nWindows*winLength]; % ending points of time windows for which to evaluate coherence [seconds]
      FCWindows = vertcat(WindowsStart, WindowsEnd);
      if strcmp(sim.FCMeasure,'Coherence')
          FC  = coherence( simResult, 1, round(FCWindows ./ (sim.dt * sim.sampling)), sim.fullFC);
      elseif strcmp(sim.FCMeasure,'ShannonEntropy')
          FC = shannonEntropy(simResult,sim.nBins,FCWindows,1,sim.fullFC);
      elseif strcmp(sim.FCMeasure,'MutualInformation')
          FC = mutualInformation(simResult,1,sim.fullFC,sim.nWindows,sim.nBins);
      elseif strcmp(sim.FCMeasure,'all')
          FC = coherence( simResult, 1, round(FCWindows ./ (sim.dt * sim.sampling)), sim.fullFC);
          simResult.SE = shannonEntropy(simResult,sim.nBins,FCWindows,1,sim.fullFC);
          %simResult.MI = mutualInformation(simResult,1,sim.fullFC,sim.nWindows,sim.nBins);
      else
          FC = zeros(66) - 66;
      end
      if ~this.params.JansenRitConnectome.storeY
          simResult = rmfield(simResult,'Y');
      else
          simResult.Y = Y_raw;
          clear Y_raw
      end
      
      % calculate FC between driven regions
      drivPosFC = zeros(length(sim.drivPos)-1,1);
      if size(FC,1) == length(sim.drivPos)
          for j=2:length(sim.drivPos)
              drivPosFC(j-1) = mean(FC(1,sim.drivPos(j),:),3);
          end
      else
          for j=2:length(sim.drivPos)
              drivPosFC(j-1) = FC(sim.drivPos(1),sim.drivPos(j));
          end
      end

      % calculate match with empirical FC data
      if this.params.JansenRitConnectome.corrSimFC
        % load resort IDs
        path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
        load(path_ResortIDs);
        resortIds = [resortIdsMarlene, resortIdsMarlene + 33];
        % load empirical FC
        path_empFCmat = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/eeg_20150114_controls_fs_funconn_lcmv_bponetrial_hilbert_3_30.mat';
        empFC_struct = load(path_empFCmat);

        % get N x N empirical FC matrix
        empFC_tmp = empFC_struct.coh_all;
        clear empFC_struct
        subjects = ones(size(empFC_tmp,1),1);
        subjects(5) = 0; subjects(14) = 0; subjects(16) = 0;
        empFC_tmp(6,2,:,:,:,:) = empFC_tmp(6,1,:,:,:,:);
        empFC_tmp(7,2,6,:,:,:) = empFC_tmp(7,1,6,:,:,:);
        empFC_tmp = squeeze(mean(empFC_tmp(subjects == 1,1,:,:,:,30), 1));
        empFC = squeeze(mean(empFC_tmp(5:6,:,:), 1));

        % reorder empFC according to marlenes IDs
        empFC = empFC(resortIds,:);
        empFC = empFC(:,resortIds);

        % calculate correlation between simulated and empirical FC
        Idx_mat = triu(ones(size(empFC)),1);
        corr_SimFC =  min(min(corrcoef(FC(Idx_mat == 1), empFC(Idx_mat == 1))));
        simResult.corr_SimFC = corr_SimFC;
      end

      simResult.drivPosFC = drivPosFC;
      simResult.FC = FC;
      
      % save results
      filename_JR = [this.params.Gridjob.jobname,num2str(this.currJobid)];
      save([this.resultpath,'/', filename_JR], 'simResult')
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

