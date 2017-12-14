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
      
      this.params.JansenRitConnectome.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
      this.params.JansenRitConnectome.k = 19; %global connection strength scaling
      this.params.JansenRitConnectome.v = 2.9; % velocity [m/s]
      this.params.JansenRitConnectome.tMax = 405; %max simulation time [seconds]
      this.params.JansenRitConnectome.dt = 0.0005; % simulation step size [seconds]
      this.params.JansenRitConnectome.sampling = 2; % sampling every x steps
      this.params.JansenRitConnectome.noiseVar = 0; % variance of white noise used to drive neural masses
      this.params.JansenRitConnectome.noiseMu = 200; % mean of white noise used to drive neural masses
      this.params.JansenRitConnectome.runningAvg = true; % if true, substract running average from input to neural masses
      this.params.JansenRitConnectome.netInp = [1 1 1]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
      this.params.JansenRitConnectome.subInp = [1 0 0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from subcortical regions
      this.params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
      this.params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
      this.params.JansenRitConnectome.drivFreq = 6; %driving frequencies
      this.params.JansenRitConnectome.drivPos = [1]; % network node to drive
      %this.params.JansenRitConnectome.drivRange = [1 1]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
      this.params.JansenRitConnectome.drivPO = [0]*pi; % phase offset of drivers
      this.params.JansenRitConnectome.drivScale = 0.; % amplitude/strength of drivers
      this.params.JansenRitConnectome.drivStart = 1; % timepoint at which to start driving [seconds]
      this.params.JansenRitConnectome.drivDur = 405; % driving duration [seconds]
      this.params.JansenRitConnectome.FCMeasure = {{'Coherence'}};
      this.params.JansenRitConnectome.fullFC = true; %whether to store full functional connectivity matrix or only FC of driven region
      this.params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated FC to empirical FC
      this.params.JansenRitConnectome.nWindows = 1; %number of timewindows over which to calculate FC
      this.params.JansenRitConnectome.nBins = 80; % only important, if information theoretic FC measure is used. Number of bins to use to discretize phase signal
      this.params.JansenRitConnectome.storeY = false; % if true, raw signal will be stored on struct aswell
      this.params.JansenRitConnectome.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
      this.params.JansenRitConnectome.C = []; % Connectivity matrix
      this.params.JansenRitConnectome.D = []; % Distance matrix
      this.params.JansenRitConnectome.nodeLesions = 0; % number of nodes to lesion along the along the connections between driven nodes
      
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
      sim.netInp = this.params.JansenRitConnectome.netInp;
      sim.subInp = this.params.JansenRitConnectome.subInp;
      sim.d = this.params.JansenRitConnectome.d;
      sim.verbose = this.params.JansenRitConnectome.verbose;
      sim.FCMeasure = this.params.JansenRitConnectome.FCMeasure;
      sim.fullFC = this.params.JansenRitConnectome.fullFC;
      sim.nWindows = this.params.JansenRitConnectome.nWindows;
      sim.nBins = this.params.JansenRitConnectome.nBins;
      sim.filterSig = this.params.JansenRitConnectome.filterSig;
      sim.nodeLesions = this.params.JansenRitConnectome.nodeLesions;
      sim.C = this.params.JansenRitConnectome.C;
      sim.D = this.params.JansenRitConnectome.D;
      
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
      %sim.drivRange = this.params.JansenRitConnectome.drivRange;
      sim.drivFreq = this.params.JansenRitConnectome.drivFreq;
      sim.drivScale = this.params.JansenRitConnectome.drivScale;
      sim.drivPO = this.params.JansenRitConnectome.drivPO;
      sim.drivStart = this.params.JansenRitConnectome.drivStart;
      sim.drivDur = this.params.JansenRitConnectome.drivDur;
      
      if length(sim.drivPO) < length(drivPos)
          sim.drivPO = [0,sim.drivPO];
      end
      
      % load C,D and anatomical indices
      if isempty(sim.C)
          homotopeScaling = 0.;
          singleHemisphere = true;
          [C,D] = getConnectome(1,sim.p,homotopeScaling,singleHemisphere);
      else
          C = sim.C;
          D = sim.D;
      end
      
      % lesion nodes in C
      for i=1:sim.nodeLesions
          [paths,~] = getDelayWeightedSWPs(C, 0.1, sim.drivPos, 5, false);
          if isempty(paths)
              break
          else
              path = paths{1,1};
          end
          if length(path) < 3
              C(path(1),path(2)) = 0;
              C(path(2),path(1)) = 0;
          else
              idx = randi([2,length(path(2:end-1))]);
              C(idx,:) = 0;
              C(:,idx) = 0;
          end
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
           drivers(drivPos(i),i) = sim.drivScale; 
      end

      % run jansen rit simulation and store PSPs of pyramidal cells
      nStatesJR = 13;
      StartStates = zeros(size(C, 1), nStatesJR, 1/sim.dt);
      [ PSPs, Driver ] = runJansenRit( StartStates, drivers, sim.drivFreq, sim.drivPO, sim.drivStart, sim.drivDur, C, D, sim.k, sim.v, sim.tMax, sim.dt, sim.d, sim.noiseVar, sim.noiseMu, sim.rAvg, sim.netInp, sim.subInp, sim.sampling, sim.verbose, JRParams);
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
      fMean = mean(freqs);
      
      % apply bandpass filter
      if sim.filterSig
          
          % bandpass filter parameters
          sigBandpass(1).Fst1 = sim.drivFreq-1.; % end stop band [Hz] 
          sigBandpass(1).Fp1 = sim.drivFreq-0.5; % start pass band [Hz]
          sigBandpass(1).Fp2 = sim.drivFreq+0.5; % [Hz] end pass band
          sigBandpass(1).Fst2 = sim.drivFreq+1.; % [Hz] start stop band
          sigBandpass(1).Ast1 = sim.drivFreq; % frequency attenuation in first stopband
          sigBandpass(1).Ap = 1; % passband ripples
          sigBandpass(1).Ast2 = sim.drivFreq; % frequency attenuation in second stopband
          
          % filter signal
          Yfiltered = filterSig(simResult.Y,Fs,1,0,sigBandpass);
          
      else
          Yfiltered = simResult.Y;
      end
      
      if this.params.JansenRitConnectome.storeY
          Y_raw = simResult.Y;
      end
      simResult.Y = Yfiltered;
      clear Yfiltered
      
      % calculate FC for phase of signal      
      winLength = floor((sim.tMax - sim.d) / sim.nWindows);
      WindowsStart = [1:winLength:sim.nWindows*winLength]; % starting points of time windows for which to evaluate coherence [seconds]
      WindowsEnd = [winLength:winLength:sim.nWindows*winLength]; % ending points of time windows for which to evaluate coherence [seconds]
      FCWindows = vertcat(WindowsStart, WindowsEnd);
      FC = getFC(simResult,sim.FCMeasure,FCWindows,1);

      if ~this.params.JansenRitConnectome.storeY
          simResult = rmfield(simResult,'Y');
      else
          simResult.Y = Y_raw;
          clear Y_raw
      end

      % calculate match with empirical FC data
      if this.params.JansenRitConnectome.corrSimFC
          
        % load empirical FC
        fTarget = round(fMean);
        if fTarget < 3
            fTarget = 3;
        elseif fTarget > 30
            fTarget = 30;
        end
        empFC = getEmpFC(1,fTarget,1);

        % calculate correlation between simulated and empirical FC
        corr_SimFC = cell(1,length(FC));
        for i=1:length(FC)
            Idx_mat = triu(ones(size(empFC)),1);
            FC_tmp = FC{1,i}(length(drivPos)+1:end,length(drivPos)+1:end);
            corr_SimFC{i} =  min(min(corrcoef(abs(FC_tmp(Idx_mat == 1)), empFC(Idx_mat == 1))));
        end
        simResult.corr_SimFC = corr_SimFC;
        
      end

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

