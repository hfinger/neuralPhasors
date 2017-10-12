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
      this.params.JansenRitConnectome.tMax = 125; %max simulation time [seconds]
      this.params.JansenRitConnectome.dt = 0.0001; % simulation step size [seconds]
      this.params.JansenRitConnectome.sampling = 10; % sampling every x steps
      this.params.JansenRitConnectome.snr = 0; % amount of noise
      this.params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
      this.params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
      this.params.JansenRitConnectome.drivFreq = 29.5; %driving frequencies
      this.params.JansenRitConnectome.drivPos = [4,31]; % network node to drive
      this.params.JansenRitConnectome.drivRange = [20,15]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
      this.params.JansenRitConnectome.drivPO = [0:0.25:1.75]*pi; % phase offset of drivers
      this.params.JansenRitConnectome.drivScale = 0.02; % amplitude/strength of drivers
      this.params.JansenRitConnectome.drivStart = 1; % timepoint at which to start driving [seconds]
      this.params.JansenRitConnectome.drivDur = 125; % driving duration [seconds]
      this.params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
      this.params.JansenRitConnectome.fullCoherence = false; %whether to store full coherence matrix or only coherence of driven region
      this.params.JansenRitConnectome.nWindows = 4; %number of timewindows over which to calculate mean coherence 
      this.params.JansenRitConnectome.crossCorr = false; % if true, calculate cross-correlation between driven regions
      
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
      sim.snr = this.params.JansenRitConnectome.snr;
      sim.d = this.params.JansenRitConnectome.d;
      sim.verbose = this.params.JansenRitConnectome.verbose;
      sim.fullCoherence = this.params.JansenRitConnectome.fullCoherence;
      sim.nWindows = this.params.JansenRitConnectome.nWindows;
      
      % Default Jansen Rit parameters
      He = 15e-3; % Average synaptic gain for excitatory synapses [V]
      Hi = 22e-3; % Average synaptic gain for inhibitory synapses [V]
      Te = 3.5e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
      Ti = 10e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
      cs = [128 128 64 64 16]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)
      JRParams.Cpe = cs(1); % Connection from pyramidal cells to excitatory interneurons
      JRParams.Cpi = cs(3); % Connection from pyramidal cells to inhibitory interneurons
      JRParams.Cep = cs(2); % Connection from excitatory interneurons to pyramidal cells
      JRParams.Cip = cs(4); % Connection from inhibitory interneurons to pyramidal cells
      JRParams.Cii = cs(5); % recurrent intrinsic connection from inhibitory interneurons to themselves
      JRParams.u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V]
      JRParams.e0 = 2.5; % determines maximum mean firing rate [1/s]
      JRParams.r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]
      JRParams.Ke = He/Te; % excitatory synaptic constant
      JRParams.Ki = Hi/Ti; % inhibitory synaptic constant
      JRParams.de1 = -2/Te; % constant for first derivative of excitatory synapses
      JRParams.de2 = -1/Te^2; % constant for second derivative of excitatory synapses
      JRParams.di1 = -2/Ti; % constant for first derivative of inhibitory synapses
      JRParams.di2 = -1/Ti^2; % constant for second derivative of inhibitory synapses
      JRParams.Ka = 1.9531; % spyke frequency adaptation rate constant [1/s]
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
      env.sigBandpass(1).Fst1 = 25.5; % end stop band [Hz] 
      env.sigBandpass(1).Fp1 = 28; % start pass band [Hz]
      env.sigBandpass(1).Fp2 = 32; % [Hz] end pass band
      env.sigBandpass(1).Fst2 = 34.5; % [Hz] start stop band
      env.sigBandpass(1).Ast1 = 30; % frequency attenuation in first stopband
      env.sigBandpass(1).Ap = 1; % passband ripples
      env.sigBandpass(1).Ast2 = 30; % frequency attenuation in second stopband
      
      %env.sigBandpass(2).Fst1 = 34; % end stop band [Hz] 
      %env.sigBandpass(2).Fp1 = 37; % start pass band [Hz]
      %env.sigBandpass(2).Fp2 = 43; % [Hz] end pass band
      %env.sigBandpass(2).Fst2 = 46; % [Hz] start stop band
      %env.sigBandpass(2).Ast1 = 40; % frequency attenuation in first stopband
      %env.sigBandpass(2).Ap = 1; % passband ripples
      %env.sigBandpass(2).Ast2 = 40; % frequency attenuation in second stopband
      
      %env.sigBandpass(3).Fst1 = 54; % end stop band [Hz] 
      %env.sigBandpass(3).Fp1 = 57; % start pass band [Hz]
      %env.sigBandpass(3).Fp2 = 63; % [Hz] end pass band
      %env.sigBandpass(3).Fst2 = 66; % [Hz] start stop band
      %env.sigBandpass(3).Ast1 = 60; % frequency attenuation in first stopband
      %env.sigBandpass(3).Ap = 1; % passband ripples
      %env.sigBandpass(3).Ast2 = 60; % frequency attenuation in second stopband
      
      % load C,D and anatomical indices
      path_SCmat = '/net/store/nbp/projects/phasesim/databases/avg_SC.mat';
      path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
      load(path_SCmat);
      load(path_ResortIDs);
      resortIds = [resortIdsMarlene, resortIdsMarlene + 33];
      
      % get 3d coordinates for each node and for eeg electrodes
      %path_ROICoordinates = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e/fs_rois/ca03_fs_rois.mat';
      %path_electrodeCoordinates = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/ca_electrodeLocations/ca03_EEGLocations.txt';
      %roiData = load(path_ROICoordinates);
      %electrodeCoordinatesTable = readtable(path_electrodeCoordinates,'Delimiter',' ','ReadVariableNames',false);
      
      %roiCoordinates = roiData.fs_rois;
      %electrodeCoordinates = table2array(electrodeCoordinatesTable(:,2:end));
      
      % reorder C, D and roi coordinates after anatomical indices
      C = C(resortIds,:);
      C = C(:,resortIds);
      D = D(resortIds,:);
      D = D(:,resortIds);
      %roiCoordinates = roiCoordinates(resortIds,:);
      
      % add homotopic connections and normalize rows to 1
      C = C + 0.1 * diag(ones(size(C,1)/2,1),size(C,1)/2) + 0.1 * diag(ones(size(C,1)/2,1),-size(C,1)/2);
      C = bsxfun(@rdivide,C,sum(C.^sim.p,2).^(1/sim.p));
      
      % calculate betweenness centrality of each node based on C
      [~, BC] = edge_betweenness_wei(1./C);
      
      % store betweenness centrality of stimulated nodes
      sim.BC = BC(drivPos);
      
      % run for all driver phase offsets
      for p=1:length(sim.drivPO)
          
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
          [ PSPs, Driver ] = runJansenRit( StartStates, drivers, sim.drivFreq, [0,sim.drivPO(p)], sim.drivStart, sim.drivDur, C, D, sim.k, sim.v, sim.tMax, sim.dt, sim.d, sim.snr, sim.sampling, sim.verbose, JRParams);
          simResult.Y = vertcat(Driver,PSPs);
          sim.drivPos = drivPos + size(Driver, 1);
          simResult.sim = sim;

          % extract peak frequencies of each node
          %Fs = 1/(sim.dt*sim.sampling);
          %Hs = spectrum.periodogram;
          %freqs = zeros(size(simResult.Y,1),1);
          %for n=1:size(simResult.Y,1)
          %  h = psd(Hs,simResult.Y(n,:),'Fs',Fs);
          %  f = h.Frequencies;
          %  d = h.Data;
          %  d(1:100) = 0;
          %  [~,target] = max(d);
          %  freqs(n) = f(target);
          %end
          %simResult.freqs = freqs;

          % apply bandpass filter
          [ simEval ] = calcEnvFC(env, simResult, 1, 0);

          % calculate coherence for each bandpass filtered signal
          Y_raw = simResult.Y;
          winLength = floor((sim.tMax - sim.d) / sim.nWindows);
          for i=1:length(simEval.phaseBP)
              simResult.Y = simEval.phaseBP{1,i};
              CohWindowsStart = [1:winLength:sim.nWindows*winLength]; % starting points of time windows for which to evaluate coherence [seconds]
              CohWindowsEnd = [winLength:winLength:sim.nWindows*winLength]; % ending points of time windows for which to evaluate coherence [seconds]
              cohWindows = vertcat(CohWindowsStart, CohWindowsEnd);
              [ Coherence_tmp ] = coherence( simResult, 1, 1, round(cohWindows ./ (sim.dt * sim.sampling)), sim.fullCoherence);
              Coherence{p,i} = Coherence_tmp;
              drivPosCoh_tmp = zeros(length(sim.drivPos)-1,1);
              if size(Coherence_tmp,1) == length(sim.drivPos)
                  for j=2:length(sim.drivPos)
                      drivPosCoh_tmp(j-1) = mean(Coherence_tmp(1,sim.drivPos(j),:),3);
                  end
              else
                  for j=2:length(sim.drivPos)
                      drivPosCoh_tmp(j-1) = Coherence_tmp(sim.drivPos(1),sim.drivPos(j));
                  end
              end
          end
          %simResult.Y = Y_raw;
          simResult = rmfield(simResult,'Y');
          drivPosCoh{p} = drivPosCoh_tmp;
          
          % calculate cross-correlation between driven regions
          if this.params.JansenRitConnectome.crossCorr
              maxLag = size(Y_raw,2)-round(size(Y_raw,2)/2);
              for i=2:size(drivPos',1)
                  crossCorr{p,i} = xcorr(Y_raw(drivPos(1),:),Y_raw(drivPos(i),:),maxLag,'unbiased');
              end
          end
          
          % calculate match with empirical FC data
          if this.params.JansenRitConnectome.corrSimFC
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
            for i=1:length(Coherence)
                simFC = Coherence{1,i};
                corr_SimFC{i} =  min(min(corrcoef(simFC(Idx_mat == 1), empFC(Idx_mat == 1))));
            end
            simResult.corr_SimFC = corr_SimFC;
          end
          
      end
      if this.params.JansenRitConnectome.crossCorr
          crossCorr = squeeze(cell2mat(crossCorr));
          simResult.crossCorr = crossCorr;
      end
      simResult.drivPosCoh = drivPosCoh;
      simResult.Coherence = Coherence;
      
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

