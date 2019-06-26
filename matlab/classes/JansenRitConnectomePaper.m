classdef JansenRitConnectomePaper < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = JansenRitConnectomePaper(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.JansenRitConnectomePaper.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
      this.params.JansenRitConnectomePaper.k = 19; %global connection strength scaling
      this.params.JansenRitConnectomePaper.v = 2.9; % velocity [m/s]
      this.params.JansenRitConnectomePaper.tMax = 405; %max simulation time [seconds]
      this.params.JansenRitConnectomePaper.dt = 0.0005; % simulation step size [seconds]
      this.params.JansenRitConnectomePaper.sampling = 2; % sampling every x steps
      this.params.JansenRitConnectomePaper.sampleNoiseEvery = 1; % sampling of noise every x steps
      this.params.JansenRitConnectomePaper.noiseVar = 0; % variance of white noise used to drive neural masses
      this.params.JansenRitConnectomePaper.noiseMu = 200; % mean of white noise used to drive neural masses
      this.params.JansenRitConnectomePaper.runningAvg = true; % if true, substract running average from input to neural masses
      this.params.JansenRitConnectomePaper.netInp = [1 1 1]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
      this.params.JansenRitConnectomePaper.subInp = [1 0 0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from subcortical regions
      this.params.JansenRitConnectomePaper.initSampRem = 5; %initial interval to remove [seconds]
      this.params.JansenRitConnectomePaper.verbose = false; % if we want to print time steps to console
      
      this.params.JansenRitConnectomePaper.fTarget = 9; % [Hz]

      c_tmp = 135;
      this.params.JansenRitConnectomePaper.cs = [c_tmp, c_tmp*0.8, c_tmp*0.25, c_tmp*0.25, 0]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)
      
      this.params.JansenRitConnectomePaper.drivFreq = 6; %driving frequencies
      this.params.JansenRitConnectomePaper.drivPosVarMatrix = []; % network node to drive
      this.params.JansenRitConnectomePaper.drivPos = [1]; % network node to drive
      %this.params.JansenRitConnectomePaper.drivRange = [1 1]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
      this.params.JansenRitConnectomePaper.drivPO = [0]*pi; % phase offset of drivers
      this.params.JansenRitConnectomePaper.drivScale = 0.; % amplitude/strength of drivers [mV]
      this.params.JansenRitConnectomePaper.drivStart = 1; % timepoint at which to start driving [seconds]
      this.params.JansenRitConnectomePaper.drivDur = 405; % driving duration [seconds]
      
      this.params.JansenRitConnectomePaper.FCMeasure = {{'Coherence'}};
      this.params.JansenRitConnectomePaper.fullFC = true; %whether to store full functional connectivity matrix or only FC of driven region
      this.params.JansenRitConnectomePaper.corrSimFC = false; % if true, compare simulated FC to empirical FC
      this.params.JansenRitConnectomePaper.nWindows = 1; %number of timewindows over which to calculate FC
      this.params.JansenRitConnectomePaper.nBins = 80; % only important, if information theoretic FC measure is used. Number of bins to use to discretize phase signal
      this.params.JansenRitConnectomePaper.storeY = false; % if true, raw signal will be stored on struct aswell
      this.params.JansenRitConnectomePaper.storeMeanAndStdOfY = false; % if true, store mean and std of raw signals.
      this.params.JansenRitConnectomePaper.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
      this.params.JansenRitConnectomePaper.C = []; % Connectivity matrix
      this.params.JansenRitConnectomePaper.D = []; % Distance matrix
      this.params.JansenRitConnectomePaper.nodeLesions = 0; % number of nodes to lesion along the along the connections between driven nodes

      this.params.JansenRitConnectomePaper.washout = 1;
      
      this.params.JansenRitConnectomePaper.u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V].. was 0 in master thesis
      this.params.JansenRitConnectomePaper.He = 3.25e-3; % Average synaptic gain for excitatory synapses [V]
      this.params.JansenRitConnectomePaper.Hi = 22e-3; % Average synaptic gain for inhibitory synapses [V]
      this.params.JansenRitConnectomePaper.Te = 10e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
      this.params.JansenRitConnectomePaper.Ti = 20e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
      
      this.params.JansenRitConnectomePaper.e0 = 2.5; % determines maximum mean firing rate [1/s]
      this.params.JansenRitConnectomePaper.r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]
      
      this.params.JansenRitConnectomePaper.calcFC_nwin1 = false;
      this.params.JansenRitConnectomePaper.subtract_S0 = false;
      this.params.JansenRitConnectomePaper.use_moran = true;
      this.params.JansenRitConnectomePaper.use_out_psp = false;
      this.params.JansenRitConnectomePaper.use_sigm_as_out = false;
      this.params.JansenRitConnectomePaper.use_inpP_as_out = false;
      this.params.JansenRitConnectomePaper.use_sigm_y0_as_out = false;
      this.params.JansenRitConnectomePaper.use_rk4 = false;
      this.params.JansenRitConnectomePaper.calcCohWithDriver = false;
      this.params.JansenRitConnectomePaper.calcKuramotoOrderParam = false;
      this.params.JansenRitConnectomePaper.saveSpectrum = false;
      
      this.params.JansenRitConnectomePaper.xCollExtended = false;
      this.params.JansenRitConnectomePaper.xCollExtendedIdxs = 1:4;

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
      
      filename_JR = [this.params.Gridjob.jobname num2str(this.currJobid) '.mat'];
      
      if exist(fullfile(this.workpath, filename_JR), 'file')
          disp('job output already exists')
          return
      end
      
      % initialze simulation parameters
      sim.p = this.params.JansenRitConnectomePaper.p;
      sim.k = this.params.JansenRitConnectomePaper.k;
      sim.v = this.params.JansenRitConnectomePaper.v;
      sim.tMax = this.params.JansenRitConnectomePaper.tMax;
      sim.dt = this.params.JansenRitConnectomePaper.dt;
      sim.sampling = this.params.JansenRitConnectomePaper.sampling;
      sim.sampleNoiseEvery = this.params.JansenRitConnectomePaper.sampleNoiseEvery;
      sim.noiseVar = this.params.JansenRitConnectomePaper.noiseVar;
      sim.noiseMu = this.params.JansenRitConnectomePaper.noiseMu;
      sim.rAvg = this.params.JansenRitConnectomePaper.runningAvg;
      sim.netInp = this.params.JansenRitConnectomePaper.netInp;
      sim.subInp = this.params.JansenRitConnectomePaper.subInp;
      sim.initSampRem = this.params.JansenRitConnectomePaper.initSampRem;
      sim.verbose = this.params.JansenRitConnectomePaper.verbose;
      sim.FCMeasure = this.params.JansenRitConnectomePaper.FCMeasure;
      sim.fullFC = this.params.JansenRitConnectomePaper.fullFC;
      sim.nWindows = this.params.JansenRitConnectomePaper.nWindows;
      sim.nBins = this.params.JansenRitConnectomePaper.nBins;
      sim.filterSig = this.params.JansenRitConnectomePaper.filterSig;
      sim.nodeLesions = this.params.JansenRitConnectomePaper.nodeLesions;
      sim.C = this.params.JansenRitConnectomePaper.C;
      sim.D = this.params.JansenRitConnectomePaper.D;
      
      % Default Jansen Rit parameters
      He = this.params.JansenRitConnectomePaper.He; %3.25e-3; % Average synaptic gain for excitatory synapses [V]
      Hi = this.params.JansenRitConnectomePaper.Hi; %22e-3; % Average synaptic gain for inhibitory synapses [V]
      Te = this.params.JansenRitConnectomePaper.Te; %10e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
      Ti = this.params.JansenRitConnectomePaper.Ti; %20e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
      
      cs = this.params.JansenRitConnectomePaper.cs;
      
      JRParams.Cpe = cs(1); % Connection from pyramidal cells to excitatory interneurons
      JRParams.Cpi = cs(3); % Connection from pyramidal cells to inhibitory interneurons
      JRParams.Cep = cs(2); % Connection from excitatory interneurons to pyramidal cells
      JRParams.Cip = cs(4); % Connection from inhibitory interneurons to pyramidal cells
      JRParams.Cii = cs(5); % recurrent intrinsic connection from inhibitory interneurons to themselves (16)
      JRParams.u0 = this.params.JansenRitConnectomePaper.u0; %0;%6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V]
      JRParams.e0 = this.params.JansenRitConnectomePaper.e0; %2.5; % determines maximum mean firing rate [1/s]
      JRParams.r = this.params.JansenRitConnectomePaper.r; %560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]
      JRParams.Ke = He/Te; % excitatory synaptic constant
      JRParams.Ki = Hi/Ti; % inhibitory synaptic constant
      JRParams.de1 = -2/Te; % constant for first derivative of excitatory synapses
      JRParams.de2 = -1/Te^2; % constant for second derivative of excitatory synapses
      JRParams.di1 = -2/Ti; % constant for first derivative of inhibitory synapses
      JRParams.di2 = -1/Ti^2; % constant for second derivative of inhibitory synapses
      JRParams.Ka = 0; % spyke frequency adaptation rate constant [1/s] (1.9531)
      JRParams.S0 = WTP(0, JRParams.e0, JRParams.u0, JRParams.r); % mean firing rate constant offset
      JRParams.use_out_psp = this.params.JansenRitConnectomePaper.use_out_psp;
      JRParams.use_sigm_as_out = this.params.JansenRitConnectomePaper.use_sigm_as_out;
      JRParams.use_inpP_as_out = this.params.JansenRitConnectomePaper.use_inpP_as_out;
      JRParams.use_sigm_y0_as_out = this.params.JansenRitConnectomePaper.use_sigm_y0_as_out;
      
      JRParams.use_rk4 = this.params.JansenRitConnectomePaper.use_rk4;
      JRParams.sampleNoiseEvery = this.params.JansenRitConnectomePaper.sampleNoiseEvery;
      
      JRParams.xCollExtended = this.params.JansenRitConnectomePaper.xCollExtended;
      JRParams.xCollExtendedIdxs = this.params.JansenRitConnectomePaper.xCollExtendedIdxs;
      
      if ~this.params.JansenRitConnectomePaper.subtract_S0
          JRParams.S0 = 0;
      end
      
      if JRParams.use_out_psp
        To = Te * 3;
        JRParams.Ko = He/To; % excitatory synaptic constant
        JRParams.do1 = -2/To; % constant for first derivative of excitatory synapses
        JRParams.do2 = -1/To^2; % constant for second derivative of excitatory synapses
      end
      
      numRepeats = 1;
      if ~isempty(this.params.JansenRitConnectomePaper.drivPosVarMatrix)
        numRepeats = size(this.params.JansenRitConnectomePaper.drivPosVarMatrix, 1);
        simResultAllTmp = cell(1,numRepeats);
      end
      
      for repeatIdx = 1:numRepeats

        if ~isempty(this.params.JansenRitConnectomePaper.drivPosVarMatrix)
          drivPos = this.params.JansenRitConnectomePaper.drivPosVarMatrix(repeatIdx, :);
        else
          drivPos = this.params.JansenRitConnectomePaper.drivPos;
        end
        
        %sim.drivRange = this.params.JansenRitConnectomePaper.drivRange;
        sim.drivFreq = this.params.JansenRitConnectomePaper.drivFreq;
        sim.drivScale = this.params.JansenRitConnectomePaper.drivScale;
        sim.drivPO = this.params.JansenRitConnectomePaper.drivPO;
        sim.drivStart = this.params.JansenRitConnectomePaper.drivStart;
        sim.drivDur = this.params.JansenRitConnectomePaper.drivDur;

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
             drivers(drivPos(i),i) = sim.drivScale * 1e-3 / 2;  % convert from mV to V and divide by 2 to use peak-to-peak amplitude for sin
        end

        % run jansen rit simulation and store PSPs of pyramidal cells
        if this.params.JansenRitConnectomePaper.use_moran
            nStatesJR = 13;
            StartStates = zeros(size(C, 1), nStatesJR, 1/sim.dt);
            [ PSPs, Driver, xCollExtended ] = runJansenRit( StartStates, drivers, sim.drivFreq, sim.drivPO, sim.drivStart+2*pi*rand(1), sim.drivDur, C, D, sim.k, sim.v, sim.tMax, sim.dt, sim.initSampRem, sim.noiseVar, sim.noiseMu, sim.rAvg, sim.netInp, sim.subInp, sim.sampling, sim.verbose, JRParams);
        else

            if JRParams.use_out_psp
              nStatesJR = 8;
            else
              nStatesJR = 6;
            end

            StartStates = zeros(size(C, 1), nStatesJR, round(1/sim.dt));
            [ PSPs, Driver, xCollExtended ] = runJansenRitOriginal( StartStates, drivers, sim.drivFreq, sim.drivPO, sim.drivStart+2*pi*rand(1), sim.drivDur, C, D, sim.k, sim.v, sim.tMax, sim.dt, sim.initSampRem, sim.noiseVar, sim.noiseMu, sim.rAvg, sim.netInp, sim.subInp, sim.sampling, sim.verbose, JRParams);
        end
        simResult.Y = vertcat(Driver,PSPs);
        sim.drivPos = drivPos + size(Driver, 1);
        simResult.sim = sim;
        clear PSPs
        clear Driver

        if this.params.JansenRitConnectomePaper.storeMeanAndStdOfY
            simResult.Y_mean = mean(simResult.Y, 2);
            simResult.Y_std = std(simResult.Y, [], 2);
        end
        
        % extract peak frequencies of each node
        Fs = 1/(sim.dt*sim.sampling);
        Hs = spectrum.periodogram;
        freqs = zeros(size(simResult.Y,1),1);
        if this.params.JansenRitConnectomePaper.saveSpectrum
            simResult.spectrumPower = cell(size(simResult.Y,1),1);
        end
        for n=1:size(simResult.Y,1)
          h = psd(Hs,simResult.Y(n,:),'Fs',Fs);
          f = h.Frequencies;
          d = h.Data;
          d(1:100) = 0;
          [~,target] = max(d);
          freqs(n) = f(target);
          if this.params.JansenRitConnectomePaper.saveSpectrum
              saveTo = find(f>30,1);
              simResult.spectrumFreq = f(1:saveTo);
              simResult.spectrumPower{n} = d(1:saveTo);
          end
        end
        simResult.freqs = freqs;
        fMean = mean(freqs(length(drivPos)+1:end));

        if strcmp(this.params.JansenRitConnectomePaper.fTarget, 'fMean')
            fTarget = fMean;
        elseif strcmp(this.params.JansenRitConnectomePaper.fTarget, 'drivFreq')
            fTarget = this.params.JansenRitConnectomePaper.drivFreq;
        else
            fTarget = this.params.JansenRitConnectomePaper.fTarget;
        end

        if fTarget < 3
            fTarget = 3;
        elseif fTarget > 30
            fTarget = 30;
        end


        % apply bandpass filter
        if sim.filterSig

            % bandpass filter parameters
            sigBandpass(1).Fst1 = fTarget-1.; % end stop band [Hz] 
            sigBandpass(1).Fp1 = fTarget-0.5; % start pass band [Hz]
            sigBandpass(1).Fp2 = fTarget+0.5; % [Hz] end pass band
            sigBandpass(1).Fst2 = fTarget+1.; % [Hz] start stop band
            sigBandpass(1).Ast1 = fTarget; % frequency attenuation in first stopband
            sigBandpass(1).Ap = 1; % passband ripples
            sigBandpass(1).Ast2 = fTarget; % frequency attenuation in second stopband

            % filter signal
            Yfiltered = filterSig(simResult.Y,Fs,1,0,sigBandpass);

        else
            Yfiltered = simResult.Y;
        end

        washout_in_sec = this.params.JansenRitConnectomePaper.washout * simResult.sim.dt * simResult.sim.sampling;

        if this.params.JansenRitConnectomePaper.storeY
            Y_raw = simResult.Y;
        end
        
        simResult.Y = Yfiltered;
        clear Yfiltered
        % calculate FC for phase of signal
        if this.params.JansenRitConnectomePaper.calcFC_nwin1 % for debugging...
            nWindows = 1;
            winLength = floor((sim.tMax - sim.initSampRem - washout_in_sec) / nWindows);
            WindowsStart = [1:winLength:nWindows*winLength]; % starting points of time windows for which to evaluate coherence [seconds]
            WindowsEnd = [winLength:winLength:nWindows*winLength]; % ending points of time windows for which to evaluate coherence [seconds]
            FCWindows = vertcat(WindowsStart, WindowsEnd);
            FC_nwin1 = getFC(simResult,sim.FCMeasure,FCWindows,this.params.JansenRitConnectomePaper.washout);
            simResult.FC_nwin1 = FC_nwin1;
        end

        nWindows = sim.nWindows;
        winLength = floor((sim.tMax - sim.initSampRem - washout_in_sec) / nWindows);
        WindowsStart = [1:winLength:nWindows*winLength]; % starting points of time windows for which to evaluate coherence [seconds]
        WindowsEnd = [winLength:winLength:nWindows*winLength]; % ending points of time windows for which to evaluate coherence [seconds]
        FCWindows = vertcat(WindowsStart, WindowsEnd);
        FC = getFC(simResult,sim.FCMeasure,FCWindows,this.params.JansenRitConnectomePaper.washout);

        % calculate kuramoto order for revision:
        if this.params.JansenRitConnectomePaper.calcKuramotoOrderParam

            % get signal
            sig = simResult.Y(length(drivPos)+1:end,2000:end-2000); % leave out 2000 samples due to filter edge effects.

            % get analytic signal
            sigHilbert = zeros(size(sig));
            for n=1:size(sigHilbert,1)
                sigHilbert(n,:) = hilbert(sig(n,:));
            end

            % get phase from analytic signal
            sigPhase = angle(sigHilbert);

            kuramotoOrderParam = abs(mean( exp(1i * sigPhase), 1 ));
            meanKuramotoOrderParam = mean(kuramotoOrderParam);

            simResult.meanKuramotoOrderParam = meanKuramotoOrderParam;
        end

        if ~this.params.JansenRitConnectomePaper.storeY
            simResult = rmfield(simResult,'Y');
        else
            simResult.Y = Y_raw;
            clear Y_raw
        end

        % calculate match with empirical FC data
        if this.params.JansenRitConnectomePaper.corrSimFC

          % load empirical FC
          empFC = getEmpFC(1,round(fTarget),1);
          simResult.empFC = empFC;

          % calculate correlation between simulated and empirical FC
          corr_SimFC = cell(1,length(FC));
          for i=1:length(FC)
              Idx_mat = triu(ones(size(empFC)),1);
              FC_tmp = FC{1,i}(length(drivPos)+1:end,length(drivPos)+1:end);
              corr_SimFC{i} =  min(min(corrcoef(abs(FC_tmp(Idx_mat == 1)), empFC(Idx_mat == 1))));
              disp(['corr_SimFC{i} = ' num2str(corr_SimFC{i})])
          end
          simResult.corr_SimFC = corr_SimFC;

        end

        simResult.FC = FC;

        if this.params.JansenRitConnectomePaper.calcCohWithDriver
            coh_of_roi_with_driver = cell(1,length(FC));
            for i=1:length(FC)
                FC_tmp = FC{1,i}(length(drivPos)+1:end,1);
                coh_of_roi_with_driver{i} = zeros(1, length(drivPos));
                for j=1:length(drivPos)
                  coh_of_roi_with_driver{i}(j) = FC_tmp(drivPos(j));
                  disp(['coh_of_roi_with_driver{i}(j) = ' num2str(coh_of_roi_with_driver{i}(j))])
                end
            end
            simResult.coh_of_roi_with_driver = coh_of_roi_with_driver;
        end

        if ~isempty(this.params.JansenRitConnectomePaper.drivPosVarMatrix)
          simResultAllTmp{repeatIdx} = simResult;
        end
        
      end

      if ~isempty(this.params.JansenRitConnectomePaper.drivPosVarMatrix)
        simResult = simResultAllTmp;
      end
        
      % save results
      if ~exist(this.workpath, 'dir')
        mkdir(this.workpath)
      end
      save([this.workpath,'/', filename_JR], 'simResult')

      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function [all_missing_job_ids] = getUnfinishedJobIds(this)
      jobDesc = load( fullfile(this.temppath,'jobDesc.mat') );
      paramComb = jobDesc.paramComb;
      numJobs = size(paramComb, 2);
      all_missing_job_ids = [];
      for j=1:numJobs
          fname = fullfile( this.workpath, [this.params.Gridjob.jobname num2str(j) '.mat']);
          if ~exist(fname, 'file')
              disp([num2str(j) ',']);
              all_missing_job_ids(end+1) = j;
          end
      end
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      disp('finishing job now...')
      
      jobDesc = load( fullfile(this.temppath,'jobDesc.mat') );
      
      paramComb = jobDesc.paramComb;
      variableParams = jobDesc.variableParams;
      numJobs = size(paramComb, 2);
      
      paramValues = cell(1, length(jobDesc.variableParams));
      for k=1:length(jobDesc.variableParams)
          newStruct = jobDesc.params;
          for f=1:length(jobDesc.variableParams{k})
              fname = jobDesc.variableParams{k}{f};
              newStruct = newStruct.(fname);
          end
          paramValues{k} = newStruct;
      end
      
      %%
      all_coh = cell(1,numJobs);
      all_FC = cell(1,numJobs);
      all_corr_SimFC = cell(1,numJobs);
      all_missing_job_ids = [];
      for j=1:numJobs
          fname = fullfile( this.workpath, [this.params.Gridjob.jobname num2str(j) '.mat']);
          if exist(fname, 'file')
              tmp = load( fname );
              
              all_FC{j} = tmp.simResult.FC;
              if jobDesc.params.JansenRitConnectomePaper.corrSimFC
                all_corr_SimFC{j} = tmp.simResult.corr_SimFC;
              end
              
              if jobDesc.params.JansenRitConnectomePaper.calcCohWithDriver
                all_coh{j} = tmp.simResult.coh_of_roi_with_driver{1};
              end
          else
              disp(['file missing: ' fname]);
              all_missing_job_ids(end+1) = j;
          end
      end
      
      disp('all_missing_job_ids:');
      disp(all_missing_job_ids);
      
      if ~exist(this.resultpath, 'dir')
        mkdir(this.resultpath)
      end
      save(fullfile( this.resultpath, 'all_coh.mat'), 'all_coh', 'paramComb', 'variableParams', 'paramValues', 'all_FC', 'all_corr_SimFC')
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

