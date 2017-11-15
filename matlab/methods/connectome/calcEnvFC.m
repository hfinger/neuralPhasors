function [ simEval ] = calcEnvFC( env, useBP, simResult, radialSig, useLP)
%CALCENVFC Function that filters a signal with given filters and calculates
%phase and envelope of the filtered signal
%
%   Input Parameters:
% 
%       env - structure that contains fields with filter parameters
%       useBP - if true, use bandpass filters on env to filter the signal
%       simResult - strucutre that has a field sim with the simulation
%                   parameters and a field Y with the signal timeseries
%       radialSig - set to true, if signal is oscillating and not
%                   monotonically increasing
%       useLP - if true, use lowpass filters on env to filter the signal
%
%   Output:
%       
%       simEval - structure that includes a) a field FC with the functional
%                 connectivity matrix (based on correlation), b) the
%                 envelope of the signal envLP and c) the phase of the
%                 signal phaseBP

%% filter signal and extract envelope and phase

% get sampling rate and rate signal
Fs=1/(simResult.sim.sampling*simResult.sim.dt);
if radialSig
    rate = simResult.Y;
else
    rate = sin(simResult.Y);
end

for b=1:length(env.sigBandpass)
  
  % apply bandpass filters
  if useBP
      bp = env.sigBandpass(b);
      dBandPass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',bp.Fst1,bp.Fp1,bp.Fp2,bp.Fst2,bp.Ast1,bp.Ap,bp.Ast2,Fs);
      HdBandPass = design(dBandPass,'butter');
      sourceBP = zeros(size(rate));
      for k=1:size(rate,1)
        sourceBP(k,:) = filter(HdBandPass,rate(k,:));
      end
  else
      sourceBP = rate;
  end
  
  % apply hilbert transform
  sigHilbert = zeros(size(sourceBP));
  for n=1:size(sigHilbert,1)
    sigHilbert(n,:) = hilbert(sourceBP(n,:));
  end
  
  % get envelope and phase from hilbert transformed signal
  envSig = abs(sigHilbert);
  phaseBP{b} = angle(sigHilbert);
  
  % apply low pass filter to envelope
  if useLP
      lp = env.envLowpass;
      dLowPass = fdesign.lowpass('Fp,Fst,Ap,Ast',lp.Fp,lp.Fst,lp.Ap,lp.Ast,Fs);
      HdLowPass = design(dLowPass,'butter');
      envLP = zeros(size(envSig));
      for k=1:size(envSig,1)
        envLP(k,:) = filter(HdLowPass,envSig(k,:));
      end
  else
      envLP = envSig;
  end
  
  % calculate correlation between nodes based on envelope
  FC{b} = corr(envLP(:,env.t_rm*Fs:end)');
  
end

%% store results on structure

simEval.FC = FC;
simEval.envLP = envLP;
simEval.phaseBP = phaseBP;

end

