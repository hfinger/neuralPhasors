function [ simEval ] = calcEnvFC( env, simResult, rateSignal, useLP)
%CALCENVFC Summary of this function goes here
%   Detailed explanation goes here

Fs=1/(simResult.sim.sampling*simResult.sim.dt);
if rateSignal
    rate = simResult.Y;
else
    rate = sin(simResult.Y);
end
for b=1:length(env.sigBandpass)
  bp = env.sigBandpass(b);
  dBandPass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',bp.Fst1,bp.Fp1,bp.Fp2,bp.Fst2,bp.Ast1,bp.Ap,bp.Ast2,Fs);
  HdBandPass = design(dBandPass,'butter');
  sourceBP = zeros(size(rate));
  for k=1:size(rate,1)
    sourceBP(k,:) = filter(HdBandPass,rate(k,:));
  end
  
  sigHilbert = zeros(size(sourceBP));
  for n=1:size(sigHilbert,1)
    sigHilbert(n,:) = hilbert(sourceBP(n,:));
  end
  envSig = abs(sigHilbert);
  phaseBP{b} = angle(sigHilbert);
  
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
  
  FC{b} = corr(envLP(:,env.t_rm*Fs:end)');
end

simEval.FC = FC;
simEval.envLP = envLP;
simEval.phaseBP = phaseBP;

end

