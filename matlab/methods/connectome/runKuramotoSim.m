function [ simResult ] = runKuramotoSim( sim, C , D )
%RUNSIM Summary of this function goes here
%   Detailed explanation goes here

startPhase = 2*pi*rand(size(C,1),1, 'double');
startState = bsxfun(@plus, startPhase, (0:sim.dt:1)*2*pi*sim.f);

Y = runKuramoto(C,D,sim.k,sim.f,sim.v,sim.t_max,sim.dt,sim.sampling,sim.sig_n,sim.d,sim.verbose,sim.approx,sim.invertSin,startState );

%%
simResult.Y = Y;
simResult.sim = sim;

end

