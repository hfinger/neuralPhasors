sim.k=12; %global connection strength scaling
sim.f=40; % frequency of oscillators [Hz]
sim.v=10; % transmission velocity [m/s]
sim.t_max=10; %max simulation time [seconds]
sim.dt=0.0001; % simulation step size [seconds]
sim.sampling=10; % sampling every x steps
sim.sig_n=0; % amount of noise
sim.d=0; %initial interval to remove [seconds]
sim.verbose=true; % if we want to print time steps to console
sim.approx=false; %if we use approximation during simulation
sim.invertSin=false;

network.weightedNetwork=false; % should we use a weighted network (weights sampled from gamma distribution)
network.binaryNetwork=true; % or a sampled binary connection matrix

network.N=5;% number neurons per populations
network.numCluster = 3; % number of populations

network.k_intraCluster = 0.9; %connection strength within population
network.delay_intraClust = 2; %connection distance within population in ms
network.k_interCluster = 0.1; %connection strength between population
network.delay_interClust = 5; %connection distance between population in ms


% post processing parameters (i.e. bandpass filter):
env.t_rm=4; %initial time interval to remove (in seconds)

% first bandpass filter:
env.sigBandpass(1).Fst1 = 5.5; % [Hz] end stop band
env.sigBandpass(1).Fp1 = 6; % [Hz] start pass band
env.sigBandpass(1).Fp2 = 22; % [Hz] end pass band
env.sigBandpass(1).Fst2 = 22.5; % [Hz] start stop band
env.sigBandpass(1).Ast1 = 40;
env.sigBandpass(1).Ap = 1;
env.sigBandpass(1).Ast2 = 40;

env.sigBandpass(2).Fst1 = 29.5;
env.sigBandpass(2).Fp1 = 30;
env.sigBandpass(2).Fp2 = 48;
env.sigBandpass(2).Fst2 = 48.5;
env.sigBandpass(2).Ast1 = 40;
env.sigBandpass(2).Ap = 1;
env.sigBandpass(2).Ast2 = 40;

env.envLowpass.Fp = 0.5;
env.envLowpass.Fst = 1;
env.envLowpass.Ap = 1;
env.envLowpass.Ast = 40;

[ C, D ] = genNetwork( network );
[ simResult ] = runKuramotoSim( sim, C , D );
[ simEval ] = calcEnvFC( env, simResult );
save('sim4');