sim.k=12;
sim.f=40;
sim.v=10;
sim.t_max=10;
sim.dt=0.0001;
sim.sampling=10;
sim.sig_n=0;
sim.d=0;
sim.verbose=true;
sim.approx=false;
sim.invertSin=false;

network.weightedNetwork=false;
network.binaryNetwork=true;
network.N=5;
network.k_intraCluster = 0.9;
network.delay_intraClust = 2; %in ms
network.numCluster = 3;
network.k_interCluster = 0.1;
network.delay_interClust = 5; %in ms

env.t_rm=4;

env.sigBandpass(1).Fst1 = 5.5;
env.sigBandpass(1).Fp1 = 6;
env.sigBandpass(1).Fp2 = 22;
env.sigBandpass(1).Fst2 = 22.5;
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
save('sim');