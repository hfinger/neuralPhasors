
conn.W = [];
conn.inputSubtract = zeros(1,1,3);
conn.inputSubtract(1,1,:) = [0.4045, 0.4165, 0.4027];

%%
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'lMLayer1';
params.ApplyWeights.inActFolder = '../20130318_AE400/labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'labelMeMeanSubtract';
params.ApplyWeights.weightsFile = conn;%'forwConn.mat';

%%
gridjobs = Gridjob(params);
start(gridjobs);


