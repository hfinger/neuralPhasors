clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'test';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegSubjIds = [1:4 7:10];
params.ConnectomeEnvelopeReduce.outDirectory = 'test';
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


