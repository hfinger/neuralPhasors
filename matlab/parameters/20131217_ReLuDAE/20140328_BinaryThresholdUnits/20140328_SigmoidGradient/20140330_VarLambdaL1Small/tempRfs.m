clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsWithActFcn';
params.Gridjob.deleteTmpFolder = true;
params.ApplyWeights.inActFolder = '../../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsWithActFcn';
params.ApplyWeights.weightsFile = cellfun(@(x) ['SigBiDAE100/' num2str(x) '/forwConn.mat'],num2cell(1:3),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn2 = @(x) (x>0.5);
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.plotColormap = 'hot';
paramsAll{1} = params;

clear params;

gridjobs = Gridjob(paramsAll);
gridjobs = start(gridjobs);


dirs = dir(fullfile(gridjobs.workpath, gridjobs.params.ApplyWeights.outActFolder));
dirs = {dirs(3:end).name};
for i=1:length(dirs)
  movefile(fullfile(gridjobs.workpath, gridjobs.params.ApplyWeights.outActFolder,num2str(i),'act1-1.png'),fullfile(gridjobs.workpath, gridjobs.params.ApplyWeights.outActFolder,['job' num2str(i) 'act.png']));
  rmdir(fullfile(gridjobs.workpath, gridjobs.params.ApplyWeights.outActFolder,num2str(i)),'s')
end
