clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsBinary';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsBinary';
params.ApplyWeights.weightsFile = cellfun(@(x) ['ReLuDAE100/' num2str(x) '/forwConn.mat'],num2cell([2 3 4 6 7 8 10 11 12]),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) (x>0);
params.ApplyWeights.plotColormap = 'hot';
paramsAll{1} = params;

% clear params;
% params.Gridjob.runLocal = true;
% params.Gridjob.requiremf = 12000;
% params.Gridjob.jobname = 'rfsNoAct';
% params.ApplyWeights.plotPdf = true;
% params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.fileid = 1;
% params.ApplyWeights.outActFolder = 'rfsNoAct';
% params.ApplyWeights.weightsFile = 'temp_ReLuDAE100/4/forwConnIter4000.mat';
% params.ApplyWeights.convType = 'same';
% params.ApplyWeights.actFcn = @(x) x; 
% params.ApplyWeights.actFcn2 = @(x) x;
% params.ApplyWeights.plotColormap = 'hot';
% paramsAll{2} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsWithActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsWithActFcn';
params.ApplyWeights.weightsFile = cellfun(@(x) ['ReLuDAE100/' num2str(x) '/forwConn.mat'],num2cell([2 3 4 6 7 8 10 11 12]),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
params.ApplyWeights.plotColormap = 'hot';
paramsAll{2} = params;

clear params;

gridjobs = Gridjob(paramsAll);
start(gridjobs);


dirs = dir(fullfile(gridjobs.workpath, gridjobs.params.ApplyWeights.outActFolder));
dirs = {dirs(3:end).name};
for i=1:length(dirs)
  movefile(fullfile(gridjobs.workpath,'rfsBinary',dirs{i},'act1-1.png'),fullfile(gridjobs.workpath, 'rfsBinary',['job' dirs{i} 'rfsBinary.png']));
  movefile(fullfile(gridjobs.workpath,'rfsWithActFcn',dirs{i},'act1-1.png'),fullfile(gridjobs.workpath, 'rfsWithActFcn',['job' dirs{i} 'rfsWithActFcn.png']));
  rmdir(fullfile(gridjobs.workpath, gridjobs.params.ApplyWeights.outActFolder,num2str(i)),'s')
end

