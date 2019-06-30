data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeOneDriver/' my_foldername]);
path_workdir = fullfile(data.workdir, ['Holger/2018_JR/ConnectomeOneDriver/' my_foldername]);

jobDesc = load( fullfile(data.workdir, 'Holger/2018_JR/ConnectomeOneDriver',my_foldername,'temp_Connectome','jobDesc.mat') );

paramComb = jobDesc.paramComb;
variableParams = jobDesc.variableParams;
params = jobDesc.params;
numJobs = size(paramComb, 2);
numRepeat = length(jobDesc.params.JansenRitConnectomePaper.drivPosVarMatrix);

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
all_coh = zeros(1,numJobs,numRepeat);
all_FC = cell(1,numJobs,numRepeat);
all_corr_SimFC = cell(1,numJobs,numRepeat);
all_meanKuramotoOrderParam = zeros(1,numJobs,numRepeat);
for j=1:numJobs
  fname = fullfile( path_workdir, ['Connectome' num2str(j) '.mat']);
  if exist(fname, 'file')
      tmp = load( fname );

      for k=1:numRepeat
        all_FC{1,j,k} = tmp.simResult{k}.FC;
        if jobDesc.params.JansenRitConnectomePaper.corrSimFC
          all_corr_SimFC{1,j,k} = tmp.simResult{k}.corr_SimFC;
        end

        if jobDesc.params.JansenRitConnectomePaper.calcCohWithDriver
          all_coh(1,j,k) = tmp.simResult{k}.coh_of_roi_with_driver{1};
        end
        
        if jobDesc.params.JansenRitConnectomePaper.calcKuramotoOrderParam
          all_meanKuramotoOrderParam(1,j,k) = tmp.simResult{k}.meanKuramotoOrderParam;
        end
        
      end
      
  else
      disp(['file missing: ' fname]);
  end
end

if ~exist(path_results)
  mkdir(path_results)
end
save(fullfile( path_results, 'all_coh.mat'), 'all_coh', 'params', 'paramComb', 'variableParams', 'paramValues', 'all_FC', 'all_corr_SimFC', 'all_meanKuramotoOrderParam')
