data = dataPaths();
[~,my_foldername] = fileparts(pwd);

my_path = ['Holger/2018_JR/ConnectomeNoDriversOpt_v_and_k/' my_foldername];
path_results = fullfile(data.resultsdir, my_path);
path_workdir = fullfile(data.workdir, my_path);

if ~exist(path_results)
  mkdir(path_results)
end

jobDesc = load( fullfile(path_workdir,'temp_Connectome','jobDesc.mat') );

paramComb = jobDesc.paramComb;
variableParams = jobDesc.variableParams;
numJobs = size(paramComb, 2);

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
all_coh = zeros(1,numJobs);
all_FC = cell(1,numJobs);
all_corr_SimFC = cell(1,numJobs);
for j=1:numJobs
  fname = fullfile( path_workdir, ['Connectome' num2str(j) '.mat']);
  if exist(fname, 'file')
      tmp = load( fname );

      all_FC{j} = tmp.simResult.FC;
      if jobDesc.params.JansenRitConnectomePaper.corrSimFC
        all_corr_SimFC{j} = tmp.simResult.corr_SimFC;
      end

      if jobDesc.params.JansenRitConnectomePaper.calcCohWithDriver
        all_coh(j) = tmp.simResult.coh_of_roi_with_driver{1};
      end
  else
      disp(['file missing: ' fname]);
  end
end

save(fullfile( path_results, 'all_coh.mat'), 'all_coh', 'paramComb', 'variableParams', 'paramValues', 'all_FC', 'all_corr_SimFC')
