data = dataPaths();
path_results = fullfile(data.resultsdir, 'rgast/JansenRit/JR_Paper/DriverVarFreq1_inp_out');

jobDesc = load( fullfile(data.workdir, 'rgast/JansenRit/JR_Paper/DriverVarFreq1_inp_out','temp_JR_Connectome_Driver','jobDesc.mat') );

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
for j=1:numJobs
    fname = fullfile( path_results, ['JR_Connectome_Driver' num2str(j) '.mat']);
    if exist(fname, 'file')
        tmp = load( fname );
        all_coh(j) = tmp.simResult.coh_of_roi_with_driver{1};
    else
        disp(['file missing: ' fname]);
    end
end

save(fullfile( path_results, 'all_coh.mat'), 'all_coh', 'paramComb', 'variableParams', 'paramValues')

