data = dataPaths();
path_results = fullfile(data.resultsdir, 'rgast/JansenRit/JR_Paper/Driver_varied5_with_sigm_out');

all_coh = zeros(1,26);
for j=1:195
   tmp = load( fullfile( path_results, ['JR_Connectome_Driver' num2str(j) '.mat']));
   all_coh(j) = tmp.simResult.coh_of_roi_with_driver{1};
end

save(fullfile( path_results, 'all_coh.mat'), 'all_coh')