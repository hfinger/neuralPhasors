load_prefix = '/net/store/nbp/projects/phasesim/results/rgast/JansenRitResults_';
save_path = '/net/store/nbp/projects/phasesim/results/Holger/JansenRitResults.mat';

%%
cohPerPhaseOffset = zeros(125,8);
Coherence = cell(125,8);
for k=1:125
  tmp = load([load_prefix num2str(k) '.mat']);
  
  cohPerPhaseOffset(k,:) = cell2mat(tmp.simResult.drivPosCohVar);
  Coherence(k,:) = tmp.simResult.Coherence;
end

%%
cohPerPhaseOffset = reshape(cohPerPhaseOffset, [5 5 5 8]);
Coherence = reshape(Coherence, [5 5 5 8]);

%%
stdCoh = std(cohPerPhaseOffset,[],4);

[I, J] = max(stdCoh(:));
[I1, I2, I3] = ind2sub([5 5 5], J);

save(save_path,'Coherence','cohPerPhaseOffset','stdCoh');