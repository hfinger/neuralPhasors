clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/SingleNode/' my_foldername]);

jobDesc = load( fullfile(data.workdir, 'Holger/2018_JR/SingleNode',my_foldername,'temp_Connectome','jobDesc.mat') );

fname = fullfile( path_results, ['Connectome1.mat']);
tmp = load( fname );

%%
figure(1)
plot(tmp.simResult.spectrumFreq, log(tmp.simResult.spectrumPower{1}))
xlabel('Frequency [Hz]')
ylabel('log power')
saveas(gcf, fullfile( path_results, 'spectrum_driver.png'))

%%
figure(2)
plot(tmp.simResult.spectrumFreq, log(tmp.simResult.spectrumPower{2}))
xlabel('Frequency [Hz]')
ylabel('log power')
saveas(gcf, fullfile( path_results, 'spectrum.png'))

%% calculate frequency with maximum power:
[maxVal,maxIdx] = max(tmp.simResult.spectrumPower{2});
freqWithMaxPower = tmp.simResult.spectrumFreq(maxIdx);
disp(freqWithMaxPower)