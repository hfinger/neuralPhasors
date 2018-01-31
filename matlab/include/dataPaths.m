function data = dataPaths( )

if ispc
    basepath = 'E:';
else
    basepath = '/net/store/nbp/projects/phasesim';
    if ~exist(basepath, 'dir')
        basepath = '/media/hofinger/OS/Users/hofinger/phasesim/';
    end
end

p = mfilename('fullpath');

if exist('strsplit')
  p = strsplit(p,filesep);
else
  p = strread(p,'%s','delimiter',[filesep]);
end

srcFolder = p{end-3};

data.HOMEIMAGES = [basepath '/databases/LabelMeImages'];
data.HOMEANNOTATIONS = [basepath '/databases/LabelMeAnnotations'];
data.pdflatex = '/net/store/nbp/phasesim/texlive/2011/bin/x86_64-linux/pdflatex';
data.labelMeToolbox = [basepath '/databases/LabelMeToolbox'];
data.cifar100matlab = [basepath '/databases/cifar-100-matlab'];
data.databases = [basepath '/databases'];
data.dti2eeg = [basepath '/dti2eeg'];

data.paramdir = [basepath '/' srcFolder '/matlab/parameters'];
data.workdir = [basepath '/workdir'];
data.resultsdir = [basepath '/results'];

data.sge_workdir = [basepath '/workdir'];
data.sge_resultsdir = [basepath '/results'];
data.sge_pathToAddScriptPaths = [basepath '/' srcFolder '/matlab'];

% data.plinkArg = '-ssh -i D:\Daten\putty\puttyUOS\putty_nopw.ppk hofinger@gate.ikw.uni-osnabrueck.de';
% use putty gateway instead:
data.plinkArg = '-ssh -P 5555 -i D:\Daten\putty\puttyUOS\putty_nopw.ppk hofinger@localhost';

[~,name] = system('whoami');
data.localTempDir = fullfile('/work',name(1:end-1));

end
