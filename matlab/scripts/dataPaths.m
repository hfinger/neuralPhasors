function data = dataPaths( )

if ispc
    basepath = 'E:';
else
    basepath = '/net/store/nbp/phasesim';
end

p = mfilename('fullpath');
p = strsplit(p,filesep);
srcFolder = p{end-3};

data.HOMEIMAGES = [basepath '/databases/LabelMeImages'];
data.HOMEANNOTATIONS = [basepath '/databases/LabelMeAnnotations'];
data.pdflatex = '/net/store/users/hofinger/texlive/2011/bin/x86_64-linux/pdflatex';
data.labelMeToolbox = [basepath '/databases/LabelMeToolbox'];
data.cifar100matlab = [basepath '/databases/cifar-100-matlab'];
data.databases = [basepath '/databases'];

data.paramdir = [basepath '/' srcFolder '/matlab/parameters'];
data.workdir = [basepath '/workdir'];
data.resultsdir = [basepath '/results'];

data.sge_workdir = ['/net/store/nbp/phasesim/workdir'];
data.sge_resultsdir = ['/net/store/nbp/phasesim/results'];
data.sge_pathToAddScriptPaths = ['/net/store/nbp/phasesim/' srcFolder '/matlab'];

% data.plinkArg = '-ssh -i D:\Daten\putty\puttyUOS\putty_nopw.ppk hofinger@gate.ikw.uni-osnabrueck.de';
% use putty gateway instead:
data.plinkArg = '-ssh -P 5555 -i D:\Daten\putty\puttyUOS\putty_nopw.ppk hofinger@localhost';

end
