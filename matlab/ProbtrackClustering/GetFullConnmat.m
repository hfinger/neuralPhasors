function  GetFullConnmat(inputPath,outputPath,subjTotal, splitPerSubject)
%GETFULLNFSCONNMAT Combine data to get connectivity matrix
%   Use raw data from probtrack for all subjects to get high res 
%   connectivity matrices and the corresponding waytotal matrices
if ~inputPath
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150804_ProbtrackXallsubjects';
end

if ~outputPath
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FullConnmat'];
if ~subjTotal
    subjTotal = 22;
end

if ~splitPerSubject
    splitPerSubject = 50;
end
for subjectNum = 1:subjTotal
    for split = 1:splitPerSubject
    connmatNum = ((subjectNum-1)*50)+split;
    connmattemp = load([inputPath '/connmat' num2str(connmatNum) '.mat']);
    waytotaltemp = load([inputPath '/waytotal' num2str(connmatNum) '.mat']);
    if split == 1
        connmat = connmattemp;
        waytotal = waytotaltemp;
    else
        connmat = vertcat(connmat, connmattemp);
        waytotal = vertcat(waytotal, waytotaltemp);
    end
    
    end
    
    save([outputPath '/connmatSubj' num2str(subjectNum)], connmat);
    save([outputPath '/waytotalSubj' num2str(subjectNum)], waytotal);
end


end

