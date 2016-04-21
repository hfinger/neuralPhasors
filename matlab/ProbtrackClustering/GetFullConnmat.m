function  GetFullConnmat(inputPath,outputPath,subjTotal, splitPerSubject)
%GETFULLNFSCONNMAT Combine data to get connectivity matrix
%   Use raw data from probtrack for all subjects to get high res 
%   connectivity matrices and the corresponding waytotal matrices
if ~inputPath
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150804_ProbtrackXallsubjects';
end

if ~outputPath
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FullConnmat'];
end
if ~subjTotal
    subjTotal = 22;
end


if ~splitPerSubject
    splitPerSubject = 50;
end
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end
    
for subjectNum = 1:subjTotal
    for split = 1:splitPerSubject
    
    connmatNum = ((subjectNum-1)*50)+split;
    connmatpath = [inputPath '/connmat' num2str(connmatNum) '.mat'];
    waytotalpath = [inputPath '/waytotal' num2str(connmatNum) '.mat'];
    if ~exist(connmatpath, 'file')
        disp([connmatpath 'for subject' subjectNum 'does not exist']);
        break
    elseif ~exist(waytotalpath, 'file')
        disp([waytotalpath 'for subject' subjectNum 'does not exist']);
        break
    end
    connmattemp = load(connmatpath);
    connmattemp = connmattemp.connmat;
    waytotaltemp = load(waytotalpath);
    waytotaltemp = waytotaltemp.waytotal;
    if split == 1
        connmat = connmattemp;
        waytotal = waytotaltemp;
    else
        connmat = vertcat(connmat, connmattemp);
        waytotal = vertcat(waytotal, waytotaltemp);
    end
    
    end
    
    save([outputPath '/connmatSubj' num2str(subjectNum)], 'connmat');
    save([outputPath '/waytotalSubj' num2str(subjectNum)], 'waytotal');
end


end


