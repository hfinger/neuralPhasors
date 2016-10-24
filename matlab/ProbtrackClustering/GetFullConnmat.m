function  GetFullConnmat(inputPath,outputPath,subjRange, splitPerSubject)
%GETFULLNFSCONNMAT Combine data to get connectivity matrix
%   Use raw data from probtrack for all subjects to get high res
%   connectivity matrices and the corresponding waytotal matrices
if ~inputPath
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150804_ProbtrackXallsubjects';
end

if ~outputPath
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FullConnmat'];
end
if ~subjRange
    subjRange = [1:4, 6:13, 15, 17:22];
end


if ~splitPerSubject
    splitPerSubject = 50;
end
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

for subjectNum = subjRange
    for split = 1:splitPerSubject
        
        if subjectNum<10
            connmatNum = ['0' num2str(subjectNum)];
        else
            connmatNum = subjectNum;
        end
        
        connmatpath  = [inputPath '/connmat_subj_' num2str(connmatNum) '_split_' num2str(split) '.mat'];
        distmatpath  = [inputPath '/distmat_subj_' num2str(connmatNum) '_split_' num2str(split) '.mat'];
        waytotalpath = [inputPath '/waytotal_subj_' num2str(connmatNum) '_split_' num2str(split) '.mat'];
        if ~exist(connmatpath, 'file')
            disp([connmatpath 'for subject' num2str(subjectNum) 'does not exist']);
            break
        elseif ~exist(waytotalpath, 'file')
            disp([waytotalpath 'for subject' num2str(subjectNum) 'does not exist']);
            break
        end
        connmattemp = load(connmatpath);
        connmattemp = connmattemp.connmat;
        waytotaltemp = load(waytotalpath);
        waytotaltemp = waytotaltemp.waytotal;
        distmattemp  = load(distmatpath);
        distmattemp  = distmattemp.distmat;
        if split == 1
            connmat  = connmattemp;
            waytotal = waytotaltemp;
            distmat  = distmattemp;
        else
            connmat  = vertcat(connmat, connmattemp);
            waytotal = vertcat(waytotal, waytotaltemp);
            distmat  = vertcat(distmat, distmattemp);
        end
        
    end
    if exist('connmat', 'var')
        column = size(connmat,2);
        connmat = connmat(1:column,:);
        save([outputPath '/connmatSubj' num2str(subjectNum)], 'connmat');
    end
    if exist('waytotal', 'var')
        waytotal = waytotal(1:column);
        save([outputPath '/waytotalSubj' num2str(subjectNum)], 'waytotal');
    end
    if exist('distmat', 'var')
        waytotal = distmat(1:column,:);
        save([outputPath '/waytotalSubj' num2str(subjectNum)], 'waytotal');
    end
    clear connmat;
    clear waytotal;
    clear connmattemp;
    clear waytotaltemp;
end


end


