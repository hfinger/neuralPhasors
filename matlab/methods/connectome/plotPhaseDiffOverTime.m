function [ phaseDiff ] = plotPhaseDiffOverTime( struct, targets, evalWindow, drivEnd, env )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(struct);
n = length(fnames);
dataTmp = struct.(fnames{1});
Fs=1/(dataTmp.sampling*dataTmp.dt);
evalWindow = evalWindow * Fs;
drivEnd = drivEnd * Fs;

% extract signal in given evalWindow from given target nodes
phaseData = zeros(n, length(targets), evalWindow(2)-evalWindow(1)+1);
for f=1:n
    data = struct.(fnames{f});
    phaseData(f,:,:) = data.Y(targets,evalWindow(1):evalWindow(2));
end

% create bandpass filter
bp = env.sigBandpass(1);
dBandPass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',bp.Fst1,bp.Fp1,bp.Fp2,bp.Fst2,bp.Ast1,bp.Ap,bp.Ast2,Fs);
HdBandPass = design(dBandPass,'butter');

% filter phaseData
sourceBP = zeros(size(phaseData));
for k=1:size(phaseData,1)
    sourceBP(k,:,:) = filter(HdBandPass,squeeze(phaseData(k,:,:)),2);
end

% apply hilbert transform and extract phase
sigHilbert = zeros(size(sourceBP));
for n=1:size(sigHilbert,1)
    sigHilbert(n,:,:) = (hilbert(squeeze(sourceBP(n,:,:))'))';
end
phaseBP = angle(sigHilbert);

% find constant lagg and shift second driver accordingly
OffsetColl = zeros(size(phaseBP, 1),1);
for l=1:size(phaseBP,1)
    [r, lags] = xcorr(squeeze(phaseBP(l,1,1:drivEnd - evalWindow(1))),squeeze(phaseBP(l,2,1:drivEnd - evalWindow(1))));
    idx = find(lags == 0);
    r(1:idx-0.01*Fs) = 0;
    r(idx+0.01*Fs:end) = 0;
    [~,idx] = max(r);
    Offset = lags(idx);
    if Offset >= 0
        for i=1:size(phaseBP, 3) - Offset
            phaseBP(l,2,i) = phaseBP(l,2,i+Offset);
        end
    else
        for i=1:size(phaseBP, 3) + Offset
            phaseBP(l,2,i-Offset) = phaseBP(l,2,i);
        end   
    end
    OffsetColl = Offset;
end

% plot phase difference
phaseDiff = mod(bsxfun(@minus,phaseBP(:,1,:),phaseBP(:,2,:)), 2*pi);
phaseDiff = squeeze(mean(phaseDiff,1));

figure()
title('Driver After Effects on Phase')
xvec = 0:size(phaseDiff, 1) - 1;
xvec = (xvec - (drivEnd - evalWindow(1)))/Fs;
maxOffset = max(OffsetColl);
minOffset = abs(min(OffsetColl));
plot(xvec(minOffset:end-maxOffset), phaseDiff(minOffset:end-maxOffset))
ylabel('phase')
xlabel('Seconds since driver was turned off')

end

