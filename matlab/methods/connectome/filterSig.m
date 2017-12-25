function [ sigFiltered ] = filterSig( sig, sr, useBP, useLP, bp, lp )
%FILTERSIG Function that applies a collection of bandpass and lowpass
%filters to a radial signal
% 
%   Input Parameters:
% 
%       sig - 2-dim array with radial signal, i.e. non-monotonic. 1. dim =
%             nodes, 2. dim = samples
%       sr - sampling rate of signal
%       useBP - if true, apply bandpass filters to signal
%       useLP - if true, apply lowpass filters to signal
%       bp - structure containing the bandpass filter parameters
%       lp - structure containing the lowpass filter parameters
% 
%   Output:
%
%       sigFiltered - Cell with arrays of same size as sig with filtered 
%                     signal

%% apply filters

sigFiltered = zeros(size(sig));

if useBP
    
    % create bandpass filter
    dBandPass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',bp.Fst1,bp.Fp1,bp.Fp2,bp.Fst2,bp.Ast1,bp.Ap,bp.Ast2,sr);
    HdBandPass = design(dBandPass,'butter');

    % apply bandpass filter
    for k=1:size(sig,1)
        sigFiltered(k,:) = filter(HdBandPass,sig(k,:));
    end
    
end

if useLP

    % create lowpass filter
    dLowPass = fdesign.lowpass('Fp,Fst,Ap,Ast',lp.Fp,lp.Fst,lp.Ap,lp.Ast,Fs);
    HdLowPass = design(dLowPass,'butter');

    % apply lowpass filter
    for k=1:size(sigFiltered,1)
        sigFiltered(k,:) = filter(HdLowPass,sigFiltered(k,:));
    end

end
    
end

