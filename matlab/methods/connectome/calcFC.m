
analyticSig = zeros(size(sourceRate));
for k=1:size(sourceRate,1)
  analyticSig(k,:) = hilbert(sourceRate(k,:));
end
phaseSig = analyticSig./abs(analyticSig);
plv = zeros(size(phaseSig,1),size(phaseSig,1));
for k=1:size(phaseSig,1)
  plv(k,:) = abs(mean( bsxfun(@times, conj(phaseSig(k,:)), phaseSig) , 2 ));
end
autospectrum = mean( conj(analyticSig) .* analyticSig , 2 );
crossspectrum = zeros(size(analyticSig,1),size(analyticSig,1));
meanImagCrossSpec = zeros(size(analyticSig,1),size(analyticSig,1));
meanAbsImagCrossSpec = zeros(size(analyticSig,1),size(analyticSig,1));
pli = zeros(size(analyticSig,1),size(analyticSig,1));
for k=1:size(analyticSig,1)
  X = bsxfun(@times, conj(analyticSig(k,:)), analyticSig);
  crossspectrum(k,:) = mean( X , 2 );
  meanImagCrossSpec(k,:) = mean(imag(X),2);
  meanAbsImagCrossSpec(k,:) = mean(abs(imag(X)),2);
  pli(k,:) = mean(sign(imag(X)),2);
end
wpli = meanImagCrossSpec./meanAbsImagCrossSpec;
coherency = crossspectrum ./ sqrt( autospectrum * autospectrum' );
coherence = abs(coherency);
icoh = imag(coherency);