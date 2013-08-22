%% create testcase:
t=0:0.01:10;
f=[2 3 4];
amplitude=[1; 2; 3];

mixingMatrix = rand(3,3)+2*eye(3);

sig = bsxfun(@times,amplitude',sin((t'*f)*2*pi));
sig = sig + 0.3*randn(size(sig));

sigMixed = sig*mixingMatrix;


%% now orthogonalize:
epsilon = 1e-2;
sigma = sigMixed' * sigMixed / size(sigMixed,1);
[u, s] = svd(sigma);
ZCAWhite = u * diag(1 ./ sqrt(diag(s) + epsilon)) * u';
sigRec = sigMixed * ZCAWhite;

%% plots:

figure(1)

subplot(3,1,1)
plot(sig)
xlim([900 1000])
title(sprintf('original independent sources. correlations: %f %f %f', ...
  corr(sig(:,1),sig(:,2)), ...
  corr(sig(:,2),sig(:,3)), ...
  corr(sig(:,3),sig(:,1))))

subplot(3,1,2)
plot(sigMixed)
xlim([900 1000])
title(sprintf('mixed signals. correlations: %f %f %f', ...
  corr(sigMixed(:,1),sigMixed(:,2)), ...
  corr(sigMixed(:,2),sigMixed(:,3)), ...
  corr(sigMixed(:,3),sigMixed(:,1))))

subplot(3,1,3)
plot(sigRec)
xlim([900 1000])
title(sprintf('orthogonalized signals. correlations: %f %f %f', ...
  corr(sigRec(:,1),sigRec(:,2)), ...
  corr(sigRec(:,2),sigRec(:,3)), ...
  corr(sigRec(:,3),sigRec(:,1))))


