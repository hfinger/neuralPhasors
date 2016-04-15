
function [cost, grad, cc] = corNorm(SC, FC, gSC, lambda, gamma, negPenalty, k, loglevel)

Q = inv(eye(size(gSC))-k*gSC);
QQT = Q*Q';

ut = triu(true(size(FC)),+1);                                               % select upper triangular matrix
sut = (size(FC,1)*(size(FC,1)-1))./2;                                       % number of elements in upper triangular matrix

m = 1;                                                                      % factor for mean:    1/sut; 
s = 1;                                                                      % factor for variance:1/(sut-1);

% regularization term, cor(SC, gSC): NOT USED, use frob instead
% comp1 = SC(ut)-mean(SC(ut));
% comp2 = gSC(ut)-mean(gSC(ut));
% nR = m*(comp1)'*(comp2);
% dR = sqrt( s*sum((comp1).^2) )*sqrt( s*sum((comp2).^2) );

% gradient term, cor(FC, SAR_cov(gSC))
comp1 = FC(ut)-mean(FC(ut));
comp2 = QQT(ut)-mean(QQT(ut));
nG = m*(comp1)'*(comp2);
dG = sqrt( s*sum((comp1).^2) )*sqrt( s*sum((comp2).^2) );

% calculate cost-function
cc.regL2 = norm(gSC(:)-SC(:))^2;                                                % nR / dR; % if cor instead of frob
cc.corr = nG / dG;
cc.neg = sum(max(-gSC(:),0));
cost = lambda*cc.regL2 - gamma*cc.corr + negPenalty*cc.neg;

if loglevel>2
  disp(['cost=' num2str(cost) ' corr=' num2str(cc.corr) ' regL2=' num2str(cc.regL2) ' neg=' num2str(cc.neg) ])
end

% calculate gradient matrix
if nargout>1
  % variables for gradient nominator
  phi1 = FC-mean(FC(triu(true(size(FC)),+1)));   
  phi1 = triu(phi1, +1);
  tmp11 = Q' * phi1 * QQT';
  tmp12 = Q' * phi1' * QQT;
  tmp1 = (-k)*(tmp11 + tmp12);                                              % (-k) for gradient ascent
  % variables for gradient denominator
  phi2 = QQT-mean(QQT(triu(true(size(QQT)),+1)));
  phi2 = triu(phi2, +1);
  tmp21 = Q' * phi2 * QQT';
  tmp22 = Q' * phi2' * QQT;
  tmp2 = (-k)*(tmp21 + tmp22);                                              % (-k) for gradient ascent

  % Regularization derivatives: NOT USED, use frob instead
  %dnR = (SC(:)-mean(SC(:)))*(1-1);
  %ddR = sqrt( sum((SC(:)-mean(SC(:))).^2) )* ...
  %      ( ((gSC(:)-mean(gSC(:)))'*(1-1)) / sqrt(sum((gSC(:)-mean(gSC(:))).^2)) );
      
  % Gradient derivatives:
  dnG = m * tmp1  - mean(phi1(ut)) * tmp1;
  ddG = sqrt( s*sum((comp1).^2) )* ...
        ( ( s * tmp2  - mean(phi2(ut)) * tmp2)...
          / sqrt( s*sum((comp2).^2)) );    
        
  reg = 2*(gSC-SC);                                                         % (dnR*dR - nR*ddR) / dR.^2; % if cor instead of frob 
  desc = (dnG*dG - nG*ddG) / dG.^2;
  negPenaltyGrad = -(gSC<0);

  grad = lambda*reg + gamma*desc + negPenalty*negPenaltyGrad;
  grad = reshape(grad, size(FC));
  grad(logical(eye(size(grad)))) = 0;                                       % remove diagonal entries
else
  grad = zeros(size(gSC));
end

end