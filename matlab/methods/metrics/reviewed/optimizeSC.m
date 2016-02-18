function [ learnSC, dev, logLearnSC, cc ] = optimizeSC( initSC, empSC, empFC, p )
%GRADDESCNEW Summary of this function goes here
%   Detailed explanation goes here

paths = dataPaths();

if ~isfield(p,'constrainSym')
  p.constrainSym = false;
end
if ~isfield(p,'useFrob')
  p.useFrob = false;
end
if ~isfield(p,'useRowRenorm')
  p.useRowRenorm = 'no'; % or 'yes' or 'allowScaling'       
end
if ~isfield(p,'useRowRenormButAllowScaling')
  p.useRowRenormButAllowScaling = false;
end
if ~isfield(p,'constrainPos')
  p.constrainPos = false;
end
if ~isfield(p,'k')
  p.k = 0.65;
end
if ~isfield(p,'lambda')
  p.lambda = 0;
end
if ~isfield(p,'gamma')
  p.gamma = 1;
end
if ~isfield(p,'maxSteps')
  p.maxSteps = 100000;
end
if ~isfield(p,'savepath')
  p.savepath = fullfile(paths.localTempDir,'gradDescResults');
end
if ~isfield(p,'unittest')
  p.unittest = false;
end
if ~isfield(p,'useRMSprop')
  p.useRMSprop = true;
end
if ~isfield(p,'useRMSprop')
  p.rmsPropDecay = 0.9;
end
if ~isfield(p,'useRMSprop')
  p.initLearnRate = 1e-4;
end
if ~isfield(p,'saveAtIters')
  p.saveAtIters = [];
end
if ~isfield(p,'breakAtCost')
  p.breakAtCost = Inf;
end
if ~isfield(p,'use_fmincon')
  p.use_fmincon = false;
end
if ~isfield(p,'loglevel')
  p.loglevel = 1;
end

learnRate = p.initLearnRate;
r = 0;
learnSC = initSC;
dev = zeros(p.maxSteps, 1);
cc = cell(p.maxSteps, 1);
logLearnSC = cell(length(p.saveAtIters),1);

if p.unittest
  
  varc = eye(numel(empSC),numel(empSC));
  costR= zeros(numel(empSC),1);
  costL= zeros(numel(empSC),1);
  epsi = 1e-8;
  
  if p.useFrob
    for j = 1:numel(learnSC)
      [costR(j)] = froNorm(empSC,empFC, learnSC+epsi*reshape(varc(:,j),size(learnSC)), p.lambda,p.gamma,p.k);
      [costL(j)] = froNorm(empSC,empFC, learnSC-epsi*reshape(varc(:,j),size(learnSC)), p.lambda,p.gamma,p.k);
    end
    [~,gradtrue, ~] = froNorm(empSC,empFC,learnSC,p.lambda,p.gamma,p.k);
  else
    for j = 1:numel(learnSC)
      [costR(j)] = corNorm(empSC,empFC, learnSC+epsi*reshape(varc(:,j),size(learnSC)), p.lambda, p.gamma, p.negPenalty, p.k, p.loglevel);
      [costL(j)] = corNorm(empSC,empFC, learnSC-epsi*reshape(varc(:,j),size(learnSC)), p.lambda, p.gamma, p.negPenalty, p.k, p.loglevel);
    end
    [~,gradtrue, ~] = corNorm(empSC,empFC,learnSC,p.lambda,p.gamma, p.negPenalty,p.k, p.loglevel);
  end
  grad = (costR - costL) ./ (2*epsi);                             % two-side gradient approximation
  grad = reshape(grad, size(learnSC,1), size(learnSC,2));                 % vector to matrix reshaping
  grad(logical(eye(size(grad)))) = 0;                             % remove diagonal entries
  disp(['corr: ' num2str(corr(grad(:),gradtrue(:)))])
  disp(['cosine similarity: ' num2str(dot(grad(:),gradtrue(:)) ./ (norm(grad(:)) *norm(gradtrue(:))))])
  
  
else
  
  if any(0==p.saveAtIters)
    logLearnSC{0==p.saveAtIters} = learnSC;
  end
  
  if p.use_fmincon
    optfun = @(learnSC) corNorm(empSC,empFC,learnSC,p.lambda,p.gamma, p.negPenalty,p.k, p.loglevel);
    
    opt = optimoptions(@fmincon,'Display','iter','GradObj','on','MaxFunEvals',20000);
    [learnSC, dev(1)] = fmincon(optfun, learnSC, [], [], [], [], -0.0000001*ones(size(empSC)), Inf(size(empSC)), [], opt);
    
%     opt = optimoptions(@fminunc,'Algorithm','trust-region','Display','iter','GradObj','on','MaxFunEvals',20000);
%     [learnSC, dev(1)] = fminunc(optfun, learnSC, opt);
    
  else
    for i = 1:p.maxSteps
    
      if p.useFrob
        [dev(i,1), grad, cc{i}] = froNorm(empSC,empFC,learnSC,p.lambda,p.gamma,p.k);
      else
        [dev(i,1), grad, cc{i}] = corNorm(empSC,empFC,learnSC,p.lambda,p.gamma, p.negPenalty,p.k, p.loglevel);
      end

      if p.constrainSym
        grad = grad + grad';
      end

      if p.useRMSprop
        r = (1-p.rmsPropDecay) * grad.^2 + p.rmsPropDecay * r;
        dSC = grad ./ sqrt(r);
        dSC(isnan(dSC)) = 0;
      else
        dSC = grad;
      end

      learnSC = learnSC - learnRate*dSC;

      if p.constrainPos
        learnSC(learnSC<0)=0;
      end

      switch p.useRowRenorm
        case 'no'

        case 'yes'
          learnSC = bsxfun(@rdivide,learnSC,sum(learnSC,2));
        case 'allowScaling'
          rowsum = sum(learnSC,2);
          learnSC = mean(rowsum)*bsxfun(@rdivide,learnSC,rowsum);
      end

      if mod(i,1000)==0
        disp(['iter=' num2str(i) ' cost=' num2str(dev(i,1),10) ' corr=' num2str(cc{i}.corr,10) ' regL2=' num2str(cc{i}.regL2,10) ' negPenalty=' num2str(cc{i}.neg,10)])
      end

      if any(i==p.saveAtIters)
        logLearnSC{i==p.saveAtIters} = learnSC;
      end

      if p.breakAtCost < cc{i}.corr
        disp(['break because corr(' num2str(i) ')=' num2str(cc{i}.corr,10) ' > ' num2str(p.breakAtCost,10)])
        break;
      end
    end
    
  end
  
  
end

end

function stop = outfun(x,optimValues,state)

disp(optimValues)

end
