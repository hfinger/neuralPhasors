
% generate data
for subject=1:10
   x(:,subject) = [1:10]+randi(30,1);
   coef(subject) = (3+rand(1));
   y(:,subject) = coef(subject)*x(:,subject)+3*randn(10,1)- 5*mean(x(:,subject));
end
% create X, Y, subject for nlmefit
Y = y(:);
X = [x(:) ones(100,1)]; 
subject = sum(kron(diag(1:10),ones(10,1)),2);
% fit the data using 'model'
model = @(Betas,X) (X*Betas(:));
[Betas,PSI,stats] = nlmefit(X,Y,subject,[],model,[1 0]);

Res = reshape(stats.ires,10,10);
for s=1:10
  Yhat(:,s) = y(:,s) - Res(:,s);
end
SSeffect = norm(Yhat(:)-mean(Yhat(:))).^2;
SStotal = norm(Y-mean(Y)).^2;
R2 = SSeffect / SStotal;
SSerror = norm(Res(:)-mean(Res(:))).^2;
df = (rank(X)-1)+9; % add number of subjects -1? 
                  % alternatively it can be nb of subjects - rank(X)
dfe = stats.dfe; 
F = (SSeffect/df) / (SSerror/dfe);
p_val = 1-fcdf(F,df,dfe);