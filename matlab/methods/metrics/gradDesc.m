function   [gSC, deltaSC, dev, glob] = gradDesc(SC, FC, mode, k1, lambda,extra)

%% define some masks
h1 = logical(diag(ones(length(SC)/2,1),+length(SC)/2));                     % get masks for LH2RH and
h2 = logical(diag(ones(length(SC)/2,1),-length(SC)/2));                     % RH2LH homotopic connections
h = h1 + h2;

ut = triu(true(size(FC)),+1);                                               % select upper triangular matrix
sut = (size(FC,1)*(size(FC,1)-1))./2;                                       % number of elements in upper triangular matrix
mut = (1:numel(FC))';                                                       % mask that can be set to upper triangular matrix

inter = logical(zeros(size(FC)));                                           % mask for interhemispheric connections
inter(1:length(FC)/2,length(FC)/2+1:end) = logical(1);
intra = logical(ut - inter);                                                % mask for intrahemispheric connections

%% load data, set parameters

paths = dataPaths();
load('/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIds.mat')
% load('C:\Users\PWJEbel\Desktop\USB goes HERE\this\databases\SC_Bastian\resortIds.mat') % get vector that sorts ROIs, eventually replace this
% resortIds = 1:66;

if (nargin<3)
  % load SC with sp = 0.5 , without added homotopic connections
  % load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics','ConnectomeMetrics','1SC.mat'));
  % SC = hSC;
  
  % load SC with sp = 0.0 (original data), without added homotopic connections
  load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
  SC = avg_ci;
  SC(isnan(SC)) = 0;
  SC = SC + SC';
  SC = normGraph(SC, avg_roi_size, 'ROIprd', false, 0);
  %SC = bsxfun(@rdivide,SC,sum(SC,2));
  
  dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat']);
  tmp = mean(dataEegTmp.coh_all([1:13 15 17:20],:,5:6,:,:,7),3);
  tmp = mean(tmp,1);
  FC = squeeze(tmp);
end

if (nargin<6)
    extra = 0;
end

dsteps = 1.0*1e-4;
nsteps = 75000;
expdec = linspace(-2.5,0,nsteps);                                           % define domain x of exp(-x), exp. decay learning rate 
% lambda = 0.5;                                                             % lambda = std([cc.c1]) ./ std([cc.c1]) , 300 for frob                                                                   
gamma = 1;
theta = 0.9;

path = strcat(paths.localTempDir,'/Results/','gradDesc/', strcat('reg', num2str(lambda),'/'),strcat('extra', num2str(extra)));
if ~exist(path,'dir')
  mkdir(path);
end

SC = bsxfun(@rdivide, SC, sum(SC,2));
gSC = SC;                                                                   % 0.5+0.5*rand(66,66); % Shuffle(Shuffle(SC)')
%gSC = bsxfun(@rdivide, gSC, sum(gSC,2));

dev = zeros(nsteps, 1);
sum(FC,2);

%% perform gradient descent

switch mode
  
  case 1
    dev = zeros(nsteps,1);
    for i = 1:nsteps
      [dev(i,1),grad, cc(i,1)] = froNorm(SC,FC,gSC,lambda,gamma,k1, true);
      gSC = gSC - max(1.0,exp(-expdec(i)))*1e-4*grad;
      gSC(gSC<0)=0;                                                         % correct for negative connectivity
      %gSC = bsxfun(@rdivide,gSC,sum(gSC,2));                               % row-renormalize again
      % dont do this -- leads to competition (and WTA?)
    end
    
  case 2
    mut = triu(true(size(FC)),+1);
    dev = zeros(nsteps,1);
    for i = 1:nsteps
      [dev(i,1),grad, cc(i,1)] = corNorm(SC,FC,gSC,lambda,gamma,k1,true);
      if ~mod(i,1000)
        disp(dev(i,1))
      end
      if(i ~= nsteps) 
        gSC = gSC - max(1.0,exp(-expdec(i)))*1e-4*grad;
        nSC = gSC;                                                           % remember negative-weight gSC for debugging purpose
        gSC(gSC<0)=0;                                                        % correct for negative connectivity
        %gSC = bsxfun(@rdivide,gSC,sum(gSC,2));                              % row-renormalize again
        %%rowsum = sum(gSC,2);
        %%gSC = mean(rowsum)*bsxfun(@rdivide,gSC,rowsum);
        % dont do this: leads to competition and WTA weight development
      end;
    end
    
    case 3 % same as case 2, but using RMS prop
        mut = triu(true(size(FC)),+1);
        dev = zeros(nsteps,1);
        r = 0;
        for i = 1:nsteps
            [dev(i,1),grad, cc(i,1)] = corNorm(SC,FC,gSC,lambda,gamma,k1,true);
            if ~mod(i,1000)
                disp(dev(i,1))
            end
            if(i ~= nsteps)                                                  % do not learn at the last step
                r = (1-theta)*grad.^2 + theta*r;
                v = dsteps*grad ./ sqrt(r); % max(1.0,exp(-expdec(i)))
                v(isnan(v)) = 0;
                gSC = gSC - v;
                %nSC = gSC;                                                  % remember negative-weight gSC for debugging purpose
                gSC(gSC<0)=0;                                                % correct for negative connectivity
                %gSC = bsxfun(@rdivide,gSC,sum(gSC,2));                      % row-renormalize again
                %%rowsum = sum(gSC,2);
                %%gSC = mean(rowsum)*bsxfun(@rdivide,gSC,rowsum);
                % dont do this: leads to competition and WTA weight development
            end;
        end
        
  case 4
    % same as case 2, but with video building
    mut = triu(true(size(FC)),+1);
    dev = zeros(nsteps,1);
    vidstep = 1;
    figure()
    for i = 1:nsteps
      [dev(i,1),grad, cc(i,1)] = corNorm(SC,FC,gSC,lambda,gamma,k1,true);
      gSC = gSC - dsteps*grad;
      gSC(gSC<0)=0;                                                         % correct for negative connectivity
      if ~mod(i,1000)
        disp(dev(i,1))
        % create video of matrix development
        imagesc(gSC);
        matdev(vidstep,1) = getframe();
        % create video of scatterplot development
        ggSC = gSC + gSC';  
        ggSC = bsxfun(@rdivide,ggSC,sum(ggSC,2)); 
        [~, corGC] = sar(ggSC,.65);   
        %hold on;
        scatter(FC(ut),corGC(ut), 'filled')
        %Yhat = FC(ut).*pGC(1) + pGC(2);
        %plot(FC(ut),Yhat,'black');       
        title('scatterplot FC vs SAR(gSC)');
        axis([0,1,0,1])
        ylabel('EEG')
        xlabel('SAR')
        scatterdev(vidstep,1) = getframe();
        vidstep = vidstep+1;
      end
      %gSC = bsxfun(@rdivide,gSC,sum(gSC,2));                               % row-renormalize again
      % dont do this -- leads to competition (and WTA?)
    end
    
  otherwise % perform unit test via finite difference approximation
    for i = 1:nsteps
      
      varc = eye(numel(SC),numel(SC));
      costR= zeros(numel(SC),1);
      costL= zeros(numel(SC),1);
      epsi = 1e-8;
      
      for j = 1:numel(gSC)
        [costR(j),~,~] = corNorm(SC,FC,gSC+epsi*reshape(varc(:,j),size(gSC)),0,1,k1,false);
        [costL(j),~,~] = corNorm(SC,FC,gSC-epsi*reshape(varc(:,j),size(gSC)),0,1,k1,false);
      end
      
      grad = (costR - costL) ./ (2*epsi);                             % two-side gradient approximation
      grad = reshape(grad, size(gSC,1), size(gSC,2));                 % vector to matrix reshaping
      grad(logical(eye(size(grad)))) = 0;                             % remove diagonal entries
      
      [~,gradtrue, ~] = corNorm(SC,FC,gSC,0,1,k1, true);
      disp(['corr: ' num2str(corr(grad(:),gradtrue(:)))])
      disp(['cosine similarity: ' num2str(dot(grad(:),gradtrue(:)) ./ (norm(grad(:)) *norm(gradtrue(:))))])
      
      gSC = gSC - dsteps*grad;                                        % after all calculations, apply grad
      
      %gSC(gSC<0)=0;                                                   % correct for negative connectivity
      %gSC = bsxfun(@rdivide,gSC,sum(gSC,2));                          % row-renormalize again
    end
    
end

%% analysis and plots

% better call normGraph here

ggSC = gSC;
%gSC = gSC + gSC';                                                          % make symmetric
%gSC = bsxfun(@rdivide,gSC,sum(gSC,2));                                     % row-normalize
%nSC = bsxfun(@rdivide,nSC,sum(nSC,2));  

% roweq = sum(gSC(1,:)); 
roweq = sum(gSC,2);                                                         % get sum of rows, equal for each row
gSC = bsxfun(@rdivide,gSC,roweq);                                           % row-normalize to 1 for comparability with SC
k2 = k1 * roweq;                                                            % global scaling parameter taking roweq into account

figure();
imagesc(gSC(resortIds,resortIds));
colorbar(); colormap(b2r(0,max([gSC(:);SC(:)])));
title('learned SC');
print(strcat(path,'/learnSC'),'-dpng')

figure();
imagesc(SC(resortIds,resortIds));
colorbar(); colormap(b2r(0,max([gSC(:);SC(:)])));
title('original SC');
print(strcat(path,'/origSC'),'-dpng')

figure();
deltaSC = gSC - SC;                                                             % calculate changes to original SC
imagesc(deltaSC(resortIds,resortIds));
% imagesc(dSC>0);
colorbar(); colormap(b2r(1*min(deltaSC(:)),1*max(deltaSC(:))));
title('changes to SC: gSC - SC');
print(strcat(path,'/SCdifference'),'-dpng')

figure()
subplot(2,1,1)
hist(SC(h==1),15);
title('weight distribution of default vs. learned homotopic connections');
xlabel('weight')
ylabel('number')
xlim([0 max([SC(h==1);gSC(h==1)])])
ylim([0 max([hist(SC(h==1),15),hist(gSC(h==1),15)])])
subplot(2,1,2)
hist(gSC(h==1),15);
xlabel('weight')
ylabel('number')
xlim([0 max([SC(h==1);gSC(h==1)])])
ylim([0 max([hist(SC(h==1),15),hist(gSC(h==1),15)])])
print(strcat(path,'/hDist'),'-dpng')

figure();
plot(dev)
ah = gca;
% location of the plot to be zoomed in, [x1 y1 x2 y2]
s_pos =[(nsteps-5000) min(0.975,dev(end)-0.02) nsteps min(1,dev(end)+0.02)];
% location of the zoom-in plot 
t_pos = [nsteps/2 0.7 nsteps 0.85];    
% generate a zoom-in plot. 
if(false) zoomPlot(ah, s_pos, t_pos); end;  
title(strcat('development of objective function', ', max.:', num2str(dev(end))));
% if lambda = 0, then this corresponds to: prediction performance corr(SAR(gSC),FC)
ylabel('regularization objective + prediction objective');
xlabel('iterations');
print(strcat(path,'/perfDev'),'-dpng')

figure();
hold on;
plot(lambda*[cc.c1])
plot(gamma*[cc.c2])
title(strcat('development of objectives, reg=',num2str(lambda)));
legend('prediction objective','regularization objective','Location','Southwest')
xlabel(strcat('iterations: [0,',num2str(nsteps), ']'));
print(strcat(path,'/objDev'),'-dpng')

% compare quality of prediction on learned gSC to default SC, 

[covSC, corSC] = sar(SC, k1);                                               % acquire SAR(SC)
[covGC, corGC] = sar(gSC,k2);                                               % acquire SAR(gSC)

tSC = triu(covSC,1);
tGC = triu(covGC,1);
tFC = triu(FC   ,1);                                                       

% (1) global prediction error (correlation)                                 % (over upper triangular matrices)
glob.CORRSC = corr(tSC(ut), tFC(ut));                                       % in case of lambda ~= 0:                                         
glob.CORRGC = corr(tGC(ut), tFC(ut));                                       % this is unequal dev(end)

glob.INTERSC = corr(tSC(inter), tFC(inter));                                % prediction performances for SC and gSC                                                
glob.INTERGC = corr(tGC(inter), tFC(inter));                                % on inter- or intrahemispherical connections     
glob.INTRASC = corr(tSC(intra), tFC(intra));                                                                             
glob.INTRAGC = corr(tGC(intra), tFC(intra));  

glob.interdSC= sum(deltaSC(inter));                                         % net change (plus or minus) of learned weights
glob.intradSC= sum(deltaSC(intra));                                         % on inter- or intrahemispherical connections

glob.globSCgSC = corr(SC(ut),gSC(ut));
glob.interSCgSC = corr(SC(inter),gSC(inter));
glob.intraSCgSC = corr(SC(intra),gSC(intra));

glob.globSCFC = corr(SC(ut),FC(ut));
glob.interSCFC = corr(SC(inter),FC(inter));
glob.intraSCFC = corr(SC(intra),FC(intra));

glob.globgSCFC = corr(FC(ut),gSC(ut));
glob.intergSCFC = corr(FC(inter),gSC(inter));
glob.intragSCFC = corr(FC(intra),gSC(intra));

% (2) local prediction error
[lSC, pSC, dSC] = fit_2D_data(FC(ut),corSC(ut), 'no');
eSC = reshape(dSC, size(FC(ut)));

[lGC, pGC, dGC] = fit_2D_data(FC(ut), corGC(ut),'no');
eGC = reshape(dGC, size(FC(ut)));

up = triu(ones(length(SC)),1);
up(~~up)= eGC;
loc.INTERGC = sum(up(inter));
loc.INTRAGC = sum(up(intra));

[~ , bins] = hist( [dSC,dGC] ,50);
figure()
hist( dSC , bins );
hd = findobj(gca,'Type','patch');
set(hd,'FaceColor','r','EdgeColor','w','facealpha',.75)
hold on;
hist( dGC , bins );
h2 = findobj(gca,'Type','patch');
set(h2,'EdgeColor','w','facealpha',0.45);
title('distribution of directed local error sizes')
legend('SAR(SC)','SAR(gSC)')
axis([-0.4,0.4,0,850])
print(strcat(path,'/SCgSCHist'),'-dpng')


figure()
up = triu(ones(length(SC)),1);
up(~~up)= eGC;
imagesc(up); colorbar(); colormap(b2r(min(eGC),max(eGC)));
title('local errors of SAR(gSC)');
print(strcat(path,'/localErr'),'-dpng')

figure()
subplot(1,2,1);
hold on;
RSS = dSC'*dSC;
TSS = sum((FC(ut)-mean(FC(ut))).^2);
r2 = 1 - RSS/TSS;
scatter(FC(ut), corSC(ut),'filled')  
Yhat = FC(ut).*pSC(1) + pSC(2); 
plot(FC(ut),Yhat,'black');
title(strcat('Rsquares=',num2str(r2)));
axis([0,1,0,1])
xlabel('EEG')
ylabel('SAR')

subplot(1,2,2);
hold on;
RSS = dGC'*dGC;
r2 = 1 - RSS/TSS;
scatter(FC(ut),corGC(ut), 'filled')   
Yhat = FC(ut).*pGC(1) + pGC(2); 
plot(FC(ut),Yhat,'black');
suptitle('scatterplot SAR(gSC) vs FC and TLS-fit');
title(strcat('Rsquares=',num2str(r2)));
axis([0,1,0,1])
xlabel('EEG')
ylabel('SAR')
print(strcat(path,'/TLSfit'),'-dpng')

% make some videos
if mode == 4
  MatVid = VideoWriter('matdev.avi');
  MatVid.FrameRate = 20;
  open(MatVid);
  hh=figure('Visible','off');
  for i=1:length(matdev)
    imshow(matdev(i,1).cdata); 
    hold on;
    %text(100,100,sprintf('Frame Number: %d',i));
    hold off;
    currFrame = getframe(hh);
    writeVideo(MatVid,currFrame);
  end
  close(MatVid);
  
  ScatterVid = VideoWriter('scatterdev.avi');
  ScatterVid.FrameRate = 20;
  open(ScatterVid);
  hh=figure('Visible','off');
  for i=1:length(scatterdev)
    imshow(scatterdev(i,1).cdata);
    hold on;
    %text(100,100,sprintf('Frame Number: %d',i));
    hold off;
    currFrame = getframe(hh);
    writeVideo(ScatterVid,currFrame);
  end
  close(ScatterVid);
end

savefilename = fullfile(path,'/');
save([savefilename 'gSC.mat'],'gSC','dSC', 'glob', 'loc', 'roweq', 'k2');


% compare changes to metrics --> might require sparseness!
sp = sum(sum(gSC==0))/numel(gSC);                       
[covSC, corSC] = sar(SC, .65); 
[covGC, corGC] = sar(gSC,.65);                    % ++regularization -> ++sp

% load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));% required to get avg_roi_size
% ggSC = normGraph(gSC, avg_roi_size,'none',false,0.7);

% sum(sum(mask .* (dSC>0))) ./ sum(sum(mask))                               % compare pos. changes to h connections ...
% sum(sum(dSC>0)) ./ numel(dSC)                                             % to pos. changes to all other connections
end

function [cost, grad, cc] = froNorm(SC, FC, gSC, lambda, gamma, k, grON)

b = inv(eye(size(gSC))-k*gSC);
bb = b*b';

% calculate cost-function
cstr = norm(gSC(:)-SC(:))^2;
cfnc = norm(FC(:) -bb(:))^2;
cc.c1 = cstr;
cc.c2 = cfnc;
cost = lambda*cc.c1 + gamma*cc.c2;

% calculate gradient matrix
if grON == true
  reg = 2*(gSC-SC);
  
  phi = FC - bb;                                                            % phi == outer derivative
  tmp1 = b' * phi * bb';
  tmp2 = b' * phi' * bb;
  des = 2*(-k)*(tmp1 + tmp2);
  
  grad = lambda*reg + gamma*des;
  grad(logical(eye(size(grad)))) = 0;                                       % remove diagonal entries
else
  grad = zeros(size(gSC));
end

end

%%

function [cost, grad, cc] = corNorm(SC, FC, gSC, lambda, gamma, k, grON)

b = inv(eye(size(gSC))-k*gSC);
bb = b*b';

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
comp2 = bb(ut)-mean(bb(ut));
nG = m*(comp1)'*(comp2);
dG = sqrt( s*sum((comp1).^2) )*sqrt( s*sum((comp2).^2) );

% calculate cost-function

cstr = norm(gSC(:)-SC(:))^2;                                                % nR / dR; % if cor instead of frob
cfnc = nG / dG;
cc.c1 = cstr;
cc.c2 = cfnc;
cost = lambda*cc.c1 + gamma*cc.c2;

% calculate gradient matrix
if grON == true
  % variables for gradient nominator
  phi1 = FC-mean(FC(triu(true(size(FC)),+1)));   
  phi1 = triu(phi1, +1);
  tmp11 = b' * phi1 * bb';
  tmp12 = b' * phi1' * bb;
  tmp1 = (-k)*(tmp11 + tmp12);                                              % (-k) for gradient ascent
  % variables for gradient denominator
  phi2 = bb-mean(bb(triu(true(size(bb)),+1)));
  phi2 = triu(phi2, +1);
  tmp21 = b' * phi2 * bb';
  tmp22 = b' * phi2' * bb;
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
  
  grad = lambda*reg + gamma*desc;
  grad = reshape(grad, size(FC));
  grad(logical(eye(size(grad)))) = 0;                                       % remove diagonal entries
else
  grad = zeros(size(gSC));
end

end

%%
% function [cov, cor] = sar(gSC, k) 
% if isscalar(k)
%   b = inv(eye(size(gSC))-k*gSC);
% else
%   b = inv(eye(size(gSC))-repmat(k,1,66).*gSC);
% end
% bb = b*b';
% 
% cov = bb;
% autocov = sqrt(diag(cov));
% cor = cov ./ (autocov * autocov');
% end