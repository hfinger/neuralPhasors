function [optTheta, cost, exitflag, output, state] = minFuncGrad( fhandle, theta, options, lastState )
      %MINFUNCRMSPROP Summary of this function goes here
      %   fhandle is a function with two arguments: 1. theta, 2. patches
      
      learnRate = options.learnrate;
      
      if nargin<4 || isempty(lastState)
        if strcmp(options.Method,'rmsprop')
          msGrad = ones(size(theta)); %mean squared gradient
        elseif strcmp(options.Method,'stdprop')
          msGrad = ones(size(theta)); %mean squared gradient
          mGrad = ones(size(theta)); %mean gradient
        elseif strcmp(options.Method,'momentum')
          mGrad = zeros([length(options.momentum) size(theta)]);
        end
      else
        if strcmp(options.Method,'rmsprop')
          msGrad = lastState.msGrad;
        elseif strcmp(options.Method,'stdprop')
          msGrad = lastState.msGrad;
          mGrad = lastState.mGrad;
        elseif strcmp(options.Method,'momentum')
          mGrad = lastState.mGrad;
        end
      end
      
      optTheta = theta;
      
      %% outputFcn init
      i=0;
      iterationType.cost = [];
      iterationType.iteration = i;
      iterationType.optout = [];
      funEvals = i;
      feval(options.outputFcn,optTheta,iterationType,i,funEvals,[],[],[],[],[],[]);
      for i=1:options.maxIter
        iterationType = struct();
        
        %% evaluate gradient:
        [cost,grad,optout] = feval(fhandle,optTheta);
        
        if strcmp(options.Method,'rmsprop')
          %% update rms of gradient:
          EMAconst = options.EMAconst;
          EMAconst = max(EMAconst,1/i); %in the beginning weight rmsGrad not less than 1/i, because just i grad's were seen so far
          msGrad = (1-EMAconst) * msGrad + EMAconst * grad.^2;
          if EMAconst==options.EMAconst
            optTheta = optTheta - learnRate * grad ./ sqrt(msGrad); %dividing gradient by root mean squared gradient
          end
        
        elseif strcmp(options.Method,'stdprop')
          %% update rms of gradient:
          EMAconst = options.EMAconst;
          EMAconst = max(EMAconst,1/i); %in the beginning weight rmsGrad not less than 1/i, because just i grad's were seen so far
          msGrad = (1-EMAconst) * msGrad + EMAconst * grad.^2;
          mGrad = (1-EMAconst) * mGrad + EMAconst * grad;
          if EMAconst==options.EMAconst
            optTheta = optTheta - learnRate * grad ./ sqrt(msGrad-mGrad.^2); %dividing gradient by root mean squared gradient
          end
        elseif strcmp(options.Method,'momentum')
          mGrad = (1-options.momentum)*grad' + bsxfun(@times,mGrad,options.momentum);
          optTheta = optTheta - sum(bsxfun(@times,options.learnrate,mGrad),1)';
        end
        
        %% outputFcn
        iterationType.cost = cost;
        iterationType.iteration = i;
        iterationType.optout = optout;
        funEvals = i;
        feval(options.outputFcn,optTheta,iterationType,i,funEvals,[],[],[],[],[],[]);
        
        %% display
        if strcmp(options.display,'on')
          if mod(i,options.displayEvery)==0
            disp(['iter ' num2str(i) ' cost: ' num2str(cost)])
          end
        end
      end
      
      if nargout>4
        if strcmp(options.Method,'rmsprop')
          state.msGrad = msGrad;
        elseif strcmp(options.Method,'stdprop')
          state.msGrad = msGrad;
          state.mGrad = mGrad;
        elseif strcmp(options.Method,'momentum')
          state.mGrad = mGrad;
        end
      end
      
      exitflag = true;
      output.exitflag = exitflag;
      
    end