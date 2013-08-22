function stats = analyzeConnShortestPath( W, maxdx )
%ANALYZECONNECTIVITY Summary of this function goes here

fids=sort(unique(W.f0));

%% format of W:
% neuron(x,y,f1) connects to neuron(x+dx,y+dy,f0)

stats.shortestPathLengths = 100*ones(max(fids),max(fids));
for fid=1:max(fids)
  disp(num2str(fid));
  
  connMat = zeros(2*maxdx+1,2*maxdx+1,max(fids)); %stores the connected neurons from neuron fid
  
  connMat(maxdx+1,maxdx+1,fid) = 1;
    stats.shortestPathLengths(fid,fid) = 0;
    
  for pathlength=1:100
%     disp(['pathlength=' num2str(pathlength)]);
    numThreads = 1;
    cutsW = round(linspace(1,length(W.w)+1,numThreads+1));
    cutsWend = cutsW(2:end)-1;
    connMatNewPar = cell(numThreads,1);
    
    for th=1:numThreads
      connMatNewPar{th} = zeros(size(connMat)); % x,y,f
      
      for i=cutsW(th):cutsWend(th)
        
        %% Shift images according to dx:
        if W.dx(i)>=0
          xCutStartIdf0 = W.dx(i);
          xCutStartIdf1 = 0;
          xCutEndIdf0 = 0;
          xCutEndIdf1 = W.dx(i);
        else
          xCutStartIdf0 = 0;
          xCutStartIdf1 = -W.dx(i);
          xCutEndIdf0 = -W.dx(i);
          xCutEndIdf1 = 0;
        end
        
        %% Shift images according to dy:
        if W.dy(i)>=0
          yCutStartIdf0 = W.dy(i);
          yCutStartIdf1 = 0;
          yCutEndIdf0 = 0;
          yCutEndIdf1 = W.dy(i);
        else
          yCutStartIdf0 = 0;
          yCutStartIdf1 = -W.dy(i);
          yCutEndIdf0 = -W.dy(i);
          yCutEndIdf1 = 0;
        end
        
        %% calculate interaction from neuron(x,y,f1) to neuron(x+dx,y+dy,f0)
        connMatNewPar{th}(1+xCutStartIdf0:end-xCutEndIdf0,1+yCutStartIdf0:end-yCutEndIdf0,W.f0(i)) = ...
          connMatNewPar{th}(1+xCutStartIdf0:end-xCutEndIdf0,1+yCutStartIdf0:end-yCutEndIdf0,W.f0(i)) + ...
          connMat(1+xCutStartIdf1:end-xCutEndIdf1,1+yCutStartIdf1:end-yCutEndIdf1,W.f1(i));
        
      end
    end
    connMat = connMat + sum(cat(4,connMatNewPar{:}),4);

    %check which features are connected:
    connectedFids = find(squeeze(connMat(maxdx+1,maxdx+1,:)));
    stats.shortestPathLengths(fid,connectedFids) = min(stats.shortestPathLengths(fid,connectedFids),pathlength);
    
    if length(connectedFids)==length(fids)
      break;
    end
    
  end
  
end

end