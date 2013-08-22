function [ results ] = phasesyncPerDist( inActFolder, inPhaseFolder, inIterations, radius )
%PHASESYNCPERDIST Summary of this function goes here
%   Example call:
% [ results ] = phasesyncPerDist( 'layer1ActRectified', 'layer1Phase', 0:5:30, 16 );

varParams = dir(inPhaseFolder);
numVarParms = length(varParams)-2;

results.dimensionLabels = {'iterations','radius','varparam'};
results.dimensionTicks = {inIterations,radius,1:numVarParms};
results.phasesync = zeros(length(inIterations),length(radius),numVarParms);
results.phasesync2 = zeros(length(inIterations),length(radius),numVarParms);
results.phasesync3 = zeros(length(inIterations),length(radius),numVarParms);

for i=1:numVarParms
  disp(num2str(i))
  
  %load activity:
  curPhaseFolder = fullfile(inPhaseFolder,num2str(i));
  [ paths, filenames ] = dirrec( curPhaseFolder , '.mat' );
  act = load(fullfile(inActFolder, paths{1}(length(curPhaseFolder)+2:end), 'act1.mat' ));
  
  %sum to 2d:
  locAct = sum( act.act, 3);
  
  for k=1:length(inIterations)
    %load phase:
    phase = load(fullfile(paths{1}, ['phaseIter' num2str(inIterations(k)) '.mat'] ) );

    %sum to 2d:
    locCoherence = sum( act.act .* exp(phase.phase * 1i), 3) ./ locAct;
    locCoherence2 = sum( act.act .* exp(2 * phase.phase * 1i), 3) ./ locAct;
    locCoherence3 = sum( act.act .* exp(3 * phase.phase * 1i), 3) ./ locAct;
    
    for j=1:length(radius)
      
      %filter the local coherence
      convFilter = getnhood(strel('disk', radius(j)-1 , 0));
      convFilter = convFilter / sum(convFilter(:));
      
      locPhasesync = abs(conv2(locCoherence,convFilter,'valid'));
      results.phasesync(k,j,i) = mean(locPhasesync(:));
      locPhasesync2 = abs(conv2(locCoherence2,convFilter,'valid'));
      results.phasesync2(k,j,i) = mean(locPhasesync2(:));
      locPhasesync3 = abs(conv2(locCoherence3,convFilter,'valid'));
      results.phasesync3(k,j,i) = mean(locPhasesync3(:));

    end
  end
end


end

