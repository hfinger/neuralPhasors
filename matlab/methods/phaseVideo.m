function phaseVideo( simfilePrefix,filename,fps,compression, iterations )
%PHASEVIDEO Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
  fps=3;
end
if nargin<5
  iterations = 1:1000;
end

for counter=1:length(iterations)
  i=iterations(counter);
  
  try
    tmp = load([simfilePrefix num2str(i) '.mat']);
  catch ex
    break;
  end
  currSeg = tmp.meanWeightedPhase;
  
  fh=sfigure(1);
  set(fh,'Color',[1 1 1]);
  
  imagesc(currSeg);
  title(['iteration ' num2str(i)]);
  colormap hsv;
  caxis([-pi pi]);
  colorbar;
  %axis off;
  axis image;
  drawnow;
  if nargin>1 && ~isempty(filename)
    movieFrames(counter)=getframe(fh);
  end
  
  lastSeg = currSeg;
  pause(1/fps);
end


if nargin>1 && ~isempty(filename)
  
  if nargin<4 || isempty(compression)
    movie2avi(movieFrames, [filename '.avi'], 'fps', fps);
  else
    %i.e. 'FFDS'
    movie2avi(movieFrames, [filename '.avi'], 'compression', compression, 'fps', fps);
  end
  linuxCompressVideo([filename '.avi']);
end
end

