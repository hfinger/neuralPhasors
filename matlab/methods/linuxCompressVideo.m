function linuxCompressVideo( filename )
%LINUXCOMPRESSVIDEO Summary of this function goes here
%   Detailed explanation goes here
if isunix
  oldLD_LIBRARY_PATH=getenv('LD_LIBRARY_PATH');
  try
    setenv('LD_LIBRARY_PATH','');
    [pathstr, ~, ~] = fileparts(filename);
    if isempty(pathstr)
      insertAtPos = 0;
    else
      insertAtPos = length(pathstr)+1;
    end
    tmpName = [filename(1:insertAtPos) 'tmp' filename(insertAtPos+1:end)];
    movefile(filename,tmpName);
    system(['ffmpeg -i ' tmpName ' -sameq ' filename]);
    delete(tmpName);
  catch ex
    setenv('LD_LIBRARY_PATH',oldLD_LIBRARY_PATH);
  end
end
end