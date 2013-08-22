function exportTikZImage( filename, climPi, genPdf )
%CORRECTTIKZ Summary of this function goes here
%   Detailed explanation goes here

if nargin<3 || isempty(genPdf)
  genPdf = false;
end

matlab2tikz( [filename '.tikz'], 'width', '\figurewidth', 'relativePngPath', '', 'interpretTickLabelsAsTex', true );

text = fileread([filename '.tikz']);

text1='width=\\figurewidth,';
text2='hide x axis, hide y axis, width=\\figurewidth, ';
text = regexprep(text, text1, text2);

text1='anchor=left of south west,\nwidth=.*, height';
text2='anchor=left of south west,\nwidth=\\colorbarwidth, height';
text = regexprep(text, text1, text2);

text1='xtick=\\empty,';
text2='xshift=\\colorbarshift,xtick=\\empty,';
text = regexprep(text, text1, text2);

if nargin>1 && ~isempty(climPi) && climPi
  text1='xtick=\\empty, yticklabel pos=right';
  text2='xtick=\\empty, yticklabel pos=right, ytick={-3.14159,0,3.14159}, yticklabels={\$-\\pi\$,\$0\$,\$\\pi\$}';
  text = regexprep(text, text1, text2);
end

outfid = fopen([filename '.tikz'], 'wt');
fwrite(outfid, text );
fclose(outfid);

if genPdf
  pathstr = fileparts(mfilename('fullpath'));
%   copyfile(fullfile(pathstr,'image.tex'),fullfile(pwd,[savefile '.tex']));
  text = fileread(fullfile(pathstr,'image.tex'));
  
  [pathstr, name, ext] = fileparts([filename '.tikz']);
  
  text1='image.tikz';
  text2=[name ext];
  text = regexprep(text, text1, text2);

  outfid = fopen([filename '.tex'], 'wt');
  fwrite(outfid, text );
  fclose(outfid);
  
  data = dataPaths( );
  tmpPwd = pwd();
  if ~isempty(pathstr)
    cd(pathstr);
  end
  system([data.pdflatex ' ' fullfile([filename '.tex'])]);
  cd(tmpPwd);
end


end

