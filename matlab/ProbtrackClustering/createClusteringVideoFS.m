%createClusteringVideo

prompt = {...
  'Enter subject range:',...
  'Enter CompSimPath:',...
  'Enter FSClusterPath:',...
  'Enter OutputPath:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {...
  '[1]',...
  '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/',...
  '/net/store/nbp/projects/phasesim/workdir/Arushi/20160707_voxelByIndFS/',...
  '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData_copy/FS_video'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
  subjRange     = str2num(answer{1});
  CompSimPath    = answer{2};
  FSClusterPath    = answer{3};
  OutputPath    = answer{4};
  
  for subjNum = subjRange
      disp(['subj' num2str(subjNum)]);
      createClusteringVideoFScalc( CompSimPath, [OutputPath '/subj' num2str(subjNum)], FSClusterPath );
  end
end
