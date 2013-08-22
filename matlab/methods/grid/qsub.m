function qsub( scriptname, arg_ids, delay )
%QSUB Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
  delay=0;
end

if system('command -v qsub >/dev/null 2>&1')
  qsubMissing = true;
else
  qsubMissing = false;
end

if delay==0 && isequal(diff(arg_ids),ones(size(diff(arg_ids))))
  if qsubMissing
    system(['rsh shaggy "cd ' pwd '; qsub -t ' num2str(arg_ids(1)) ':' num2str(arg_ids(end)) ' ' scriptname '"']);
  else
    system(['qsub -t ' num2str(arg_ids(1)) ':' num2str(arg_ids(end)) ' ' scriptname]);
  end
else
  for i=1:length(arg_ids)
    if qsubMissing
      system(['rsh shaggy "cd ' pwd '; qsub -t ' num2str(arg_ids(i)) ' ' scriptname '"']);
    else
      system(['qsub -t ' num2str(arg_ids(i)) ' ' scriptname]);
    end
    pause(delay);
  end
end
end
