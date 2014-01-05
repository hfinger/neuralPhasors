function qstat()
%QSTAT Summary of this function goes here
%   Detailed explanation goes here

%check if qstat command is missing
if ispc
  paths = dataPaths();
  systemCmd = ['plink.exe ' paths.plinkArg ' "qstat"'];
  system(systemCmd)
elseif system('command -v qstat >/dev/null 2>&1')
  system('ssh isonoe "qstat"');
else
  system('qstat');
end

end

