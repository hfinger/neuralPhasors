function qstat()
%QSTAT Summary of this function goes here
%   Detailed explanation goes here

%check if qstat command is missing
if system('command -v qstat >/dev/null 2>&1')
  system('ssh isonoe "qstat"');
else
  system('qstat');
end

end

