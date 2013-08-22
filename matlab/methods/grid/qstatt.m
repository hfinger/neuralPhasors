function qstatt(  )
%QSTATT Summary of this function goes here
%   Detailed explanation goes here

%check if qstat command is missing
if system('command -v qstat -u \* >/dev/null 2>&1')
  system('ssh shaggy "qstat"');
else
  system('qstat -u \*');
end

end

