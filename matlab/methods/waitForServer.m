function [ S ] = waitForServer( varargin )
%SAFELOAD same arguements as the matlab build in load function
%   uses exponential back-off strategy for loading. is testing bandwidth by uploading random file and then downloading again.



%                    after sec | required bandwith [ns]
afterSecReqNs(1,:) = [0          250000000];
afterSecReqNs(2,:) = [120        270000000];
afterSecReqNs(3,:) = [240        300000000];
afterSecReqNs(4,:) = [480        600000000];
afterSecReqNs(5,:) = [600        1000000000];
afterSecReqNs(6,:) = [1200 Inf];

tic;
if ~ispc
  counter = 1;
  waittime = 1;
  while true
    [status, result] = system('/net/store/nbp/projects/phasesim/src_Arushi/matlab/methods/bandwidthUpDown.sh 1');
    if status
      disp('could not find bandwidthUpDown.sh... load anyway...')
      break;
    end
    ns = str2double(result);
    
    totalTimeSoFar = toc;
    id = find(afterSecReqNs(:,1) <= totalTimeSoFar,1,'last');
    
    if ns < afterSecReqNs(id,2)
      if counter>1
        disp(['now load (ns=' num2str(ns) ' < ' num2str(afterSecReqNs(id,2)) ')... waited in total ' num2str(totalTimeSoFar) ' seconds...']);
      end
      break;
    else
      disp(['netstore is overloaded (ns=' num2str(ns) ' > ' num2str(afterSecReqNs(id,2)) ')... wait for ' num2str(waittime) ' seconds...']);
      pause(waittime);
      counter = counter + 1;
      
      [status, result] = system('cat /dev/urandom | tr -cd 0-9 | head -c 3'); 
      result=str2num(result); 
      sysrnd = result/1000;
      
      waittime = waittime*(sysrnd+1.5);
    end
  end
else
end

end

