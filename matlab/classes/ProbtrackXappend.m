classdef ProbtrackXappend < Gridjob
  %ProbtrackX Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ProbtrackXappend(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ProbtrackXappend.split = num2cell(1:1);
      this.params.ProbtrackXappend.numberPerSplit = 1000;
      this.params.ProbtrackXappend.numberRegions = 52442;
      this.params.ProbtrackXappend.divisions = 53;
         
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
      
      disp('this excecutes before the parallel job is started');
      
      %%%% END EDIT HERE:                                        %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      disp('some job parameters:')
      disp(this.workpath);
      disp(this.temppath);
      disp(this.resultpath);
      disp(this.currJobid);
      datapaths = dataPaths();
      finalpath = [datapaths.workdir '/Arushi/20150422_ProbtrackX/'];

      finalconnmat = zeros(this.params.ProbtrackXappend.numberRegions,this.params.ProbtrackXappend.numberRegions);
      finalwaytotal = zeros(this.params.ProbtrackXappend.numberRegions,1);
     
      for i = 1:this.params.ProbtrackXappend.divisions
         if i == this.params.ProbtrackXappend.divisions
         maxIdx = this.params.ProbtrackXappend.numberRegions;
         else
             
         maxIdx = i*this.params.ProbtrackXappend.numberPerSplit;
         end
         minIdx = (i-1)*this.params.ProbtrackXappend.numberPerSplit+1;
         
         connmat = load([finalpath 'connmat' num2str(i) '.mat']);
         
         if i == this.params.ProbtrackXappend.divisions
            finalconnmat(minIdx:maxIdx,:) = connmat.connmat(1:this.params.ProbtrackXappend.numberRegions-((i-1)*1000),:);
         else
            finalconnmat(minIdx:maxIdx,:) = connmat.connmat;
         end
         clear connmat;
         
       
         waytotal = load([finalpath 'waytotal' num2str(i) '.mat']);
          if i == this.params.ProbtrackXappend.divisions
            finalwaytotal(minIdx:maxIdx) = waytotal.waytotal(1:this.params.ProbtrackXappend.numberRegions-((i-1)*1000));
         else
             finalwaytotal(minIdx:maxIdx) = waytotal.waytotal;
         end
         clear waytotal;
      end
       save([finalpath 'aggregate/' 'rawconnmat.mat'], 'finalconnmat', '-v7.3');
       save([finalpath 'aggregate/' 'finalwaytotal.mat'],'finalwaytotal');
             
      
      disp(result);
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      % TODO: do this in another job, because it requires much more RAM...
      
%       connmat = zeros(this.params.ProbtrackX.numberRegions,this.params.ProbtrackX.numberRegions);
%       waytotal = zeros(this.params.ProbtrackX.numberRegions,1);
%       
%       for i = 1:53
%         minIdx = (i-1)*this.params.ProbtrackX.numberPerSplit+1;
%         maxIdx = i*this.params.ProbtrackX.numberPerSplit;
%         
%         temp = load([this.workpath '/connmat' i '.mat']);
%         connmat(minIdx:maxIdx,:) = connmat(minIdx:maxIdx,:) + temp.connmat;
%         save([this.workpath '/connmat'] , 'connmat','-v7.3')
%         
%         tempwaytotal = load([this.workpath '/waytotal' i '.mat']);
%         waytotal(minIdx:maxIdx) = waytotal(minIdx:maxIdx) + tempwaytotal.waytotal;
%         save([this.workpath '/waytotal'],'waytotal');
%       end
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

