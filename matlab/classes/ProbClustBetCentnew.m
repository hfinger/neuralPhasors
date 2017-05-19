classdef ProbClustBetCentnew < Gridjob
    %ProbtrackX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = ProbClustBetCentnew(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            this.params.ProbClustBetCentnew.split = num2cell([2:100:1000]);
            this.params.ProbClustBetCentnew.cosText = 'conn'; %cos
            this.params.ProbClustBetCentnew.splitType = 'NonRec'; %'NonRec'
            this.params.ProbClustBetCentnew.threshRange =  [1];
            this.params.ProbClustBetCentnew.WeighFacRange = ([0.5]);
            
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('ProbClustBetCentnew:this excecutes before the parallel job is started');
            
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
            
            [~, username] = system('whoami');
            userfolder = fullfile('/work', username(1:end-1));
            
            datapaths = dataPaths();
            datapaths.workdir
            
            disp ([num2str(this.params.ProbClustBetCentnew.split) ...
                   this.params.ProbClustBetCentnew.cosText ...
                   this.params.ProbClustBetCentnew.splitType ...
            num2str(this.params.ProbClustBetCentnew.threshRange)...
            num2str(this.params.ProbClustBetCentnew.WeighFacRange)]);
            
            
            
            waitForServer();
            subjNum = 1;
            WeighingFactor = this.params.ProbClustBetCentnew.WeighFacRange;
            threshRange = this.params.ProbClustBetCentnew.threshRange;
            cosText = this.params.ProbClustBetCentnew.cosText;
            splitType = this.params.ProbClustBetCentnew.splitType;
           PostProcessPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/';
            decayParam = -1;
            normBy = 'sum';
            clustRange = this.params.ProbClustBetCentnew.split: this.params.ProbClustBetCentnew.split + 99;
            
            
         
            clusterPath = [PostProcessPath splitType '/decay' num2str(decayParam)...
                'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh100'...
                 '/' cosText  '/' 'Subj'...
                num2str(subjNum) '/'];
            
            %      clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/Rec/ClustOrgW0/detailsClust'];
            %     OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/Rec/ClustOrgW0/StatandMet/'];
            
            betcentNew(clustRange, clusterPath, clusterPath);
            
            
            
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

