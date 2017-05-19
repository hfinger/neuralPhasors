classdef ProbClustmetricsnew < Gridjob
    %ProbtrackX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = ProbClustmetricsnew(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            this.params.ProbClustmetricsnew.split = num2cell(0.5);
            this.params.ProbClustmetricsnew.cosText = 'cos'; %cos
            this.params.ProbClustmetricsnew.splitType = 'Rec'; %'NonRec'
            this.params.ProbClustmetricsnew.threshRange =  [1];
            
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('Probclustmetricsnew:this excecutes before the parallel job is started');
            
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
            
            
            
            
            waitForServer();
            subjNum = 1;
            WeighingFactor = this.params.ProbClustmetricsnew.split;
            threshRange = this.params.ProbClustmetricsnew.threshRange;
            cosText = this.params.ProbClustmetricsnew.cosText;
            splitType = this.params.ProbClustmetricsnew.splitType;
            PostProcessPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_ClusteringPostprocessing/';
            decayParam = -1;
            normBy = 'sum';
            clustRange = 2:1000;
            
            
            for threshFactor = threshRange
                if exist([PostProcessPath splitType '/decay' num2str(decayParam)...
                        'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                        num2str(threshFactor) '/' cosText  '/' 'Subj'...
                        num2str(subjNum) '/'], 'dir')
                    threshFound = 1;
                    break;
                end
            end
            
            if ~threshFound
                error('No cluster matrix found for subj %i for weighingFactor %i',subjNum, WeighingFactor);
            end
            clusterPath = [PostProcessPath splitType '/decay' num2str(decayParam)...
                'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                num2str(threshFactor) '/' cosText  '/' 'Subj'...
                num2str(subjNum) '/'];
            
            %      clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/Rec/ClustOrgW0/detailsClust'];
            %     OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/Rec/ClustOrgW0/StatandMet/'];
            
            clusterStatsNew(clustRange, clusterPath, clusterPath);
            
            
            
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

