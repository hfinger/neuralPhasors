classdef GenerateClustConnmat < Gridjob
    %ProbtrackX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = GenerateClustConnmat(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            this.params.GenerateClustConnmat.WeighFactor = num2cell(0:0.1:1);
            this.params.GenerateClustConnmat.clusterCount = (2:1000);
            this.params.GenerateClustConnmat.subjNum     = 1;
            this.params.GenerateClustConnmat.decayParam    = -1;
            this.params.GenerateClustConnmat.useCosineSim  = true;
            this.params.GenerateClustConnmat.recursiveSplit = true;
            this.params.GenerateClustConnmat.normBy        = 'sum';
            this.params.GenerateClustConnmat.WholeNormText = 'WholeMax';
            this.params.GenerateClustConnmat.CompSimPath   = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/';
            this.params.GenerateClustConnmat.GraclusPath   =  '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_GraclusCut/';
            this.params.GenerateClustConnmat.PostProcessPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_ClusteringPostprocessing/';
            this.params.GenerateClustConnmat.threshRange   = [1];
            
            
            
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('GenerateClustConnmat:this excecutes before the parallel job is started');
            
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
            
            if this.params.GenerateClustConnmat.useCosineSim
                cosText = 'cos';
            else
                cosText = 'conn';
            end
            if this.params.GenerateClustConnmat.recursiveSplit
                splitType = 'Rec';
            else
                splitType = 'NonRec';
            end
            GraclusPath = this.params.GenerateClustConnmat.GraclusPath;
            decayParam = this.params.GenerateClustConnmat.decayParam;
            WeighingFactor = this.params.GenerateClustConnmat.WeighFactor;
            subjNum =  this.params.GenerateClustConnmat.subjNum;
            normBy =     this.params.GenerateClustConnmat.normBy;
            threshRange = this.params.GenerateClustConnmat.threshRange;
            CompSimPath = this.params.GenerateClustConnmat.CompSimPath;
            PostProcessPath = this.params.GenerateClustConnmat.PostProcessPath;
            clustRange = this.params.GenerateClustConnmat.clusterCount;
            useCosineSim = this.params.GenerateClustConnmat.useCosineSim;
            recursiveSplit = this.params.GenerateClustConnmat.recursiveSplit;
            
            if useCosineSim
                cosText = 'cos';
            else
                cosText = 'conn';
            end
            if recursiveSplit
                splitType = 'Rec';
            else
                splitType = 'NonRec';
            end
            InputGraclusPath = [GraclusPath cosText '/' splitType '/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor) '/'...
                num2str(subjNum) '/normby' normBy 'thresh' ];
            
            for threshFactor = threshRange
                if exist([InputGraclusPath num2str(threshFactor)], 'dir')
                    threshFound = 1;
                    break;
                end
            end
            
            if ~threshFound
                error('No cluster matrix found for subj %i for weighingFactor %i',subjNum, WeighingFactor);
            end
              switch splitType
                case 'Rec'
                    clusterPath = [InputGraclusPath  num2str(threshFactor) '/graclusResultClust1000.mat' ];
                case 'NonRec'
                    clusterPath = [InputGraclusPath  num2str(threshFactor) '/graclusResultClustcollected.mat' ];
              end
            
                  CoordPath = [CompSimPath 'FinalCoord/' 'FinalCoord' num2str(subjNum)];
            
                    SPath = [CompSimPath 'SSymDiag0/SSymDiag0subj' num2str(subjNum)];
            
            
            OutputPath = [PostProcessPath splitType '/decay' num2str(decayParam)...
                'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                num2str(threshFactor) '/' cosText '/' 'Subj'...
                num2str(subjNum) '/'];
            
            GenerateClusterConnmat(clustRange, clusterPath, CoordPath, SPath, recursiveSplit, OutputPath)
                       
            
            
            
            
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

