classdef Graclusjob < Gridjob
    %CompSimMatCalc Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = Graclusjob(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            
            this.params.Graclusjob.clusterCount = num2cell(2:40:1000);
            this.params.Graclusjob.WeighFactor = num2cell([0 0.2 0.4 0.5 0.6 0.8 1]);
            this.params.Graclusjob.subjNum       = 1;
            this.params.Graclusjob.decayParam    = -1;
            this.params.Graclusjob.threshFactor  = 1;
            this.params.Graclusjob.useCosineSim  = false;
            this.params.Graclusjob.normBy        = 'sum';
            this.params.Graclusjob.WholeNormText = 'WholeMax';
            this.params.Graclusjob.CompSimPath   = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/';
            this.params.Graclusjob.GraclusPath   =  '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_GraclusCut/';
            
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('Graclusjob:this excecutes before the parallel job is started');
            
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
            datapaths.workdir

            disp(this.params.Graclusjob.WeighFactor);
            disp(this.params.Graclusjob.clusterCount);
            
            WeighingFactor = this.params.Graclusjob.WeighFactor;
            subjNum = this.params.Graclusjob.subjNum;
            decayParam = this.params.Graclusjob.decayParam ;
            threshFactor = this.params.Graclusjob.threshFactor;
            useCosineSim = this.params.Graclusjob.useCosineSim;
            normBy = this.params.Graclusjob.normBy;
            WholeNormText =  this.params.Graclusjob.WholeNormText;
            CompSimPath = this.params.Graclusjob.CompSimPath ;
            GraclusPath = this.params.Graclusjob.GraclusPath ;
            recursiveSplit = false; % only keeping this for non-recursive
            %because recursive doesn't make sense
            %on grid job as each cut builds on the
            %previous one
            
            if recursiveSplit
                clustRange = clustRange(end);
                splitType = 'Rec';
            else
                splitType = 'NonRec';
            end
            if useCosineSim
                cosText = 'cos';
            else
                cosText = 'conn';
            end
           
            for clusterCount = this.params.Graclusjob.clusterCount:(this.params.Graclusjob.clusterCount+39)
                
                InputPath = [CompSimPath  cosText '/decay'...
                    num2str(decayParam) 'weigh' num2str(WeighingFactor)];
                OutputPath = [GraclusPath  cosText '/' splitType '/decay'...
                    num2str(decayParam) 'weigh' num2str(WeighingFactor) '/' ...
                   num2str(subjNum) '/normby' normBy 'thresh' num2str(threshFactor)];
                 disp( [cosText 'WEIGHINGFACTOR:' num2str(WeighingFactor) 'Thresh' num2str(threshFactor)]);
                
                compSimMat = load([InputPath '/compSimMat' WholeNormText '/' 'compSimMat' WholeNormText  normBy num2str(subjNum)]);
                
                compSimMat = compSimMat.compSimMat;
               
                            disp(['do clustering subj:' num2str(subjNum) 'clustercount:' num2str(clusterCount) ])
                            
                            [clusterIdPerVoxel, clusterIdPerVoxelCurrent, largestClusterId, cutValue]...
                                = applyClustering( compSimMat, clusterCount, recursiveSplit );
                           
                
                stepWiseClusters = clusterIdPerVoxel;
                allClusters = clusterIdPerVoxelCurrent;
                
                
                if ~exist(OutputPath, 'dir')
                    mkdir(OutputPath);
                end
                
                
                save([OutputPath '/graclusResultClust'  num2str(clusterCount)],...
                    'stepWiseClusters', 'allClusters', 'largestClusterId', 'cutValue');
                
            end
            
        
            
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
            %         minIdx = (i-1)*this.params.ProbtrackX.numberVoxThisSplit+1;
            %         maxIdx = i*this.params.ProbtrackX.numberVoxThisSplit;
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

