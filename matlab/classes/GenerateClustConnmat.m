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
            this.params.GenerateClustConnmat.split = num2cell(2:1000);
            this.params.GenerateClustConnmat.WeighingFactor = num2cell(0:0.1:0.9);

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
            clusterNum = this.params.GenerateClustConnmat.split;
            SPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/SSymDiag0/SSymDiag0subj1.mat'];
            CoordPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/FinalCoord/FinalCoord1.mat';
            S = load(SPath);
            S = S.S;
            FinalCoord = load(CoordPath);
            FinalCoord = FinalCoord.FinalCoord;
            WeighingFactor = this.params.GenerateClustConnmat.WeighingFactor;

            clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/decay-1weigh' num2str(WeighingFactor) 'connnormbysumthresh10/graclusResultnormBysumthresh10clust1000subj1.mat'];
            OutputPath =  ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW' num2str(WeighingFactor)];
          



% zeroBothVoxelIdx = load(ZeroVoxelIdPath);
% zeroBothVoxelIdx = zeroBothVoxelIdx.zeroBothVoxelIdx;
% 
% S(zeroBothVoxelIdx,:) = [];
% S(:,zeroBothVoxelIdx) = [];


Cluster = load(clusterPath);
stepWiseClusters = Cluster.stepWiseClusters;

% clusterCenter = cell(clusterRange(end));
% clusterConnmat = cell(clusterRange(end),1);
% normClusterConnmat = cell(clusterRange(end),1);
% voxelIndByCluster = cell(clusterRange(end));
% voxelCoordByCluster = cell(clusterRange(end));

    clustConnmat = zeros(clusterNum);
    normClustConnmat = zeros(clusterNum);
    
    
            voxelIndByCluster = arrayfun(@(x) find(stepWiseClusters(:,clusterNum) == x), 1:clusterNum, 'UniformOutput', false)';
            voxelCoordByCluster = cellfun(@(x) FinalCoord(x,:), voxelIndByCluster, 'UniformOutput', false);
            clusterCenter = cellfun(@(x) mean(x,1), voxelCoordByCluster, 'UniformOutput', false);
            disp(num2str(clusterNum));
        for singleCluster = 1:clusterNum
            
            for otherCluster = 1:clusterNum
                a = S(voxelIndByCluster{singleCluster},voxelIndByCluster{otherCluster});
                clustConnmat(singleCluster, otherCluster) = sum(a(:));
                normClustConnmat(singleCluster, otherCluster) = clustConnmat(singleCluster, otherCluster)/(size(voxelIndByCluster{singleCluster},1) *size(voxelIndByCluster{otherCluster},1));
            end
        end
        
        if ~exist([OutputPath], 'dir')
            mkdir(OutputPath);
        end
        
       save([OutputPath '/detailsClust' num2str(clusterNum) ],  'clusterCenter','clustConnmat',...
    'normClustConnmat', 'voxelIndByCluster', 'voxelCoordByCluster', '-v7.3')

            
            
    
  


if ~exist(OutputPath, 'dir');
    mkdir(OutputPath);
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

