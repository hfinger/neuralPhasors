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
            this.params.ProbClustmetricsnew.split = num2cell(2:1000);
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
                        clusterNum = this.params.ProbClustmetricsnew.split;

            clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/NonRec/Clustorg/detailsClust' num2str(clusterNum) '.mat'];
            OutputPath =  ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/NonRec/StatandMet/'];
            Cluster = load(clusterPath);
            normClusterConnmat = Cluster.normClustConnmat;
            clusterCoords = Cluster.voxelCoordByCluster;
                
                %    betCentr = nan(clusterNum,1);
                %
                clustConnmat = normClusterConnmat;
                %     includeClusterIds = find(~isnan(sum(clustConnmat,2)));
                %     clustConnmat = clustConnmat(includeClusterIds,includeClusterIds);
                %
                L = 1./(clustConnmat);
                %
                %     betCentr(includeClusterIds) = betweenness_wei(L);
                %     betwCent{clusterNum} = betCentr;
                %
                %
                D = distance_wei(L);
                meanShortestDist1 = nan(clusterNum,1);
                meanShortestDist2 = nan(clusterNum,1);
                clustering_coef = nan(clusterNum,1);
                
                %     meanShortestDist1(includeClusterIds) = mean(D,1);
                meanShortestDist1 = mean(D,1);
                %     meanShortestDist2(includeClusterIds) = mean(D,2);
                
                meanShortestDist2 = mean(D,2);
                
                [lambda, efficiency] = charpath(D);
                
                %     clustering_coef(includeClusterIds) = clustering_coef_wu(clustConnmat);
                clustering_coef = clustering_coef_wu(clustConnmat);
                
                         
                clusterSize = cellfun(@(x) length(x), clusterCoords(:), 'UniformOutput', false);
                clusterSecondMoment = cellfun(@(x) mean(sqrt(sum(bsxfun(@minus,x,mean(x,1)).^2,2))), clusterCoords(:), 'UniformOutput', false);
                
                if ~exist(OutputPath, 'dir')
                    mkdir(OutputPath);
                end
        
            save ([OutputPath '/clusterstatandmetClust' num2str(clusterNum)], 'D', ...
                'meanShortestDist1', 'meanShortestDist2', 'clustering_coef',...
                'lambda', 'efficiency', 'clusterSize', 'clusterSecondMoment');
            
            
            
       
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

