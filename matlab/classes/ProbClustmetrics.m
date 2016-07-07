classdef ProbClustmetrics < Gridjob
    %ProbtrackX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = ProbClustmetrics(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            this.params.ProbClustmetrics.split = num2cell(1:1276);
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('ProbtrackXall:this excecutes before the parallel job is started');
            
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
            subjNum = ceil(this.params.ProbClustmetrics.split/58);
            rangeflag = mod(this.params.ProbClustmetrics.split,58);
            if rangeflag == 0
                rangeflag = 58;
            end
            if rangeflag == 1
                clustRange = 2:201;
            elseif rangeflag == 2
                clustRange = 202:351;
            elseif (rangeflag >= 3) && (rangeflag <= 8)
                clustRange = 352+(25*(rangeflag-3)):376+(25*(rangeflag-3));
            elseif (rangeflag >= 9) && (rangeflag <=58)
                clustRange = 502+(10*(rangeflag-9)):511+(10*(rangeflag-9));
            end
            clusterPath = [datapaths.workdir '/Arushi/20150824Allsubjectrecursivencut/'];
            outputPath = [datapaths.workdir '/Arushi/20150920StatsAndMetrics/'];
            if clustRange
            for clusterTypenum = 1:3
                switch clusterTypenum
                    case 1
                        clusterType = 'fscos';
                    case 2
                        clusterType = 'fullconn';
                    case 3
                        clusterType = 'fsconn';
                    case 4
                        clusterType = 'fullcos';
                end
                
                for clusterNum = clustRange
                    %             clusterNum = mod(this.params.ProbClustbetcent.split, 1000) + 1;
                    %             if clusterNum == 1
                    %                 clusterNum = (clusterNum * 1000) + 1;
                    %             end
                    
                    if ((strcmp(clusterType, 'fsconn')) && (clusterNum > 66))
                        pass = 1;
                    elseif ((strcmp(clusterType, 'fscos')) && (clusterNum > 66))
                        pass = 1;
                    elseif (strcmp(clusterType,'fullcos'))
                        pass = 0;
                    elseif (strcmp(clusterType, 'fullconn'))
                        pass = 1;
                    else
                        pass = 0;
                    end
                    
                    if pass
                        if exist([clusterPath num2str(subjNum)], 'dir')
                            disp(['subj' num2str(subjNum) clusterType num2str(clusterNum)]);

                            subjClusterPath = [clusterPath num2str(subjNum) '/' clusterType '/' clusterType '/postprocessing/'];
                            
                            disp(['subj' num2str(subjNum) 'cluster' num2str(clusterNum)]);
                            filepath = [subjClusterPath 'clusterConnmat' num2str(clusterNum) '.mat'];
                            clusterConnmat = load(filepath);
                            clusterConnmat = clusterConnmat.clusterConnmat;
                            clusterConnmat = bsxfun(@rdivide, clusterConnmat, sum(clusterConnmat,2));
                            
                            includeClusterIds = find(~isnan(sum(clusterConnmat,2)));
                            
                            connSym = clusterConnmat + clusterConnmat';
                            connSym = connSym(includeClusterIds,includeClusterIds);
                            connSym(logical(eye(size(connSym)))) = 0;
                            L = 1./(connSym);
                            
                            clustering_coef = nan(size(clusterConnmat,1),1);
                            meanShortestDist1 = nan(size(clusterConnmat,1),1);
                            meanShortestDist2 = nan(size(clusterConnmat,1),1);

                            clustering_coef(includeClusterIds) = clustering_coef_wu(connSym);
                            D = distance_wei(L);
                            meanShortestDist1(includeClusterIds) = mean(D,1);
                            meanShortestDist2(includeClusterIds) = mean(D,2);
                            [lambda, efficiency] = charpath(D);
           
                            newOutputPath = [outputPath num2str(subjNum) '/' clusterType '/'];
                            if ~exist(newOutputPath, 'dir')
                                mkdir(newOutputPath);
                            end
           
                            save([newOutputPath clusterType num2str(clusterNum) '.mat'], 'meanShortestDist1', 'meanShortestDist2', 'lambda', 'efficiency', 'clustering_coef');
                                  
                        end
                    end
                    %%%% END EDIT HERE:                                %%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            end
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

