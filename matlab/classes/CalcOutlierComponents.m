classdef CalcOutlierComponents < Gridjob
    %ProbtrackX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = CalcOutlierComponents(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            this.params.CalcOutlierComponents.WeighingFactor = num2cell([0:0.1:0.4 0.6:0.1:1]);
            this.params.CalcOutlierComponents.recursiveSplit = 'NonRec';
            this.params.CalcOutlierComponents.cosText        = 'conn';
            this.params.CalcOutlierComponents.clustRange     = 2:1000;
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('CalcOutlierComponents:this excecutes before the parallel job is started');
            
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
            RecursiveText = this.params.CalcOutlierComponents.recursiveSplit ;
            cosText = this.params.CalcOutlierComponents.cosText;
            clustRange =  this.params.CalcOutlierComponents.clustRange;
            for WeighingFactor = this.params.CalcOutlierComponents.WeighingFactor
                load('/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/FinalCoord/FinalCoord1.mat');
                clustPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/',...
                    cosText '/' RecursiveText '/decay-1weigh' num2str(WeighingFactor) '/1/normbysumthresh100/'];
                switch RecursiveText
                    case 'Rec'
                        load([clustPath 'graclusResultClust1000.mat']);
                        number_of_outliers = zeros(1000,1000);
                        number_of_connectedComponents = zeros(1000,1000);
                        number_of_outliers(:,1) = NaN;
                        number_of_connectedComponents(:,1) = NaN;
                        for clusterNum = clustRange
                            a = stepWiseClusters(:,clusterNum);
                            number_of_outliers(clusterNum+1:end,clusterNum) = NaN;
                            number_of_connectedComponents(clusterNum+1:end,clusterNum) = NaN;
                            
                            
                            %                 detClust = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW'...
                            %                     num2str(WeighingFactor) '/detailsClust' num2str(clusterNum) '.mat']);
                            
                            for singleCluster = 1:clusterNum
                                %                     voxel_coords = voxelCoordByCluster{singleCluster};
                                voxel_coords = FinalCoord((a==singleCluster), :);
                                [number_of_outliers(singleCluster, clusterNum), number_of_connectedComponents(singleCluster, clusterNum),clusterSecondMoment] = connectedComp(voxel_coords);
                            end
                        end
                        OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_ClusteringPostprocessing/' RecursiveText '/decay-1'...
                            'weigh' num2str(WeighingFactor) '/normbysumthresh1'...
                             '/' cosText '/' 'Subj1'...
                             '/'];
                        save([OutputPath 'components.mat'], 'number_of_outliers', 'number_of_connectedComponents', '-v7.3');
                        
                    case 'NonRec'
                        number_of_outliers = zeros(1000,1000);
                        number_of_connectedComponents = zeros(1000,1000);
                        number_of_outliers(:,1) = NaN;
                        number_of_connectedComponents(:,1) = NaN;
                        load([clustPath 'graclusResultClustcollected' '.mat']);
                        for clusterNum = clustRange
                            a = stepWiseClusters(:,clusterNum);
                            number_of_outliers(clusterNum+1:end,clusterNum) = NaN;
                            number_of_connectedComponents(clusterNum+1:end,clusterNum) = NaN;
                            
                            for singleCluster = 1:clusterNum
                                voxel_coords = FinalCoord((a==singleCluster), :);
                                [number_of_outliers(singleCluster, clusterNum), number_of_connectedComponents(singleCluster, clusterNum),clusterSecondMoment] = connectedComp(voxel_coords);
                                
                            end
                            
                            
                        end
                         OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/' RecursiveText '/decay-1'...
                            'weigh' num2str(WeighingFactor) '/normbysumthresh100'...
                             '/' cosText '/' 'Subj1'...
                             '/'];
                         if ~exist(OutputPath, 'dir')
                             mkdir(OutputPath);
                         end
                        save([OutputPath 'components.mat'], 'number_of_outliers', 'number_of_connectedComponents', '-v7.3');
                end
            end
            %             for WeighingFactor = this.params.CalcOutlierComponents.WeighingFactor
            %             OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW' num2str(WeighingFactor) '/Outliers'];
            %                 detClust = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW'...
            %                     num2str(WeighingFactor) '/detailsClust' num2str(clusterNum) '.mat']);
            %                 number_of_outliers = zeros(clusterNum,1);
            %                 number_of_connectedComponents = zeros(clusterNum,1);
            %                 for singleCluster = 1:clusterNum
            %                     voxel_coords = detClust.voxelCoordByCluster{singleCluster};
            %                     [number_of_outliers(singleCluster), number_of_connectedComponents(singleCluster),clusterSecondMoment] = connectedComp(voxel_coords);
            %                 end
            %                 if ~exist(OutputPath, 'dir');
            %                 mkdir(OutputPath);
            %             end
            %
            %             save([OutputPath '/outlier' num2str(clusterNum)], 'number_of_outliers','number_of_connectedComponents');
            %
            %             end
            
            
            
            
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

