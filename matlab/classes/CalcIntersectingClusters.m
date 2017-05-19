classdef CalcIntersectingClusters < Gridjob
    %ProbtrackX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        %% Subclass Constructor: initialize standard parameters:
        function this = CalcIntersectingClusters(varargin)
            
            % Call superclass constructor:
            this = this@Gridjob(varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: define standard parameters for the job %%%%
            this.params.CalcIntersectingClusters.split = num2cell(2:100:1000);
            this.params.CalcIntersectingClusters.cosText = 'conn';
            this.params.CalcIntersectingClusters.splitType = 'NonRec';
            this.params.CalcIntersectingClusters.threshRange = [100];
            this.params.CalcIntersectingClusters.subjRange = 1;
            this.params.CalcIntersectingClusters.WeighFacRange = num2cell(0:0.1:1);
            
            %%%% END EDIT HERE:                                          %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this = this.init(varargin{:});
            
        end
        
        %% Start: is executed before all individual parameter jobs are started
        function startJob(this)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
            
            disp('CalcIntersectingClusters:this excecutes before the parallel job is started');
            
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
             clustRange = this.params.CalcIntersectingClusters.split : this.params.CalcIntersectingClusters.split +99 ;
            threshRange = this.params.CalcIntersectingClusters.threshRange;
            cosText = this.params.CalcIntersectingClusters.cosText;
            splitType = this.params.CalcIntersectingClusters.splitType;
            subjRange = this.params.CalcIntersectingClusters.subjRange;
%             FSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160606_voxelByIndFS/';
            PostProcessPathNon = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/';
            PostProcessPathRec = '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_ClusteringPostprocessing/';
            OutputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_IntersectingClusters/';
            decayParam = '-1';
            normBy = 'sum';
            WeighFacRange = this.params.CalcIntersectingClusters.WeighFacRange;
            
            
%                 RecText = splitType;
                        
                        
            for subjNum = subjRange
                if subjNum >= 10
                    caNum = [num2str(subjNum)];
                else
                    caNum = ['0' num2str(subjNum)];
                end
                
              
                    OverlapText = 'RemoveOverlap';
              
%                 FinalFSPath = [FSPath caNum 'voxelByIndFS' OverlapText '.mat'];
%                 voxelByIndFS = load(FinalFSPath);
%                 voxelByIndFS = voxelByIndFS.voxelByInd;
                
                
                for WeighFactor = WeighFacRange
                    
                    
%                         FinalPostProcessPathNon = [PostProcessPathNon RecText '/decay' decayParam 'weigh'...
                        FinalPostProcessPathNon = [PostProcessPathNon 'NonRec/decay' decayParam 'weigh'...
                            num2str(WeighFactor) '/normby'  normBy 'thresh100' ...
                            '/' cosText '/Subj' num2str(subjNum) '/detailsClust.mat'];
                        
%                         FinalPostProcessPathRec = [PostProcessPathRec RecText '/decay' decayParam 'weigh'...
                        FinalPostProcessPathRec = [PostProcessPathRec 'Rec/decay' decayParam 'weigh'...
                            num2str(WeighFactor) '/normby'  normBy 'thresh1' ...
                             '/' cosText '/Subj' num2str(subjNum) '/detailsClust2to1000.mat'];
                        
                        
                    
                    
                    detailsClustNon = load(FinalPostProcessPathNon);
                    voxelIndAllNon = detailsClustNon.voxelIndByCluster;
                    clear detailsClustNon;
                    disp('x');
                    
                     detailsClustRec = load(FinalPostProcessPathRec);
                    voxelIndAllRec = detailsClustRec.voxelIndByCluster;
                    clear detailsClustRec;
                    disp('x');
                    
                    for clusterNum = clustRange
                        
                        if clusterNum > 1000
                            break;
                        end
                        disp(num2str(clusterNum));
                        
                        voxelIndByCluster1 = voxelIndAllNon(:,clusterNum);
                        voxelIndByCluster2 = voxelIndAllRec(:,clusterNum);
%                         voxelIndByCluster2 = voxelByIndFS;
                        
                        ret = calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2);
                        
                       
                        finalret.MaxIntRatio{clusterNum,1} = ret.MaxIntRatio1;
                        finalret.MaxIntRatio{clusterNum,2} = ret.MaxIntRatio2;
                       
                        clear ret;
                    end
%                     ComparisonText = 'FSCompare';
                    ComparisonText = 'RecNonRecCompare';
%                     FinalOutputPath = [OutputPath ComparisonText OverlapText '/Subj' num2str(subjNum)...
                      FinalOutputPath = [OutputPath ComparisonText '/Subj' num2str(subjNum)...
                        '/' cosText  '/decay' decayParam 'weigh' num2str(WeighFactor) '/'];
                    
                    if ~exist(FinalOutputPath, 'dir')
                        mkdir(FinalOutputPath);
                    end
                    save([FinalOutputPath 'IntersectRange' num2str(clustRange(1)) 'to' num2str(clustRange(end))],...
                        'finalret', '-v7.3');
                    disp(['saved subj' num2str(subjNum) 'weigh' num2str(WeighFactor)]);
                    
                    
                end
            end
            %             clusterNum = this.params.CalcIntersectingClusters.split;
            %             if clusterNum ~=1000
            %                 OutputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/RECNonREC/IntersectingClusters';
            %              voxelIndByClusterREC = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/Clustorg/detClust'...
            %                  num2str(clusterNum)]);
            %              voxelIndByClusterREC = voxelIndByClusterREC.voxelIndByCluster;
            %              voxelIndByClusterNonREC = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/NonRec/Clustorg/detailsClust'...
            %                  num2str(clusterNum)]);
            %              voxelIndByClusterNonREC = voxelIndByClusterNonREC.voxelIndByCluster;
            %
            %            calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2);
            %             end
            %
            
            
            
            
            
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

