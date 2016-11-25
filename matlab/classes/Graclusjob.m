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
            
            this.params.Graclusjob.WeighFactor = num2cell(0:0.1:0.9);
            this.params.Graclusjob.clusterCount = num2cell(2:50:200);
            
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
            
            %call program to wait
            waitForServer();
            disp('do clustering')
            for clusterCount = this.params.Graclusjob.clusterCount:this.params.Graclusjob.clusterCount + 19
                WeighingFactor = this.params.Graclusjob.WeighFactor;
                recursiveSplit = false;
                subjRange = [1];%, 6:13, 15, 17:22];%1
                decayParam = -1; % -0.5, -1
                normBy = 'sum';
                if recursiveSplit
                    splitType = '';
                else
                    splitType = 'NonRec';
                end
                InputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/decay'...
                    num2str(decayParam) 'weigh' num2str(WeighingFactor)]; %'/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSubwithoutnormbymean/decay-0.25weigh0.5'; %
                OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/decay'...
                    num2str(decayParam) 'weigh' num2str(WeighingFactor) 'conn' splitType ];
                thresholdingFactor = 10;
                disp( ['WEIGHINGFACTOR:' num2str(WeighingFactor)]);
                graclusNcut(recursiveSplit, clusterCount, subjRange, InputPath, OutputPath, normBy, thresholdingFactor)
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

