function collectstatsandmetrics(filepath, subjRange)
    for subjNum = subjRange
        for mode = 1:4
            switch mode
                case 1
                    modeText = 'fscos';
                    clusterStart = 67;
                case 2
                    modeText = 'fsconn';
                    clusterStart = 67;
                case 3
                    modeText = 'fullcos';
                    clusterStart = 2;
                case 4
                    modeText = 'fullconn';
                    clusterStart = 2;
            end
        
            finalfilepath = [filepath num2str(subjNum) '/' modeText '/'];
            statfilepath = [finalfilepath 'stats.mat'];
            savepath = [filepath '/' num2str(subjNum) '/' modeText 'statandmet.mat'];
            
            if ~exist(statfilepath, 'file')
                continue
            end
            statandmet.(modeText) = load(statfilepath);
            
            for clusterNum = clusterStart:1001
                disp(['subj' num2str(subjNum) modeText num2str(clusterNum)]);
                othermetfilepath = [finalfilepath modeText num2str(clusterNum) '.mat'];
                betcentfilepath = [finalfilepath 'betCentr' num2str(clusterNum) '.mat'];
                othermet = load(othermetfilepath);
                betcent = load(betcentfilepath);
                statandmet.(modeText).betcent{clusterNum} = betcent.betCentr;
                statandmet.(modeText).clustcoef{clusterNum} = othermet.clustering_coef;
                statandmet.(modeText).efficiency{clusterNum} = othermet.efficiency;
                statandmet.(modeText).lambda{clusterNum} = othermet.lambda;
                statandmet.(modeText).meanShortestDist1 = othermet.meanShortestDist1;
                statandmet.(modeText).meanShortestDist2 = othermet.meanShortestDist2;
            end
            
            save(savepath, 'statandmet');
            clear statandmet;

        end
    end
end