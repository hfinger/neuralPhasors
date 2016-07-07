function plotstatsandmetrics(filepath, subjRange)
for subjNum = subjRange
    for mode = 1:4
        switch mode
            case 1
                modeText = 'fsconn';
            case 2
                modeText = 'fscos';
            case 3
                modeText = 'fullconn';
            case 4
                modeText = 'fullcos';
        end
        
        finalfilepath = [filepath num2str(subjNum) '/' modeText 'statandmet.mat'];
        statandmettemp = load(finalfilepath);
        statandmettemp = statandmettemp.statandmet;
        statandmettemp = statandmettemp.(modeText);
        
        statandmet.clusterSecondMomentMean.(modeText) = cellfun(@nanmean, statandmettemp.clusterSecondMomentAll);
        statandmet.clusterSecondMomentStd.(modeText) = cellfun(@nanstd, statandmettemp.clusterSecondMomentAll);
        statandmet.clusterSizeMean.(modeText) = cellfun(@nanmean, statandmettemp.clusterSizeAll);
        statandmet.betcentMean.(modeText) = cellfun(@nanmean, statandmettemp.betcent);
        statandmet.clustcoefmean.(modeText) = cellfun(@nanmean, statandmettemp.clustcoef);
        %             statandmet.clustcoef.(modeText) = statandmet.clustcoef.(modeText)(~isnan(statandmet.clustcoef));
        statandmet.efficiency.(modeText) = cell2mat(statandmettemp.efficiency);

        statandmet.lambda.(modeText) = cell2mat(statandmettemp.lambda);
                        statandmet.(modeText) = statandmettemp;

    end
    save([filepath num2str(subjNum) '/statandmet.mat'], 'statandmet');
    cmap = lines(4);
    %% clusterSecondMomentMean
    yname = 'clusterSecondMomentMean';
    figure(5)
    hold on;
    for mode = 1:4
        switch mode
            case 1
                modeText = 'fsconn';
            case 2
                modeText = 'fscos';
            case 3
                modeText = 'fullconn';
            case 4
                modeText = 'fullcos';
        end
        plot(2:1001, statandmet.(yname).(modeText), 'color', cmap(mode,:))
    end
    ylabel('mean cluster second moment')
    xlabel('resolution [#ROI]')
    set(gca,'ylim',[0 15])
    legend(mode(:))
    saveas(gcf, ['figures/mean_' yname],'png');
    
    
    %% efficiency
    yname = 'efficiency';
    figure(2)
    hold on;
    for mode = 1:4
        switch mode
            case 1
                modeText = 'fsconn';
            case 2
                modeText = 'fscos';
            case 3
                modeText = 'fullconn';
            case 4
                modeText = 'fullcos';
        end
        plot(2:1001, statandmet.(yname).(modeText), 'color', cmap(mode,:))
    end
    ylabel(yname)
    xlabel('resolution [#ROI]')
    set(gca,'ylim',[0.006 0.02])
    legend(mode(:))
    saveas(gcf, ['figures/mean_' yname],'png');
    
%% lambda
yname = 'lambda';
figure(3)
hold on;
for mode = 1:4
        switch mode
            case 1
                modeText = 'fsconn';
            case 2
                modeText = 'fscos';
            case 3
                modeText = 'fullconn';
            case 4
                modeText = 'fullcos';
        end
  plot(2:1001, statandmet.(yname).(modeText), 'color', cmap(mode,:))
end
ylabel(yname)
xlabel('resolution [#ROI]')
% set(gca,'ylim',[0.006 0.02])
set(gca,'YScale','log');
legend(mode(:))
saveas(gcf, ['figures/mean_' yname],'png');

%% clustering_coef
yname = 'clustering_coef';
figure(4)
hold on;
for mode = 1:4
    switch mode
            case 1
                modeText = 'fsconn';
            case 2
                modeText = 'fscos';
            case 3
                modeText = 'fullconn';
            case 4
                modeText = 'fullcos';
        end
  plot(2:1001, statandmet.(yname).(modeText), 'color', cmap(mode,:))
end
ylabel('clustering coef')
xlabel('resolution [#ROI]')
set(gca,'ylim',[1e-4 1e-2])
set(gca,'YScale','log');
legend(mode(:))
saveas(gcf, ['figures/mean_' yname],'png');

    
end
end
