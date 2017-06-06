OutputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160613_IntersectingClusters/';
ComparisonText = 'FSCompare';
OverlapText = 'RemoveOverlap';
subjNum     = 1;
cosText     = 'conn';
RecText     = 'Rec';
decayParam  = '-1';
WeighRange = (0.5);
clustRange  = (2:100:1000);


for WeighFactor = WeighRange
    
    for clustNum = clustRange
        
        
        disp([num2str(WeighFactor) num2str(clustNum)]);

FinalOutputPath = [OutputPath ComparisonText OverlapText '/Subj' num2str(subjNum)... 
                '/' cosText '/' RecText '/decay' decayParam 'weigh' num2str(WeighFactor) '/'];
            
            tmp = load([FinalOutputPath 'IntersectRange' num2str(clustNum) 'to' num2str(clustNum + 99) '.mat']);
            
            tmp = tmp.finalret;
            
            if clustNum == clustRange(end)
           
                finalret.IntersectRatio(clustNum:1000, :)   = tmp.IntersectRatio(clustNum:1000, :);
                finalret.MaxIntRatio(clustNum:1000, :)      = tmp.MaxIntRatio(clustNum:1000, :);
                finalret.IndMaxIntRatio(clustNum:1000, :)   = tmp.IndMaxIntRatio(clustNum:1000, :);
            else
            
                finalret.IntersectRatio(clustNum:clustNum +99, :)   = tmp.IntersectRatio(clustNum:clustNum +99, :);
                finalret.MaxIntRatio(clustNum:clustNum +99, :)      = tmp.MaxIntRatio(clustNum:clustNum +99, :);
                finalret.IndMaxIntRatio(clustNum:clustNum +99, :)   = tmp.IndMaxIntRatio(clustNum:clustNum +99, :);
            end
    end
    
    save([FinalOutputPath 'finalret'], 'finalret' , '-v7.3');
end
            
            