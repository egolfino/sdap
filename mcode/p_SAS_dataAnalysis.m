function p_SAS_dataAnalysis()
%% Config
% subjIDs={'AS04','AS05','AS06','AS07','AS09','AS10','AS12','AS14','AS15','AS16','AS17','AS18','AS19','AS20','AS21','AS22','AS23'};
% subjIDs={'AS04', 'AS05', 'AS06', 'AS07', 'AS09', 'AS10', 'AS12', 'AS15', 'AS16', 'AS18', 'AS19', 'AS20', 'AS22'};
subjIDs={'AS06','AS07','AS09','AS10','AS15','AS16','AS17','AS18'};
rawDataDir='E:\DATA\SAS';

nReps.start=8;
nReps.ramp=2;
nReps.stay=8;
nReps.end=8;

nTrialsPerRep.start=8;
nTrialsPerRep.ramp=8;
nTrialsPerRep.stay=8;
nTrialsPerRep.end=8;

colors.up='r';
colors.down='b';

MVA_REPS=3;

%%
nTrialsTot=nReps.start*nTrialsPerRep.start+nReps.ramp*nTrialsPerRep.ramp+...
    nReps.stay*nTrialsPerRep.stay+nReps.end*nTrialsPerRep.end;
nRepsTot=nReps.start+nReps.ramp+nReps.stay+nReps.end;

mat_F1Mid.up=nan(0,nTrialsTot);         mat_F1Mid.down=nan(0,nTrialsTot);
mat_repMean_F1Mid.up=nan(0,nRepsTot);   mat_repMean_F1Mid.down=nan(0,nRepsTot);
dm_repMean_F1Mid.up=nan(0,nRepsTot);    dm_repMean_F1Mid.down=nan(0,nRepsTot);

for i1=1:numel(subjIDs)
    pdataNewFN=fullfile(rawDataDir,subjIDs{i1},'pdata_new.mat');
    exptFN=fullfile(rawDataDir,subjIDs{i1},'expt.mat');
    load(pdataNewFN);
    load(exptFN);
    
    if isequal(expt.subject.shiftDirection,'F1Up');
        sDir='up';
    else
        sDir='down';
    end
    
    t_F1Mid=[];
    t_repMean_F1Mid=[];
    flds=fields(nReps);
    for i2=1:numel(flds)
        fld=flds{i2};
        
        idxStage=fsic(pdata.stage,fld);
        for i3=1:nReps.(fld)
            idxRep=find(pdata.repNum(idxStage)==i3);
            idx=idxStage(idxRep);
            
            array_F1Mid=pdata.F1Mid(idx);
            if length(idx)<nTrialsPerRep.(fld)
                array_F1Mid=[array_F1Mid,nan(1,nTrialsPerRep.(fld)-length(idx))];
            end
            
            t_F1Mid=[t_F1Mid,array_F1Mid];
            t_repMean_F1Mid(end+1)=nanmean(array_F1Mid);
        end
    end
    
    t_repMean_F1Mid(1:nReps.start)=mva(t_repMean_F1Mid(1:nReps.start),MVA_REPS,'Hamming');
    t_repMean_F1Mid(nReps.start+1:nReps.start+nReps.ramp+nReps.stay)=mva(t_repMean_F1Mid(nReps.start+1:nReps.start+nReps.ramp+nReps.stay),MVA_REPS,'Hamming');
    t_repMean_F1Mid(nReps.start+nReps.ramp+nReps.stay+1:nReps.start+nReps.ramp+nReps.stay+nReps.end)=...
        mva(t_repMean_F1Mid(nReps.start+nReps.ramp+nReps.stay+1:nReps.start+nReps.ramp+nReps.stay+nReps.end),MVA_REPS,'Hamming');
    
    mat_F1Mid.(sDir)=[mat_F1Mid.(sDir); t_F1Mid];
    mat_repMean_F1Mid.(sDir)=[mat_repMean_F1Mid.(sDir); t_repMean_F1Mid];
    dm_repMean_F1Mid.(sDir)=[dm_repMean_F1Mid.(sDir); t_repMean_F1Mid-mean(t_repMean_F1Mid(1:nReps.start))];
end

mean_repMean_F1Mid.up=mean(mat_repMean_F1Mid.up);
mean_repMean_F1Mid.down=mean(mat_repMean_F1Mid.down);
sem_repMean_F1Mid.up=std(mat_repMean_F1Mid.up)/sqrt(size(mat_repMean_F1Mid.up,1));
sem_repMean_F1Mid.down=std(mat_repMean_F1Mid.down)/sqrt(size(mat_repMean_F1Mid.down,1));

% dm_repMean_F1Mid.up=mean_repMean_F1Mid.up-mean(mean_repMean_F1Mid.up(1:nReps.start));
% dm_repMean_F1Mid.down=mean_repMean_F1Mid.down-mean(mean_repMean_F1Mid.down(1:nReps.start));

% comp_repMean_F1Mid=[-dm_repMean_F1Mid.up;dm_repMean_F1Mid.down];
meanDM_repMean_F1Mid.up=mean(dm_repMean_F1Mid.up);
meanDM_repMean_F1Mid.down=mean(dm_repMean_F1Mid.down);
semDM_repMean_F1Mid.up=std(dm_repMean_F1Mid.up)/sqrt(size(dm_repMean_F1Mid.up,1));
semDM_repMean_F1Mid.down=std(dm_repMean_F1Mid.down)/sqrt(size(dm_repMean_F1Mid.down,1));

comp_repMean_F1Mid=[-dm_repMean_F1Mid.up;dm_repMean_F1Mid.down];
meanComp_repMean_F1Mid=nanmean(comp_repMean_F1Mid);
semComp_repMean_F1Mid=std(comp_repMean_F1Mid)/sqrt(size(comp_repMean_F1Mid,1));

%% Visualization
figure;
errorbar(1:nRepsTot,mean_repMean_F1Mid.up,sem_repMean_F1Mid.up,'color',colors.up);
hold on;
errorbar(1:nRepsTot,mean_repMean_F1Mid.down,sem_repMean_F1Mid.down,'color',colors.down);
set(gca,'XLim',[0,nRepsTot+1]);
ylabel('F1 prod (Hz, mean\pm1 SEM)');

%%
figure;
set(gca,'XLim',[0,nRepsTot+1]);
xs=get(gca,'XLim');
plot(xs,[0,0],'-','color',[0.5,0.5,0.5]);
hold on;
errorbar(1:nRepsTot,meanDM_repMean_F1Mid.up,semDM_repMean_F1Mid.up,'color',colors.up);
errorbar(1:nRepsTot,meanDM_repMean_F1Mid.down,semDM_repMean_F1Mid.down,'color',colors.down);
set(gca,'XLim',[0,nRepsTot+1]);

%%
ylabel('F1 prod change (Hz, mean\pm1 SEM)');
figure;
set(gca,'XLim',[0,nRepsTot+1]);
xs=get(gca,'XLim');
plot(xs,[0,0],'-','color',[0.5,0.5,0.5]);
hold on;
errorbar(1:nRepsTot,meanComp_repMean_F1Mid,semComp_repMean_F1Mid,'color','k');
ylabel('Compensation (Hz, mean\pm1 SEM)');
set(gca,'XLim',[0,nRepsTot+1]);

%% 
figure;
set(gca,'XLim',[0,nRepsTot+1]);
xs=get(gca,'XLim');
plot(xs,[0,0],'-','color',[0.5,0.5,0.5]);
hold on;
compStep=[0,diff(meanComp_repMean_F1Mid)];
compStep(1:nReps.start)=0;
plot(compStep,'ko-');
set(gca,'XLim',[0,nRepsTot+1]);
return