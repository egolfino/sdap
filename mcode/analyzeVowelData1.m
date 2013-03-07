function analyzeVowelData1(expDir)
%%
mvaWinWidth=24;

lim1=0.3;
lim2=0.7;
%%

load(fullfile(expDir,'expt.mat'));
load(fullfile(expDir,'micRMS_100dBA.mat'));

trainVowel.F1=struct;
trainVowel.F2=struct;
testVowel.F1=struct;
testVowel.F2=struct;

pdata.F1=[];
pdata.F2=[];
pdata.tVoiceBegin=[];       
pdata.tVoiceEnd=[];    
pdata.voiceDur=[];      
pdata.lv=[];
pdata.stage=cell(1,0);
pdata.file=cell(1,0);
pdata.rmsThresh=[];
pdata.x=[];

if expt.subject.designNum==1
    mdata=pdata;
end

for stg_={'start','ramp','stay','end'}
	stg=stg_{1};
	trainVowel.F1.(stg)=[];
	trainVowel.F2.(stg)=[];
	testVowel.F1.(stg)=[];
	testVowel.F2.(stg)=[];
	for n=1:expt.script.(stg).nReps
		d=dir(fullfile(expDir,stg,['rep',num2str(n)],['trial-*-1.mat']));
        if expt.subject.designNum==1
            dnm=dir(fullfile(expDir,stg,['rep',num2str(n)],['trial-*-2.mat']));
            d=[d;dnm];
        end
        
		for k=1:length(d)            
			load(fullfile(expDir,stg,['rep',num2str(n)],d(k).name));
            disp(['Loading data: ',fullfile(expDir,stg,['rep',num2str(n)],d(k).name)]);
            fn1=strrep(d(k).name,'.mat','');
            trialType=str2num(fn1(end));
            
			[i1,i2,jnk1,jnk2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));

            rms1=calcAWeightedRMS(data.signalIn(iv1*data.params.frameLen:iv2*data.params.frameLen-1),data.params.sr);
            
            trainVowel.F1.(stg)=[trainVowel.F1.(stg),mean(data.fmts(round(iv1+lim1*(iv2-iv1)):round(iv1+lim2*(iv2-iv1)),1))];
            trainVowel.F2.(stg)=[trainVowel.F2.(stg),mean(data.fmts(round(iv1+lim1*(iv2-iv1)):round(iv1+lim2*(iv2-iv1)),2))];            

            if isempty(iv1) || isempty(iv2)
                fprintf('\tWarning: no vowel was detected.\n');
                iv1=NaN;
                iv2=NaN;
            end
            
            
            if trialType==1
                pdata.tVoiceBegin=[pdata.tVoiceBegin,iv1*data.params.frameLen/data.params.sr];
                pdata.tVoiceEnd=[pdata.tVoiceEnd,iv2*data.params.frameLen/data.params.sr];
                pdata.voiceDur=[pdata.voiceDur,(iv2-iv1)*data.params.frameLen/data.params.sr];
                pdata.rmsThresh=[pdata.rmsThresh,data.params.rmsThresh];
                pdata.F1=[pdata.F1,trainVowel.F1.(stg)(end)];
                pdata.F2=[pdata.F2,trainVowel.F2.(stg)(end)];
                pdata.lv=[pdata.lv,100+20*log10(rms1/micRMS_100dBA)];
                pdata.stage{length(pdata.stage)+1}=stg;
                pdata.file{length(pdata.file)+1}=fullfile(stg,['rep',num2str(n)],d(k).name);
                if (isempty(pdata.x))
                    pdata.x=1;
                else
                    pdata.x=[pdata.x,pdata.x(end)+1];
                end
            elseif trialType==2
                mdata.tVoiceBegin=[mdata.tVoiceBegin,iv1*data.params.frameLen/data.params.sr];
                mdata.tVoiceEnd=[mdata.tVoiceEnd,iv2*data.params.frameLen/data.params.sr];
                mdata.voiceDur=[mdata.voiceDur,(iv2-iv1)*data.params.frameLen/data.params.sr];
                mdata.rmsThresh=[mdata.rmsThresh,data.params.rmsThresh];
                mdata.F1=[mdata.F1,trainVowel.F1.(stg)(end)];
                mdata.F2=[mdata.F2,trainVowel.F2.(stg)(end)];
                mdata.lv=[mdata.lv,100+20*log10(rms1/micRMS_100dBA)];
                mdata.stage{length(mdata.stage)+1}=stg;
                mdata.file{length(mdata.file)+1}=fullfile(stg,['rep',num2str(n)],d(k).name);
                if (isempty(mdata.x))
                    mdata.x=1;
                else
                    mdata.x=[mdata.x,mdata.x(end)+1];
                end
            end
            
            clear('data');
        end
        
		d=dir(fullfile(expDir,stg,['rep',num2str(n)],['trial-*-2.mat']));
		if ~isempty(d)
			for k=1:length(d)
				load(fullfile(expDir,stg,['rep',num2str(n)],d(k).name));
				[i1,i2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
				f1=data.fmts(i1:i2,1);
				f2=data.fmts(i1:i2,2);
				rms1=data.rms(i1:i2,1);
% 				rms1=mva(rms1,mvaWinWidth);
				[jnk,idxmax]=max(rms1);
				testVowel.F1.(stg)=[testVowel.F1.(stg),mean(f1(idxmax-5:idxmax+5))];
				testVowel.F2.(stg)=[testVowel.F2.(stg),mean(f2(idxmax-5:idxmax+5))];
				clear('data');
			end
		end
	end
end

%% Manual data inspection
items={'F1','F2','lv','tVoiceBegin','tVoiceEnd','voiceDur'};
hf1=figure('Position',[10,100,1200,600]);

if expt.subject.designNum==1
    data_list={'pdata','mdata'};
else
    data_list={'pdata'};
end

for k=1:length(data_list)
    eval(['idata=',data_list{k},';']);
    
    for n=1:length(items)        
        hs1=subplot(211);
        plot(idata.x,idata.(items{n}),'bo-'); hold on;
        ylabel([items{n},' (',data_list{k},')']);
        set(gca,'XLim',[0,idata.x(end)+1]);
        xs=get(gca,'XLim');

        x=ginput(1); x=x(1);
        while(x>=xs(1) && x<=xs(2))
            [jnk,idx]=min(abs(idata.x-x));
            plot(idata.x(idx),idata.(items{n})(idx),'bs-','LineWidth',1);

            load(fullfile(expDir,idata.file{idx}));
            subplot(212);
            cla;
            
        %     plot(data.fmts(:,1:2));
            [i1,i2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));	
            sigIn=data.signalIn;

            p=data.params;
            p.rmsThresh=idata.rmsThresh(idx);
            fs=testTSM2(data,p);
            data.fmts=fs;

            fs=data.params.sr;
            flen=data.params.frameLen;
            taxis1=0:(1/fs*flen):(1/fs*flen)*(size(data.fmts,1)-1);
            [s,f,t]=spectrogram(sigIn,128,96,1024,fs);
            imagesc(t,f,10*log10(abs(s))); hold on;
            axis xy;
            plot(taxis1,data.fmts(:,1:2),'w-','LineWidth',1);
            ys=get(gca,'ylim');
            [i1,i2,jnk1,jnk2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
            plot(repmat(taxis1(round(iv1+lim1*(iv2-iv1))),1,2),ys,'w-','LineWidth',1);
            plot(repmat(taxis1(round(iv1+lim2*(iv2-iv1))),1,2),ys,'w-','LineWidth',1);    

            plot(taxis1,data.rms(:,1)*1e4,'g-','LineWidth',1);
        % 	set(gca,'XLim',[taxis1(i1),taxis1(i2)]);
            set(gca,'YLim',[0,2500]);
            

            disp(['rmsThresh = ',num2str(idata.rmsThresh(idx))]);
            rmsThresh_new=input('New rmsThresh = ','s');
            bDone=0;
            if isempty(rmsThresh_new)
                bDone=1;
            end

            while ~bDone
                if isequal(rmsThresh_new,'a')   % accept
                    idata.F1(idx)=mean(data.fmts(round(iv1+lim1*(iv2-iv1)):round(iv1+lim2*(iv2-iv1)),1));
                    idata.F2(idx)=mean(data.fmts(round(iv1+lim1*(iv2-iv1)):round(iv1+lim2*(iv2-iv1)),2));
                    idata.tVoiceBegin(idx)=iv1*data.params.frameLen/data.params.sr;
                    idata.tVoiceEnd(idx)=iv2*data.params.frameLen/data.params.sr;
                    idata.voiceDur(idx)=(iv2-iv1)*data.params.frameLen/data.params.sr;
                    idata.rmsThresh(idx)=p.rmsThresh;
                    disp(['New: F1 = ',num2str(idata.F1(idx)),'; F2 = ',num2str(idata.F2(idx)),...
                        '; tVoiceBegin = ',num2str(idata.tVoiceBegin(idx)),'; tVoiceEnd = ',num2str(idata.tVoiceEnd(idx)),...
                        '; voiceDur = ',num2str(idata.voiceDur(idx))]);

                    bDone=1;
                    break;
                elseif isequal(rmsThresh_new,'q')
                    bDone=1;
                    disp('Changes discarded.');
                    break;                
                elseif ~isnan(str2num(rmsThresh_new))
                    p=data.params;
                    p.rmsThresh=str2num(rmsThresh_new);
                    fs=testTSM2(data,p);
                    f1=fs(:,1); f2=fs(:,2);
                    plot(taxis1,fs,'m-');
                    [i1,i2,jnk1,jnk2,iv1,iv2]=getFmtPlotBounds(f1,f2);
                    plot(repmat(taxis1(round(iv1+lim1*(iv2-iv1))),1,2),ys,'m-','LineWidth',1);
                    plot(repmat(taxis1(round(iv1+lim2*(iv2-iv1))),1,2),ys,'m-','LineWidth',1);
                end

                rmsThresh_new=input('New rmsThresh = ','s');
            end

            set(gcf,'CurrentAxes',hs1);
            cla;
            plot(idata.x,idata.(items{n}),'bo-'); hold on;
            set(gca,'XLim',[0,idata.x(end)+1]);
            ylabel(items{n});

            if exist('data')
                clear('data');
            end
            set(gcf,'CurrentAxes',hs1);
            x=ginput(1); x=x(1);
        end
        close(hf1);
    end
    
    eval([data_list{k},'=idata;']);
end

%% Write data to disk
toWrite=1;
if isfile(fullfile(expDir,'pdata.mat'))
    a=input([fullfile(expDir,'pdata.mat'),' already exists. Overwrite? (0 - no, 1 - yes): ']);
    if a==0
        toWrite=0;
    else
        toWrite=1;
    end
end
if toWrite
    if expt.subject.designNum==1
        save(fullfile(expDir,'pdata.mat'),'pdata','mdata');
    else
        save(fullfile(expDir,'pdata.mat'),'pdata');
    end
    disp([fullfile(expDir,'pdata.mat'),' saved.']);
end
return