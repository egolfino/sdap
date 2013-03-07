function analyzeVowelData2(expDir)
%%
mvaWinWidth=24;

lim1=0.3;
lim2=0.7;

Y_LIM=[0,3000];
MID_LIM1=0.4;
MID_LIM2=0.6;
%%
pdataFN=fullfile(expDir,'pdata_new.mat');

load(fullfile(expDir,'expt.mat'));
load(fullfile(expDir,'micRMS_100dBA.mat'));

trainVowel.F1=struct;
trainVowel.F2=struct;
testVowel.F1=struct;
testVowel.F2=struct;

pdata.MID_LIM1=MID_LIM1;
pdata.MID_LIM2=MID_LIM2;

pdata.F1=[];
pdata.F2=[];
pdata.F1Mid=[];
pdata.F2Mid=[];
pdata.F1_old=[];
pdata.F2_old=[];
pdata.tVoiceBegin=[];
pdata.tVoiceEnd=[];
pdata.vowelEndForce=[];
pdata.stage={};
pdata.repNum=[];
pdata.trialNum=[];
pdata.voiceDur=[];
pdata.lv=[];
% pdata.stage=cell(1,0);
pdata.file=cell(1,0);
pdata.rmsThresh=[];
pdata.x=[];
pdata.rating=[];    % 0 - 3 integer. 0 means discard; 3 (default) means perfect.
pdata.comments={};
pdata.words={};
pdata.trajF1=[];
pdata.trajF2=[];

pdata.bProdErr=[];

pdata.rmsThresh=[];
pdata.fn1=[];
pdata.fn2=[];
pdata.aFact=[];
pdata.bFact=[];
pdata.gFact=[];

pdata.rmsThresh_orig=[];
pdata.fn1_orig=[];
pdata.fn2_orig=[];
pdata.aFact_orig=[];
pdata.bFact_orig=[];
pdata.gFact_orig=[];

%%
if isfile(pdataFN)
    bContinue=input(sprintf('pdata_new file already exists. Continue? (0/1): '));
    if bContinue==0
        delete(pdataFN);
    else
        load(pdataFN);
    end    
end

if expt.subject.designNum==1
    mdata=pdata;
end

figure('Name','SAS behav data screening','Unit','Normalized','Position',[0.1,0.3,0.8,0.5],'ToolBar','none');
h1=subplot('Position',[0.05,0.1,0.9,0.85]);

persistent_rmsThresh=NaN;

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
            dataFN=fullfile(expDir,stg,['rep',num2str(n)],d(k).name);
            
            if ~isempty(fsic(pdata.file,dataFN))
                fprintf('%s already processed. Skipped.\n',dataFN)
                continue;
            end
            
			load(dataFN);
            dataOrig=data;
                        
            this_utter.rmsThresh_orig=data.params.rmsThresh;
            if isnan(persistent_rmsThresh)
                this_utter.rmsThresh=data.params.rmsThresh;
            else
                this_utter.rmsThresh=data.params.rmsThresh;
%                 this_utter.rmsThresh=persistent_rmsThresh;
            end
            this_utter.fn1_orig=data.params.fn1;
            this_utter.fn1=data.params.fn1;
            this_utter.fn2_orig=data.params.fn2;
            this_utter.fn2=data.params.fn2;
            this_utter.aFact_orig=data.params.aFact;
            this_utter.aFact=data.params.aFact;
            this_utter.bFact_orig=data.params.bFact;
            this_utter.bFact=data.params.bFact;
            this_utter.gFact_orig=data.params.gFact;
            this_utter.gFact=data.params.gFact;            
            
            this_utter.stage=stg;
            this_utter.file=fullfile(expDir,stg,['rep',num2str(n)],d(k).name);
            this_utter.repNum=n;
            this_utter.trialNum=str2num(strrep(strrep(d(k).name,'trial-',''),'-1.mat',''));
            this_utter.word=data.params.name;
                       
            
            [i1,i2,jnk1,jnk2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
            if isempty(iv1) || isempty(iv2)
                fprintf('\tWarning: no vowel was detected.\n');
                iv1=NaN;
                iv2=NaN;               
                
                pdata.bProdErr(end+1)=1;
                
                pdata.F1(end+1)=NaN;
                pdata.F2(end+1)=NaN;
                pdata.F1Mid(end+1)=NaN;
                pdata.F2Mid(end+1)=NaN;
                pdata.F1_old(end+1)=NaN;
                pdata.F2_old(end+1)=NaN;
                pdata.tVoiceBegin(end+1)=NaN;
                pdata.tVoiceEnd(end+1)=NaN;
                pdata.vowelEndForce(end+1)=NaN;
                pdata.trajF1{end+1}=[];
                pdata.trajF2{end+1}=[];

                pdata.voiceDur(end+1)=NaN;
                pdata.lv(end+1)=NaN;

                pdata.stage{end+1}=this_utter.stage;
                pdata.repNum(end+1)=this_utter.repNum;
                pdata.trialNum(end+1)=this_utter.trialNum;
                pdata.file{end+1}=this_utter.file;
                pdata.words{end+1}=this_utter.word;

                pdata.rating(end+1)=NaN;
                pdata.comments{end+1}='';

                pdata.rmsThresh(end+1)=this_utter.rmsThresh;
                pdata.fn1(end+1)=this_utter.fn1;
                pdata.fn2(end+1)=this_utter.fn2;
                pdata.aFact(end+1)=this_utter.aFact;
                pdata.bFact(end+1)=this_utter.bFact;
                pdata.gFact(end+1)=this_utter.gFact;
                
                pdata.rmsThresh_orig(end+1)=this_utter.rmsThresh_orig;
                pdata.fn1_orig(end+1)=this_utter.fn1_orig;
                pdata.fn2_orig(end+1)=this_utter.fn2_orig;
                pdata.aFact_orig(end+1)=this_utter.aFact_orig;
                pdata.bFact_orig(end+1)=this_utter.bFact_orig;
                pdata.gFact_orig(end+1)=this_utter.gFact_orig;
                
                save(pdataFN,'pdata');
                
                continue;
            end
            
            this_utter.F1_old=[trainVowel.F1.(stg),mean(data.fmts(round(iv1+lim1*(iv2-iv1)):round(iv1+lim2*(iv2-iv1)),1))];
            this_utter.F2_old=[trainVowel.F2.(stg),mean(data.fmts(round(iv1+lim1*(iv2-iv1)):round(iv1+lim2*(iv2-iv1)),2))];
            
            this_utter.rating=3;
            this_utter.ratingComments='';
            this_utter.vowelEndForce=0;
            
            toRepeat=1;
            while toRepeat==1
                data=reprocAPSTVData(dataOrig,'rmsThresh',this_utter.rmsThresh,'fn1',this_utter.fn1,'fn2',this_utter.fn2,...
                    'aFact',this_utter.aFact,'bFact',this_utter.bFact,'gFact',this_utter.gFact);

%                 disp(['Loading data: ',fullfile(expDir,stg,['rep',num2str(n)],d(k).name)]);
                fn1=strrep(d(k).name,'.mat','');
                trialType=str2num(fn1(end));

                [i1,i2,jnk1,jnk2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));

                

                f1v=data.fmts(:,1);
                f2v=data.fmts(:,2);
                
                
                frameDur=data.params.frameLen/data.params.sr;
                taxis1=0:(frameDur):(frameDur*(length(f1v)-1));

                [s,f,t]=spectrogram(data.signalIn,128,96,1024,data.params.sr);
                st_spec=[];
                for k1=iv1:iv2
                    [foo,idx_t]=min(abs(t-taxis1(k1)));
                    [foo,idx_f1]=min(abs(f-f1v(k1)));
                    [foo,idx_f2]=min(abs(f-f2v(k1)));

            %         st_spec(end+1)=sqrt(abs(s(idx_f1,idx_t))^2+abs(s(idx_f2,idx_t))^2);
                    st_spec(end+1)=log10(abs(s(idx_f2,idx_t)));
                end
                st_spec=mva(st_spec,mvaWinWidth*2,'hamming');
                [foo,idx_end]=min(diff(st_spec));

                if this_utter.vowelEndForce==0
                    vowelEnd=iv1+idx_end;
                else
                    vowelEnd=this_utter.vowelEnd;
                end
%                 this_utter.vowelEnd=vowelEnd;
                
                this_utter.trajF1=f1v(iv1:vowelEnd);
                this_utter.trajF2=f2v(iv1:vowelEnd);
                
                rms1=calcAWeightedRMS(data.signalIn(iv1*data.params.frameLen:vowelEnd*data.params.frameLen-1),data.params.sr);
                
                this_utter.F1=mean(f1v(iv1:vowelEnd));
                this_utter.F2=mean(f2v(iv1:vowelEnd));
                
                mid_1=round(iv1+(vowelEnd-iv1)*MID_LIM1);
                mid_2=round(iv1+(vowelEnd-iv1)*MID_LIM2);
                this_utter.F1Mid=mean(f1v(mid_1:mid_2));
                this_utter.F2Mid=mean(f2v(mid_1:mid_2));

                cla;
                imagesc(t,f,10*log10(abs(s))); hold on;
                axis xy;            

                plot(taxis1,f1v,'w-','LineWidth',1);
                hold on;
                plot(taxis1,f2v,'w-','LineWidth',1);

                ylim=get(gca,'YLim');
                plot(repmat(taxis1(iv1),1,2),ylim,'b--','LineWidth',1);
                plot(repmat(taxis1(iv2),1,2),ylim,'b--','LineWidth',1);
                plot(repmat(taxis1(vowelEnd),1,2),ylim,'b-','LineWidth',1);
                
                plot(repmat(taxis1(mid_1),1,2),ylim,'y-','LineWidth',1);
                plot(repmat(taxis1(mid_2),1,2),ylim,'y-','LineWidth',1);
                
                xlim=[taxis1(iv1)-0.2,taxis1(iv2)+0.2];
                
                set(gca,'XLim',xlim,'YLim',Y_LIM);
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                
                title(sprintf('%s - rep #%d - trial #%d: %s',this_utter.stage,this_utter.repNum,this_utter.trialNum,this_utter.word));
                
                % Buttons
                set(gcf,'CurrentAxes',h1);
                xs=get(gca,'XLim'); ys=get(gca,'YLim');        
                x1=xs(1); x2=xs(2); y1=ys(1); y2=ys(2);
                rx=range(xs); ry=range(ys);
                
                cmd.items={'Play sigIn','Play sigOut','vowelEnd','rmsThresh'...
                    'fn1','fn2','aFact','bFact','gFact','Rating'};
                cmd.x_left=repmat(x1+0.9*rx,1,length(cmd.items));
                cmd.y_bottom=y2-0.05*(1:length(cmd.items))*ry;
                cmd.x_width=repmat(0.1*rx,1,length(cmd.items));
                cmd.y_height=repmat(0.05*ry,1,length(cmd.items));
                
                for i3=1:length(cmd.items)
                    rectangle('Position',[cmd.x_left(i3),cmd.y_bottom(i3),cmd.x_width(i3),cmd.y_height(i3)],'FaceColor','w','EdgeColor','k');
                    text(cmd.x_left(i3)+0.01*rx,cmd.y_bottom(i3)+0.025*ry,cmd.items{i3},'Color','k');
                end
                % ~Butttons
                
                coord=ginput(1);
                inBox=0;
                for i3=1:length(cmd.items)
                    if (coord(1)>=cmd.x_left(i3) && coord(1)<=cmd.x_left(i3)+cmd.x_width(i3) && coord(2)>=cmd.y_bottom(i3) && coord(2)<=cmd.y_bottom(i3)+cmd.y_height(i3))
                        inBox=i3;
                        break;
                    end
                end
                
                if inBox==0
                    xs=get(gca,'XLim'); ys=get(gca,'YLim');
                    if ((coord(1)<xs(1) || coord(1)>xs(2)) || (coord(2)<ys(1) || coord(2)>ys(2)))
                        toRepeat=0;
                    end
                else
                    if isequal(cmd.items{inBox},'Play sigIn')
                        soundsc(dataOrig.signalIn,data.params.sr);
                    elseif isequal(cmd.items{inBox},'Play sigOut')
                        soundsc(dataOrig.signalOut,data.params.sr);
                    elseif isequal(cmd.items{inBox},'vowelEnd')            
                        title('Set vowelEnd','Color','b','FontSize',13,'FontWeight','Bold');
                        pnt=ginput(1);
                        this_utter.vowelEnd=round(pnt(1)/frameDur);
                        this_utter.vowelEndForce=1;
                    elseif isequal(cmd.items{inBox},'rmsThresh')
                        this_utter.rmsThresh=input(sprintf('[rmsThresh_orig = %.4f; rmsThresh_curr = %.4f] rmsThresh = ',...
                            this_utter.rmsThresh_orig,this_utter.rmsThresh));
                    elseif isequal(cmd.items{inBox},'fn1')
                        this_utter.fn1=input(sprintf('[fn1_orig = %.1f; fn1_curr = %.1f] fn1 = ',...
                            this_utter.fn1_orig,this_utter.fn1));
                    elseif isequal(cmd.items{inBox},'fn2')
                        this_utter.fn2=input(sprintf('[fn2_orig = %.4f; fn2_curr = %.4f] fn2 = ',...
                            this_utter.fn2_orig,this_utter.fn2));
                    elseif isequal(cmd.items{inBox},'aFact')
                        this_utter.aFact=input(sprintf('[aFact_orig = %.4f; aFact_curr = %.4f] aFact = ',...
                            this_utter.aFact_orig,this_utter.aFact));
                    elseif isequal(cmd.items{inBox},'bFact')
                        this_utter.bFact=input(sprintf('[bFact_orig = %.4f; bFact_curr = %.4f] bFact = ',...
                            this_utter.bFact_orig,this_utter.bFact));
                    elseif isequal(cmd.items{inBox},'gFact')
                        this_utter.gFact=input(sprintf('[gFact_orig = %.4f; gFact_curr = %.4f] gFact = ',...
                            this_utter.gFact_orig,this_utter.gFact));
                    elseif isequal(cmd.items{inBox},'Rating')
                            this_utter.rating=input(['[Old rating = ',num2str(this_utter.rating),'] New rating = ']);
                            this_utter.ratingComments=input(['[e.g., (i1) = problematic i1.] Rating comments = '],'s');
%                     elseif isequal(cmd.items{inBox},'Discard utter')
%                         if this_utter.bDiscard==0
%                             title('Utter discarded','Color','r','FontSize',13,'FontWeight','Bold');
%                             this_utter.bDiscard=1;
%                         else
%                             title('Utter un-discarded','Color','g','FontSize',13,'FontWeight','Bold');
%                             this_utter.bDiscard=0;
%                         end                        
%                         pause(0.5);
%                     elseif isequal(cmd.items{inBox},'Save')
%                         save(dacacheFN,'pdata');
%                         fprintf('%s saved\n',dacacheFN);                 
                    end
                end
            end
            
            pdata.bProdErr(end+1)=0;
            
            pdata.F1(end+1)=this_utter.F1;
            pdata.F2(end+1)=this_utter.F2;
            pdata.F1Mid(end+1)=this_utter.F1Mid;
            pdata.F2Mid(end+1)=this_utter.F2Mid;
            pdata.F1_old(end+1)=this_utter.F1_old;
            pdata.F2_old(end+1)=this_utter.F2_old;
            pdata.tVoiceBegin(end+1)=iv1;
            pdata.tVoiceEnd(end+1)=vowelEnd;
            pdata.vowelEndForce(end+1)=this_utter.vowelEndForce;
            pdata.trajF1{end+1}=this_utter.trajF1;
            pdata.trajF2{end+1}=this_utter.trajF2;
            
            pdata.voiceDur(end+1)=frameDur*(vowelEnd-iv1);
            pdata.lv(end+1)=100+20*log10(rms1/micRMS_100dBA);
            
            pdata.stage{end+1}=this_utter.stage;
            pdata.repNum(end+1)=this_utter.repNum;
            pdata.trialNum(end+1)=this_utter.trialNum;
            pdata.file{end+1}=this_utter.file;
            pdata.words{end+1}=this_utter.word;
            
            pdata.rating(end+1)=this_utter.rating;
            pdata.comments{end+1}=this_utter.ratingComments;
            
            pdata.rmsThresh(end+1)=this_utter.rmsThresh;
            pdata.fn1(end+1)=this_utter.fn1;
            pdata.fn2(end+1)=this_utter.fn2;
            pdata.aFact(end+1)=this_utter.aFact;
            pdata.bFact(end+1)=this_utter.bFact;
            pdata.gFact(end+1)=this_utter.gFact;
            
            pdata.rmsThresh_orig(end+1)=this_utter.rmsThresh_orig;
            pdata.fn1_orig(end+1)=this_utter.fn1_orig;
            pdata.fn2_orig(end+1)=this_utter.fn2_orig;
            pdata.aFact_orig(end+1)=this_utter.aFact_orig;
            pdata.bFact_orig(end+1)=this_utter.bFact_orig;
            pdata.gFact_orig(end+1)=this_utter.gFact_orig;
            
            persistent_rmsThresh=this_utter.rmsThresh;
            
            save(pdataFN,'pdata');
        end
        
	end
end

%% Sort the trials according to the real order

order_code=nan(size(pdata.file));
nRepTrials=numel(expt.script.start.rep1.trialOrder);
stages={'start','ramp','stay','end'};

for i1=1:length(pdata.file)
    order_code(i1)=1000*fsic(stages,pdata.stage{i1})+(pdata.repNum(i1)-1)*nRepTrials+pdata.trialNum(i1);
end
[foo,idx_sort]=sort(order_code);

flds=fields(pdata);
for i1=1:numel(flds)
    fld=flds{i1};
    if numel(pdata.(fld))==numel(pdata.file)
        pdata.(fld)=pdata.(fld)(idx_sort);
    end
end

save(pdataFN,'pdata');
fprintf('pdata_new file %s saved.\n',pdataFN);
return