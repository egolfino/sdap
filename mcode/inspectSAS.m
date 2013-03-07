function inspectSAS(subjID)
%%
hostName=getHostName;

if isequal(hostName,'smcg-w510')
    rawDataDir='E:\DATA\SAS';
elseif isequal(hostName,'glossa')
    rawDataDir='G:\DATA\SAS';
end

%%
expDir=fullfile(rawDataDir,subjID);
if ~isdir(expDir)
    fprintf('ERROR: expDir %s does not exist.\n',expDir);
    return
end

d1=dir(fullfile(expDir,'stay/rep*'));

figure;
h1=subplot(2,1,1);
h2=subplot(2,1,2);
for i1=1:numel(d1)
    repDir=fullfile(expDir,'stay',d1(i1).name);
    d2=dir(fullfile(expDir,'stay',d1(i1).name,'trial-*-1.mat'));
    for i2=1:numel(d2)
        dataFN=fullfile(repDir,d2(i2).name);
        load(dataFN);
        
        
        
        fs=data.params.sr;
        frameDur=data.params.frameLen/data.params.sr;
        taxis1=0:frameDur:(frameDur*(length(data.fmts(:,1))-1));
        sigIn=data.signalIn;
        sigOut=data.signalOut;
        f1v=data.fmts(:,1);
        f2v=data.fmts(:,2);
        f1s=data.sfmts(:,1);
        [j1,j2]=getFmtPlotBounds(f1v,f2v);
        
        set(gcf,'CurrentAxes',h1);
        cla;
        [s,f,t]=spectrogram(sigIn,128,96,1024,fs);
        imagesc(t,f,10*log10(abs(s))); hold on;
        axis xy;
        plot(taxis1,f1v,'w');
        plot(taxis1,f2v,'w');
        plot(taxis1,f1s,'g');
%         set(gca,'XLim',taxis1([j1,j2]),'YLim',[0,3000]);
        set(gca,'YLim',[0,3000]);
        
        title(strrep(dataFN,'\','/'));
        
        set(gcf,'CurrentAxes',h2);
        cla;
        [s,f,t]=spectrogram(sigOut,128,96,1024,fs);
        imagesc(t,f,10*log10(abs(s))); hold on;
        axis xy;
        plot(taxis1,f1v,'w');
        plot(taxis1,f2v,'w');
        plot(taxis1,f1s,'g');
%         set(gca,'XLim',taxis1([j1,j2]),'YLim',[0,3000]);
        set(gca,'YLim',[0,3000]);
        
        drawnow;
        soundsc(data.signalIn,fs);
        soundsc(data.signalOut,fs);
        
        ginput(1);
    end
end
return