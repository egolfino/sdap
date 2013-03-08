function [dataOut, bRmsGood, bSpeedGood] = checkData(data,handles)
% This function is called after stoprec is executed. It checks the data and
% displays the rms and transition length . If the rms and speed are in range
%, the data is stored in handles.dataOut.

t0 = 1/data.params.sr;
taxis = 0:t0:t0*(length(data.signalIn)-1);

dBrange = handles.dBRange; %SC-Mod(2008/01/05). One side tolerance: 0.4*dBrange 
rmsval = 50;
speedval = 0;
dataOut = data;

frameDur = data.params.frameLen / data.params.sr;

% if (handles.bMon==1)
% [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
% [k1,k2]=detectVowel(data.fmts(:,1),data.fmts(:,2),iv1,iv2,'eh','rms',data.rms(:,1));
if ~isempty(data)
	vocaLen = handles.vocaLen;
	lenRange = handles.lenRange;
	rmsTransTarg = handles.rmsTransTarg;

    idxAbvThr = find(data.fmts(:, 1));
    vocaLenNow = length(idxAbvThr); 
% 	t1=k1*data.params.frameLen;
% 	t2=k2*data.params.frameLen;

% 	vocaLenNow = k2 - k1 + 1;   %SC-Mod(2008/04/06)

	%SC-Mod(2008/04/06): Look at the rms during the transition, instead of
	%   during the entire vocal part.
% 	if (isnan(t1) | isnan(t2) | isempty(t1) | isempty(t2) | t1>=t2)
% 		rmsTrans=0;
% 		rmsBGNoise=0;
%         if ~isempty(data.signalIn)
%             rmsBGNoise=calcAWeightedRMS(data.signalIn(1:round(0.2*data.params.sr)),data.params.sr);
%         end
% 	else
% % 		rmsTrans=sqrt(mean(data.signalIn(t1:t2).^2));
% % 		rmsBGNoise=sqrt(mean(data.signalIn(1:round(0.2*data.params.sr)).^2));
% 		rmsTrans=calcAWeightedRMS(data.signalIn(t1:t2),data.params.sr);
% 		rmsBGNoise=calcAWeightedRMS(data.signalIn(1:round(0.2*data.params.sr)),data.params.sr);
% 	end

    if isempty(idxAbvThr)
        rmsTrans = 0;
    else
        rmsTrans = sqrt(mean(data.rms(idxAbvThr, 1) .^ 2));
    end
    
    if ~isempty(data.signalIn)
        rmsBGNoise = sqrt(mean(data.signalIn(1 : round(0.2*data.params.sr)).^2));
    else
        rmsBGNoise = 0;   
    end
    

	rmsval   = round(100/dBrange*max(0,min(dBrange,dBrange/2+10*log10(rmsTrans/rmsTransTarg))));   %SC-Mod(2007/12/29)
	speedval = round(100/lenRange*max(0,min(lenRange,lenRange/2+(vocaLen-vocaLenNow)/2)));   
end
%--------------------------------------------------------------------------
%SC Set the volume/speed indicator
if (handles.vumeterMode==1)
    vumeter=handles.skin.vumeter;
elseif (handles.vumeterMode==2);
    vumeter=handles.skin.vumeter2;
end
vumeter0=vumeter*0.5;
% vubounds=handles.skin.vubounds; %SC(12/11/2007)


mask=zeros(size(vumeter));
mask0=zeros(size(vumeter));

rmsval1=floor(rmsval/10)*10;
if (rmsval1+10>size(vumeter,2))
	rmsval1=size(vumeter,2)-10;
end
if (rmsval1<0)
	rmsval1=0;
end
if (handles.trialType==3 | handles.trialType==4 | handles.trialType==5)
	rmsval1=40+rand*10;
	if (rmsval1>45) rmsval1=50;
	else rmsval1=40;
	end
end

if (handles.vumeterMode==1)
    mask(:,rmsval1+1:rmsval1+10,:) = 1;   %SC-Commented(12/11/2007)
    mask0=1-mask;
elseif (handles.vumeterMode==2)
    if (rmsval1<30)
        mask(:,1:30,:) = 1;
    elseif (rmsval1>=30 && rmsval1<70)
        mask(:,31:70,:) = 1;
    else
        mask(:,70:100,:) = 1;                
    end
    mask0=1-mask;    
end

set(handles.rms_imgh, 'Cdata', vumeter .* mask + vumeter0 .* mask0);

mask = zeros(size(vumeter));
mask0 = zeros(size(vumeter));

speedval1=floor(speedval/10)*10;
if (speedval1+10>size(vumeter,2))
	speedval1=size(vumeter,2)-10;
end
if (speedval1<0)
	speedval1=0;
end
% if (handles.trialType==3 || handles.trialType==4 || handles.trialType==5)
% 	speedval1=40+rand*10;
% 	if (speedval1>45) speedval1=50;
% 	else speedval1=40;
% 	end	
% end

if (handles.vumeterMode==1)
    mask(:,speedval1+1:speedval1+10,:) = 1;         %SC(2008/01/05)
    mask0=1-mask;
elseif (handles.vumeterMode==2)
    if (speedval1<30)
        mask(:,1:30,:) = 1;
    elseif (speedval1>=30 && speedval1<70)
        mask(:,31:70,:) = 1;
    else
        mask(:,70:100,:) = 1;          
    end
    mask0=1-mask;
end

% if (handles.trialType==1 | handles.trialType==2)
set(handles.speed_imgh, 'Cdata', vumeter.*mask + vumeter0.*mask0);
    
% end
%SC ~Set the volume/speed indicator
%--------------------------------------------------------------------------
drawnow

msg1='';    msg2='';
instr1='';  instr2='';
if (rmsval < 70 && rmsval > 30)
	bRmsGood=1;
else
	bRmsGood=0;
	if (rmsval >= 70)   %SC(2008/01/05)
		msg1='Softer';
		instr2='Loud';
	else % Then rmsval <= 30
		msg1='Louder';
		instr2='Soft';
	end
end

if (speedval < 70 && speedval > 30) %SC-Mod(2008/01/05) Used to be speedval > 20
	bSpeedGood=1;
else
	bSpeedGood=0;
	if (speedval >= 70) %SC (2008/01/05)
		msg2='Slower';
		instr1='Fast';
	else
		msg2='Faster';
		instr1='Slow';
	end
end

% if (handles.trialType==3 | handles.trialType==4 | handles.trialType==5)
% 	bRmsGood=1;
% 	bSpeedGood=1;
% end

%SC(2008/01/05)
if (~bRmsGood || ~bSpeedGood)
	if (~bRmsGood && ~bSpeedGood)
%         msgc=['Please speak ',msg1,' and ',msg2,'.'];
		msgc=[msg1,' and ',lower(msg2),' please!'];
	elseif (~bRmsGood)
%         msgc=['Please speak ',msg1,'.'];
		msgc=[msg1,' please!'];        
	elseif (~bSpeedGood)
		msgc=['Please speak ',msg2,'.'];
		msgc=[msg2,' please!'];
	end

	if (handles.showTextCue)
		set(handles.msgh,'string',{'';msgc});
% 		pause(1);
	end
end

% SCai: update the data monitor window
set(0,'CurrentFigure',handles.figIdDat(1));
set(gcf,'CurrentAxes',handles.figIdDat(2));
cla;
plot(taxis,data.signalIn);      hold on;
set(gca,'XLim',[taxis(1);taxis(end)]);
set(gca,'YLim',[-1,1]);
ylabel('Wave In');

set(gcf,'CurrentAxes',handles.figIdDat(3));
cla;
taxis=0:t0:t0*(length(data.signalOut)-1);
plot(taxis,data.signalOut*data.params.dScale);     hold on;
set(gca,'XLim',[taxis(1);taxis(end)]);
set(gca,'YLim',[-1,1]);
xlabel('Time (s)');
ylabel('Wave Out');



% [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
% [k1,k2]=detectVowel(data.fmts(:,1),data.fmts(:,2),iv1,iv2,'eh','rms',data.rms(:,1));
% if (~isnan(i1) && ~isnan(i2) && ~isempty(i1) && ~isempty(i2) && k2 >= k1)
% t1=k1*data.params.frameLen;
% t2=k2*data.params.frameLen;
% tv1=min(find(data.fmts(:,1)>0));
% tv2=max(find(data.fmts(:,1)>0));

% idx1=max([1,tv1-50]);
% idx2=min([tv2+50,length(data.signalIn)]);
% 
% %     wavInGain=0.13827;  % w/Pa
% p0=20e-6;           % Pa
% % 	tRMSIn=sqrt(mean((data.signalIn(t1:t2)).^2));
% 
% set(gcf,'CurrentAxes',handles.figIdDat(2));
% xs=get(gca,'XLim'); ys=get(gca,'YLim');
% 
% if (~isnan(t1) && ~isnan(t2) && t1>0 && t2>0 && t2>t1)
%     plot([taxis(t1),taxis(t1)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
%     plot([taxis(t2),taxis(t2)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
%     tRMSOut=calcAWeightedRMS(data.signalOut(t1:t2),data.params.sr);
% else 
%     tRMSOut=0;
% end
% 
% %         load calibMic;  % gets micGain: wav rms/ Pa rms (Pa^-1)
% load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
% %     text(xs(1)+0.05*range(xs),ys(2)-0.1*range(ys),['RMS(In)=',num2str(tRMSIn)]);
% text(xs(1)+0.05*range(xs),ys(2)-0.2*range(ys),...
%     ['soundLevel=',num2str(100+20*log10((rmsTrans/micRMS_100dBA))),' dBA SPL'],'FontSize',11);
% dataOut.vowelLevel=100+20*log10((rmsTrans/micRMS_100dBA));
% 
% text(xs(1)+0.05*range(xs),ys(2)-0.25*range(ys),...
%     ['BGNoiseLevel=',num2str(100+20*log10((rmsBGNoise/micRMS_100dBA))),' dBA SPL'],'FontSize',11);
% 
% text(xs(1)+0.05*range(xs),ys(2)-0.3*range(ys),...
%     ['SNR=',num2str(20*log10(rmsTrans/rmsBGNoise))],'FontSize',11);
% 
% %         load calibOutput;   
% % gives 'freq' and 'voltGains', measured at 'shanqing' M-audio configuration 
% % and -1.65 Phone volume knob.
% %         mvg=mean(voltGains) * sqrt(2);    % mean voltage gain (V_rms / wavAmp_rms)
% 
% soundLvOut=20*log10(tRMSOut*data.params.dScale/(dBSPL2WaveAmp(0,1000)/sqrt(2)));  % dBA SPL
% 
% set(gcf,'CurrentAxes',handles.figIdDat(3));
% xs=get(gca,'XLim'); ys=get(gca,'YLim');
% if (~isnan(t1) && ~isnan(t2) && t1>0 && t2>0 && t2>t1)	
%     plot([taxis(t1),taxis(t1)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
%     plot([taxis(t2),taxis(t2)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
% end
% text(xs(1)+0.05*range(xs),ys(2)-0.15*range(ys),...
%     ['dScale=',num2str(data.params.dScale)],'FontSize',11);    
% text(xs(1)+0.05*range(xs),ys(2)-0.2*range(ys),...
%     ['soundLevel=',num2str(soundLvOut),' dBA SPL'],'FontSize',11);
% 
% 
% set(gcf,'CurrentAxes',handles.figIdDat(4));
% cla;
% if (data.params.frameLen*idx1>=1 & data.params.frameLen*idx2<=length(taxis) & idx1>=1 & idx2 <= size(data.fmts,1))
%     plot(taxis(data.params.frameLen*(idx1:idx2)),data.fmts(idx1:idx2,1),'k-','LineWidth',1.5);   hold on;
%     plot(taxis(data.params.frameLen*(idx1:idx2)),data.fmts(idx1:idx2,2),'k-','LineWidth',1.5); 
%     plot(taxis(data.params.frameLen*(idx1:idx2)),data.sfmts(idx1:idx2,1),'b-','LineWidth',1.5);
%     plot(taxis(data.params.frameLen*(idx1:idx2)),data.sfmts(idx1:idx2,2),'b-','LineWidth',1.5);
%     set(gca,'XLim',taxis(data.params.frameLen*([idx1,idx2])));
%     set(gca,'YLim',[0,3000]);
%     xs=get(gca,'XLim'); ys=get(gca,'YLim');
%     plot(taxis(data.params.frameLen*([k1,k1])),[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
%     plot(taxis(data.params.frameLen*([k2,k2])),[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;    
%     xlabel('Time (sec)');
%     ylabel('Formant freqs (Hz)');
% else
%     cla;
% end
% 
% set(gcf,'CurrentAxes',handles.figIdDat(5));
% cla;
% if (~isnan(k1) && ~isnan(k2) && k1>0 && k2>0 && t2>t1)
%     plot(data.fmts(k1:k2,1),data.fmts(k1:k2,2),'b-','LineWidth',1.5);   hold on;
%     plot(data.sfmts(k1:k2,1),data.sfmts(k1:k2,2),'b-','LineWidth',1.5);   hold off;
% end
% set(gca,'XLim',[0,2000]);
% set(gca,'YLim',[0,3000]);
% grid on;
% xlabel('F1 (Hz)');
% ylabel('F2 (Hz)');
% 
% set(gcf,'CurrentAxes',handles.figIdDat(6));
% cla;
% 
% set(gcf,'CurrentAxes',handles.figIdDat(7));
% cla;
% if (~isnan(k1) && ~isnan(k2) && k1>0 && k2>0 && t2>t1)
%     plot(hz2mel(data.fmts(k1:k2,1)),hz2mel(data.fmts(k1:k2,2)),'b-','LineWidth',1.5);   hold on;
%     plot(hz2mel(data.sfmts(k1:k2,1)),hz2mel(data.sfmts(k1:k2,2)),'b-','LineWidth',1.5);   hold off;
% end
% set(gca,'XLim',[0,1000]);
% set(gca,'YLim',[0,2000]);
% grid on;	
% xlabel('F1 (mel)');
% ylabel('F2 (mel)');    

% else
% 	pause(1e-3);
% end
% ~SCai: update the data monitor window
% --------------------------------------------------------------------------

% handles.time2 = clock;
% if (handles.trigByScanner == 0)
% 	timeToPause=handles.ITI-etime(handles.time2,handles.time1)-0.1; % 0.1 is the safety margin
% 	if (timeToPause>0)
% 		pause(timeToPause);
%     end
% else
%     pause(0.25);
% end

% Spectrogram and RMS of input 
set(gcf,'CurrentAxes',handles.figIdDat(4));
cla;
[s, f, t]=spectrogram(data.signalIn, 128, 96, 1024, data.params.sr);
imagesc(t, f, 10 * log10(abs(s))); hold on;
axis xy;

tAxis1 = 0 : frameDur : frameDur * (length(data.rms(:, 1)) - 1);
plot(tAxis1, data.rms(:, 1) / max(data.rms(:, 1)) * 2e3, 'b-', 'LineWidth', 1.5);
plot([tAxis1(1), tAxis1(end)], repmat(data.params.rmsThresh / max(data.rms(:, 1)) * 2e3, 1, 2), 'b--', 'LineWidth', 1.5);

set(gca, 'XLim', [t(1), t(end)]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Spectrogram of output
set(gcf,'CurrentAxes',handles.figIdDat(5));
cla;
[s, f, t]=spectrogram(data.signalOut, 128, 96, 1024, data.params.sr);
imagesc(t, f, 10 * log10(abs(s))); hold on;
axis xy;
    
set(gca, 'XLim', [t(1), t(end)]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');


set(handles.msgh, 'string', '');
if (handles.vumeterMode==1)
    vumeter=handles.skin.vumeter;
elseif (handles.vumeterMode==2)
    vumeter=handles.skin.vumeter2;
end
mask=0.5*ones(size(vumeter));
set(handles.rms_imgh, 'Cdata', vumeter.*mask);
set(handles.speed_imgh, 'Cdata', vumeter.*mask);