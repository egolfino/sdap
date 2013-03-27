function runExperiment(varargin)
%%
TRIAL_TYPES = {'NORMAL', 'DOWN_SHIFT', 'UP_SHIFT', 'MASKED', 'BASELINE'};

%%

DEBUG=0;
[ret,hostName]=system('hostname');

expt_config_fn = '../expt/expt_config.txt';
fprintf('Reading experiment configuration from "%s" ... ', expt_config_fn);
expt_config=read_parse_expt_config(expt_config_fn);
fprintf('Done.\n');

nRuns = expt_config.NUM_RUNS;
nTrialsPerRun = expt_config.TRIALS_PER_RUN;

expt_design_fn = '../expt/expt_design.xls';
fprintf('Reading experiment design from "%s" ... ', expt_design_fn);
expt_design = read_expt_design(expt_design_fn, nRuns, nTrialsPerRun);
fprintf('Done.\n');

%% 
expt.expt_config         = expt_config; 
expt.expt_design         = expt_design;
% expt.name				= 'TS_20101029_1'; 
% expt.sex					= 'male';  % male / female
% expt.shiftDirection		= 'F1Down';  %SC F1Up / F1Down
% expt.shiftRatio			= 0.3;    %
% expt.mouthMicDist        =6;   % cm

expt.closedLoopGain      = 14 - 20 * log10(10 / expt.expt_config.MOUTH_MIC_DIST);

% expt.dBRange1            =20;        % is the one-sided dBRange1*0.4
% expt.dBRange2            =12;         % Tightened level range after the initial pract1 training. 
% expt.trialLen            =2.2;

expt.hostName            =deblank(hostName);
% if isequal(expt.hostName,'smcg-w510')
%     expt.dataDir             ='E:\DATA\APSTV2';
% else
    % expt.dataDir				='C:\CS_2004\PROJECTS\SAP-FMRI\';
dataDir				= expt.expt_config.DATA_DIR;
% end

expt.trigByScanner		=1;
TR                  = expt.expt_config.MRI_TR;
TA                  = expt.expt_config.MRI_TA;
stimDelay           = expt.expt_config.STIM_DELAY;
recDelay           = expt.expt_config.REC_DELAY;

expt.vumeterMode         =2;     % 1: 10 ticks; 2: 3 ticks;

if expt.trigByScanner == 1
	expt.showProgress		= 0;
	expt.showPlayButton      = 0;
else
	expt.showProgress		= 1;
	expt.showPlayButton      = 1;
end

% expt.designNum			= 2;
% expt.lvNoise             = 75; % dBA SPL. The level of noise for completely masking speech (mode trialType = 2 or 3).

%%
expt.date				= clock;

% if (~isempty(findStringInCell(varargin, 'subject')))
%     clear('subject');
%     subject=varargin{findStringInCell(varargin,'subject')+1};
% end

% if (~isfield(subject,'pcrKnob'))
%     subject.pcrKnob=input('Phone / Ctrl Room knob = ');
% end
expt.pcrKnob=0;

%%

bNew=true;

dirname=fullfile(dataDir, num2str(expt.expt_config.SUBJECT_ID));

if (~isempty(findStringInCell(varargin,'dirname')))
    clear('dirname');
    dirname=varargin{findStringInCell(varargin,'dirname')+1};
end

if isdir(dirname)
    messg={'The specified directory already contains a previously recorded experiment'
        ''
        'Continue experiment, overwrite  or cancel ?'};
    button1 = questdlg(messg,'DIRECTORY NOT EMPTY','Continue','Overwrite','Cancel','Continue');
    switch button1
        case 'Overwrite'
            button2 = questdlg({'Are you sure you want to overwrite experiment'} ,'OVERWRITE EXPERIMENT ?');
            switch button2
                case 'Yes',
                    rmdir(dirname,'s')
                otherwise,
                    return
            end
        case 'Continue'
            bNew=false;

        otherwise,
            return

    end
end


if bNew % set up new experiment
    mkdir(dirname)
    
    p = getTSMDefaultParams(expt.expt_config.SUBJECT_GENDER, ...
                            'DOWNSAMP_FACT', expt.expt_config.DOWNSAMP_FACT, ...
                            'FRAME_SIZE', expt.expt_config.FRAME_SIZE / expt.expt_config.DOWNSAMP_FACT);
    state.run = 1;
%     state.trial = 1;
    state.params = p;
    rmsPeaks = [];
    
    save(fullfile(dirname, 'expt.mat'), 'expt');
    save(fullfile(dirname, 'state.mat'), 'state');
else % load expt
    load(fullfile(dirname, 'state.mat'));
    load(fullfile(dirname, 'expt.mat'));            
    p = state.params;
%     nPeaks = length(expt.trainWords);
%     if state.phase>1
%         rmsPeaks=ones(length(expt.trainWords),1)*p.rmsMeanPeak;    %SC ***Bug!!***
%     end
%     subject = expt.subject;
end


%% initialize algorithm
MexIO('init',p);      %SC Set the initial (default) parameters

TransShiftMex(1);
pause(0.5);
TransShiftMex(2);
[~] = TransShiftMex(4);

if ((p.frameShift-round(p.frameShift)~=0) || (p.frameShift>p.frameLen))
    uiwait(errordlg(['Frameshift = ' num2str(p.frameShift) ' is a bad value. Set nWin and frameLen appropriately. Frameshift must be an integer & Frameshift <= Framelen'],'!! Error !!'))
    return
else
    fprintf('\n  \n')
    TransShiftMex(0);           %SC Gives input/output device info, and serves as an initialization.

    fprintf('\nSettings : \n')
    fprintf('DMA Buffer    = %i samples \n',p.frameLen) %SC Buffer length after downsampling
    fprintf('Samplerate    = %4.2f kHz \n',p.sr/1000)   %SC sampling rate after downsampling
    fprintf('Analysis win  = %4.2f msec \n',p.bufLen/p.sr*1000)
    fprintf('LPC  window   = %4.2f msec \n',p.anaLen/p.sr*1000)

    fprintf('Process delay = %4.2f msec \n',p.nDelay*p.frameLen/p.sr*1000)
    fprintf('Process/sec   = %4.2f \n',p.sr/p.frameShift)

end

%% Load the multi-talker babble noise
[x_mtb, fs_mtb]=wavread('mtbabble48k.wav');
% lenMTB=round(2.5*fs_mtb);

gainMTB_fb=dBSPL2WaveAmp(expt.expt_config.MASK_NOISE_LV, 1000)/sqrt(2)/calcMaskNoiseRMS;
% gainMTB_fb3=dBSPL2WaveAmp(subject.lvNoise3,1000,subject.pcrKnob)/sqrt(2)/calcMaskNoiseRMS;
% x_mtb=cell(1,3);
% x_mtb{1}=x(1:lenMTB);               x_mtb{1}=x_mtb{1}-mean(x_mtb{1});
% x_mtb{2}=x(lenMTB+1:lenMTB*2);      x_mtb{2}=x_mtb{2}-mean(x_mtb{2});
% x_mtb{3}=x(lenMTB*2+1:lenMTB*3);    x_mtb{3}=x_mtb{3}-mean(x_mtb{3});
% TransShiftMex(3,'datapb',x_mtb{1});

%% expt
      
figIdDat=makeFigDataMon;

% wordList=expt.words;

% allPhases = expt.allPhases;
% recPhases = expt.recPhases;
% nWords=length(wordList);

hgui=UIRecorder('figIdDat',figIdDat);

% if (expt.subject.designNum == 2)
%     expt.script=addFaceInfo(expt.script,hgui.skin.dFaces);
%     expt.dFaces=hgui.skin.dFaces;
% end

hgui.pcrKnob = expt.pcrKnob;
hgui.TR = expt.expt_config.MRI_TR;
hgui.TA = expt.expt_config.MRI_TA;
hgui.stimDelay = expt.expt_config.STIM_DELAY;
hgui.recDelay = expt.expt_config.REC_DELAY;

hgui.exptDir = dirname;

hgui.trigByScanner = 1;

hgui.dBRange = expt.expt_config.SPL_RANGE;
hgui.trialLen = expt.expt_config.REC_LEN;
% hgui.skin.faceOrder=randperm(length(hgui.skin.dFaces));
% hgui.skin.facePnt=1;

hgui.vumeterMode = expt.vumeterMode;

hgui.rmsTransTarg_spl = getSPLTarg(expt.expt_config.MOUTH_MIC_DIST);
load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));

fprintf('\n');
disp(['Mouth-microphone distance = ', num2str(expt.expt_config.MOUTH_MIC_DIST), ' cm']);
disp(['hgui.rmsTransTarg_spl = ', num2str(hgui.rmsTransTarg_spl), ' dBA SPL']);
fprintf('\n');

hgui.vocaLen=round(300*p.sr/(p.frameLen*1000)); % 300 ms, 225 frames
hgui.lenRange=round(250*p.sr/(p.frameLen*1000));  % single-sided tolerance range: 0.4*250 = 100 ms
disp(['Vowel duration range: [',num2str(300-0.4*250),',',num2str(300+0.4*250),'] ms.']);

hgui.debug=DEBUG;
% hgui.trigKey='equal';

if isempty(findStringInCell(varargin,'twoScreens'))
% 	set(hgui.UIrecorder,...
% 		'position', [0    5.0000  250.6667   65.8750],...
% 		'toolbar','none');  %SC Set the position of the expt window, partially for the use of multiple monitors.
else
	if (expt.trigByScanner == 1)
		ms=get(0,'MonitorPosition');
		set(hgui.UIrecorder,'Position',[ms(2,1),ms(1,4)-ms(2,4),ms(2,3)-ms(2,1)+1,ms(2,4)+20],'toolbar','none','doublebuffer','on','renderer','painters');
		pos_win=get(hgui.UIrecorder,'Position');
		pos_strh=get(hgui.strh,'Position');
		pos_axes_pic=get(hgui.axes_pic,'Position');
		pos_rms_axes=get(hgui.rms_axes,'Position');
		pos_speed_axes=get(hgui.speed_axes,'Position');
		pos_rms_label=get(hgui.rms_label,'Position');
		pos_rms_too_soft=get(hgui.rms_too_soft,'Position');
		pos_rms_too_loud=get(hgui.rms_too_loud,'Position');
		pos_speed_label=get(hgui.speed_label,'Position');
		pos_speed_too_slow=get(hgui.speed_too_slow,'Position');
		pos_speed_too_fast=get(hgui.speed_too_fast,'Position');
		set(hgui.strh,'Position',[(pos_win(3)-pos_strh(3))/2 + 5, (pos_win(4)-pos_strh(4))/2 - 48, pos_strh(3), pos_strh(4)]);
		set(hgui.axes_pic,'Position',[(pos_win(3)-pos_axes_pic(3))/2,(pos_win(4)-pos_axes_pic(4))/2,pos_axes_pic(3),pos_axes_pic(4)]);
		set(hgui.rms_axes,'Position',[(pos_win(3)-pos_rms_axes(3))/2,pos_rms_axes(2),pos_rms_axes(3),pos_rms_axes(4)]);
		set(hgui.rms_label,'Position',[(pos_win(3)-pos_rms_label(3))/2,pos_rms_label(2),pos_rms_label(3),pos_rms_label(4)]);
		set(hgui.rms_too_soft,'Position',[(pos_win(3)-pos_rms_axes(3))/2,pos_rms_too_soft(2),pos_rms_too_soft(3),pos_rms_too_soft(4)]);
		set(hgui.rms_too_loud,'Position',[(pos_win(3)-pos_rms_axes(3))/2+pos_rms_axes(3)-pos_rms_too_loud(3),pos_rms_too_loud(2),pos_rms_too_loud(3),pos_rms_too_loud(4)]);
		set(hgui.speed_axes,'Position',[(pos_win(3)-pos_speed_axes(3))/2,pos_speed_axes(2),pos_speed_axes(3),pos_speed_axes(4)]);		
		set(hgui.speed_label,'Position',[(pos_win(3)-pos_speed_label(3))/2,pos_speed_label(2),pos_speed_label(3),pos_speed_label(4)]);
		set(hgui.speed_too_slow,'Position',[(pos_win(3)-pos_speed_axes(3))/2,pos_speed_too_slow(2),pos_speed_too_slow(3),pos_speed_too_slow(4)]);
		set(hgui.speed_too_fast,'Position',[(pos_win(3)-pos_speed_axes(3))/2+pos_speed_axes(3)-pos_speed_too_fast(3),pos_speed_too_fast(2),pos_speed_too_fast(3),pos_speed_too_fast(4)]);
        set(hgui.msgh, 'FontSize', 30);
	else
		set(hgui.UIrecorder,'Position',[-1400,180,1254,857],'toolbar','none');
    end
    if isequal(expt_config.BACKGROUND_COLOR, 'BLACK')
        set(hgui.UIrecorder, 'Color', 'k');
        warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        jframe=get(hgui.UIrecorder,'javaframe');
        jIcon=javax.swing.ImageIcon('C:\Users\egolfino\Desktop\newicon.gif');
        jframe.setFigureIcon(jIcon);
        set(hgui.msgh, 'BackgroundColor', 'k');
        set(hgui.strh, 'BackgroundColor', 'k');

        hgui.skin.fixation = imread(fullfile(pwd, 'graphics', 'fixation_blackBG.bmp'));
    end
    
    if isequal(expt_config.TEXT_COLOR, 'WHITE')
        set(hgui.msgh, 'ForegroundColor', 'w');
        set(hgui.strh, 'ForegroundColor', 'w');
    end
    
    
	
end

if (expt.showProgress)
	set(hgui.progress_axes,'visible','on');
	set(hgui.progress_imgh,'visible','on');
	progress_meter=0.5*ones(1,100,3);
	progress_mask=zeros(1,100,3);
	set(hgui.progress_imgh,'Cdata',progress_meter.*progress_mask);        
else
	set(hgui.progress_axes,'visible','off');
	set(hgui.progress_imgh,'visible','off');
end

rProgress=0;
startRun = state.run; %SC For the purpose of resumed experiments
% startRep=state.rep;     %SC For the purpose of resumed experiments
n = startRun;

while n <= expt.expt_config.NUM_RUNS
% for n = startRun : expt.expt_config.NUM_RUNS
    fprintf('Entering run #%d.\n', n);
    str = input('Press Enter to continue, or enter ("run x") to jump to run #x > ', 's');
    if ~isempty(deblank(str))
        n = str2num(strrep(lower(deblank(str)), 'run ', ''));
        fprintf('INFO: jumping to run #%d...\n', n);
        continue;
    end
    
    state.run = n;
    save(fullfile(dirname, 'state.mat'), 'state');
    
%     state.rep = 1;
%     thisphase = allPhases{1, n};
    thisrun = sprintf('run%d', n);
    subdirname=fullfile(dirname, thisrun);
    if ~isdir(subdirname)
        mkdir(subdirname);
    else
        d1 = dir(fullfile(dirname, sprintf('run%d_old*', n)));
        subdirname_cp = fullfile(dirname, sprintf('%s_old%d', thisrun, length(d1) + 1));
        mkdir(subdirname_cp);
        movefile(subdirname, subdirname_cp);
        fprintf('Backed up the content of old directory %s at %s\n', subdirname, subdirname_cp);
        
        mkdir(subdirname);
    end
    
    hgui.run = thisrun;
    
    % Adjust the number of reps
%     if (~isequal(thisphase,'ramp') && ~isequal(thisphase,'stay'))
%         disp(['--- Coming up: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
%             '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
%         nRepsNew=input('(Enter to skip) nRepsNew = ','s');
%         nRepsNew=str2num(nRepsNew);
%         if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
%             expt.script.(thisphase).nReps=nRepsNew;
%             expt.script.(thisphase)=genPhaseScript(thisphase, expt.script.(thisphase).nReps,...
%                 expt.trialTypes,expt.trainWords,expt.testWords,expt.pseudoWords,...
%                 expt.trialOrderRandReps,expt.subject.designNum);
%             disp(['Changed: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
%                 '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
%             save(fullfile(dirname,'expt.mat'),'expt');
%             disp(['Saved ',fullfile(dirname,'expt.mat')]);
%         end
%     elseif isequal(thisphase,'ramp')
%         disp(['--- Coming up: ','ramp','. nReps = ',num2str(expt.script.ramp.nReps),...
%             '; nTrials = ',num2str(expt.script.ramp.nTrials),' ---']);
%         disp(['--- Coming up: ','stay','. nReps = ',num2str(expt.script.stay.nReps),...
%             '; nTrials = ',num2str(expt.script.stay.nTrials),' ---']);
%         disp(['--- Ramp+Stay: nReps = ',num2str(expt.script.ramp.nReps+expt.script.stay.nReps),...
%             '; nTrials = ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials)]);
%         nRepsNew=input('(Enter to skip) Ramp: nRepsNew = ','s');
%         nRepsNew=str2num(nRepsNew);
%         if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
%             expt.script.ramp.nReps=nRepsNew;
%             expt.script.ramp=genPhaseScript('ramp',expt.script.ramp.nReps,...
%                 expt.trialTypes,expt.trainWords,expt.testWords,expt.pseudoWords,...
%                 expt.trialOrderRandReps,expt.subject.designNum);
%             disp(['Changed: ramp. ','nReps = ',num2str(expt.script.ramp.nReps),...
%                 '; nTrials = ',num2str(expt.script.ramp.nTrials),' ---']);
%             save(fullfile(dirname,'expt.mat'),'expt');
%             disp(['Saved ',fullfile(dirname,'expt.mat')]);
%         end
%         nRepsNew=input('(Enter to skip) Stay: nRepsNew = ','s');
%         nRepsNew=str2num(nRepsNew);
%         if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
%             expt.script.stay.nReps=nRepsNew;
%             expt.script.stay=genPhaseScript('stay',expt.script.stay.nReps,...
%                 expt.trialTypes,expt.trainWords,expt.testWords,expt.pseudoWords,...
%                 expt.trialOrderRandReps,expt.subject.designNum);
%             disp(['Changed: stay. ','nReps = ',num2str(expt.script.stay.nReps),...
%                 '; nTrials = ',num2str(expt.script.stay.nTrials),' ---']);
%             save(fullfile(dirname,'expt.mat'),'expt');
%             disp(['Saved ',fullfile(dirname,'expt.mat')]);
%         end
%         disp(['--- Ramp+Stay: nReps = ',num2str(expt.script.ramp.nReps+expt.script.stay.nReps),...
%             '; nTrials = ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials)]);
%     end
    % Adjust the number of reps
    
    %%% TODO: adjust the run number %%%
    
    nTrials = expt.expt_design.(thisrun).nTrials;
%     if ~isequal(thisphase,'stay')
%         phaseTrialCnt=1;
%     end

%     expt.script.(thisphase).startTime=clock;

    if expt.expt_config.USE_SPL_TARGET
        switchWord = 'on';
    else
        switchWord = 'off';
    end
    set(hgui.rms_axes, 'visible', switchWord);
    set(hgui.rms_imgh, 'visible', switchWord);
    set(hgui.rms_label, 'visible', switchWord);
    set(hgui.rms_too_soft, 'visible', switchWord);
    set(hgui.rms_too_loud, 'visible', switchWord);
    
    if expt.expt_config.USE_VOWEL_LEN_TARGET
        switchWord = 'on';
    else
        switchWord = 'off';
    end
    set(hgui.speed_axes, 'visible', switchWord);
    set(hgui.speed_imgh, 'visible', switchWord);
    set(hgui.speed_label, 'visible', switchWord);
    set(hgui.speed_too_slow, 'visible', switchWord);
    set(hgui.speed_too_fast, 'visible', switchWord);
    
    hgui.bSpeedRepeat=0;
    hgui.bRmsRepeat=0;
	
	if (expt.showPlayButton==0)
		set(hgui.play, 'visible', 'off');
	end
    
%     switch(thisphase)
%         case 'pre'
%             set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             p.bDetect=0;
%             p.bShift = 0;   %SC No shift in the pre phase
% 			
%             hgui.bRmsRepeat=0;
%             hgui.bSpeedRepeat=0;			
%             
% % 			set(hgui.rms_axes,'visible','off');
% % 		    set(hgui.rms_imgh,'visible','off');
% % 		    set(hgui.rms_label,'visible','off');
% % 			set(hgui.rms_too_soft,'visible','off');
% % 			set(hgui.rms_too_loud,'visible','off');	
% % 		    set(hgui.speed_axes,'visible','off');
% % 		    set(hgui.speed_imgh,'visible','off');
% % 		    set(hgui.speed_label,'visible','off');
% % 			set(hgui.speed_too_slow,'visible','off');
% % 		    set(hgui.speed_too_fast,'visible','off');
% 
%             if (~isempty(rmsPeaks))
%                 p.rmsMeanPeak=mean(rmsPeaks);
%                 p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
%             end
% 
%             hgui.showTextCue=0;
%         case 'pract1'           
%             set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             if (hgui.vumeterMode==1)
%                 vumeter=hgui.skin.vumeter;
%             elseif (hgui.vumeterMode==2)
%                 vumeter=hgui.skin.vumeter2;
%             end
%             mask=0.5*ones(size(vumeter));
%             % mask(1:50,:,:) = 1;           %SC-Commented(12/11/2007)
%             set(hgui.rms_imgh,'Cdata',vumeter.*mask);
%             p.bDetect=0;
%             p.bShift = 0;       %SC No shift in the practice-1 phase
%             hgui.bRmsRepeat=1;
%             hgui.bSpeedRepeat=0;          
% 		    set(hgui.speed_axes,'visible','off');
% 		    set(hgui.speed_imgh,'visible','off');
% 		    set(hgui.speed_label,'visible','off');
% 			set(hgui.speed_too_slow,'visible','off');
% 		    set(hgui.speed_too_fast,'visible','off');
% 
%             if (~isempty(rmsPeaks))
%                 p.rmsMeanPeak=mean(rmsPeaks);
%                 p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
%             end
% 
%             hgui.showTextCue=1;
%             
%             subjProdLevel=[];
%          case 'pract2'
%             subjProdLevel=subjProdLevel(find(~isnan(subjProdLevel)));
%             
%             if (~isempty(subjProdLevel))
%                 hgui.rmsTransTarg_spl=mean(subjProdLevel);
%                 load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
%                 hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));
%             end
%             
%             fprintf('\n');
%             disp(['Target level set as subject mean production level: ',num2str(hgui.rmsTransTarg_spl),' dBA SPL']);
%             fprintf('\n');            
%              
%             set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             if (hgui.vumeterMode==1)
%                 vumeter=hgui.skin.vumeter;
%             elseif (hgui.vumeterMode==2)
%                 vumeter=hgui.skin.vumeter2;
%             end
%             mask=0.5*ones(size(vumeter));
%             set(hgui.speed_imgh,'Cdata',vumeter.*mask);           
%             mask=0.5*ones(size(vumeter));
%             set(hgui.speed_imgh,'Cdata',vumeter.*mask);
%             p.bDetect=0;
%             p.bShift = 0;
%             hgui.bRmsRepeat=1;  %1 
%             hgui.bSpeedRepeat=1;        %SC Make the speed monitor visible %1
%             
%             if (~isempty(rmsPeaks))
%                 p.rmsMeanPeak=mean(rmsPeaks);
%                 p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
%             end 
% 
%             hgui.showTextCue=1;
%         case 'start'
%             % SC(2008/06/10) Manually determine the optimum tracking params
%             % Warning: for consistency, don't change nDelay			
% 			set(hgui.msgh,'visible','on');
% %             set(hgui.msgh_imgh,'CData',CDataMessage.ftparampicking,'visible','on');
% 			drawnow;
% 			
% 			set(hgui.msgh,'string',{'Please stand by...'},'visible','on');
%             
%             if (hgui.debug==0)
%                 [vowelF0Mean,vowelF0SD]=getVowelPitches(dirname);
%                 disp(['Vowel meanF0 = ',num2str(vowelF0Mean),' Hz: stdF0 = ',num2str(vowelF0SD),' Hz']);
%                 disp(['Recommended cepsWinWidth = ',num2str(round(p.sr/vowelF0Mean*0.54))]);
%                 [vowelF1Mean,vowelF2Mean]=getVowelMeanF1F2(dirname);
%                 disp(['Vowel meanF1 = ',num2str(vowelF1Mean),' Hz; meanF2 = ',num2str(vowelF2Mean),' Hz']);
%                 optimFTParams=compareFormTrackParams(dirname);
%                 p.nLPC=optimFTParams.nLPC;
%                 p.nDelay=optimFTParams.nDelay;
%                 p.bufLen=(2*p.nDelay-1)*(p.frameLen);
%                 p.anaLen=p.frameShift+2*(p.nDelay-1)*p.frameLen;
%                 p.avgLen=optimFTParams.avgLen;
%                 p.bCepsLift=optimFTParams.bCepsLift;
%                 p.cepsWinWidth=optimFTParams.cepsWinWidth;
%                 p.fn1=optimFTParams.fn1;
%                 p.fn2=optimFTParams.fn2;       
%                 p.aFact=optimFTParams.aFact;
%                 p.bFact=optimFTParams.bFact;
%                 p.gFact=optimFTParams.gFact;
%             end
%             % ~SC(2008/06/10) Manually determine the optimum tracking
%                    
%             set(hgui.msgh,'string',{''},'visible','on'); 
%             p.rmsMeanPeak=mean(rmsPeaks);
%             p.rmsThresh=p.rmsMeanPeak/4;    %SC !! Adaptive RMS threshold setting. Always updating
% % 			if (p.rmsThresh>0.015)
% % 				p.rmsThresh=0.015;
% % 				disp('********* Warning: rms too high! Limited at 0.015. *********');
% % 			end
%             p.bDetect=0;
%             p.bShift=0;
%             set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             hgui.showTextCue=1;
%         case 'ramp'             %SC !! Notice that adaptive RMS threshold updating is no longer done here.           			
%             p.bDetect = 1;
% 			p.bShift = 1;
%             
% 			hgui.showTextCue=1;
% %             if doPlot
% %                 uiwait(gcf,10);
% % 			end
%         case 'stay'         	
% 			set(hgui.msgh,'visible','on');
%             p.bDetect = 1;
%             p.bShift = 1;
%             hgui.showTextCue=1;
%         case 'end'       			
% 			set(hgui.msgh,'visible','on');         
%             p.bDetect = 0;
%             p.bShift = 0;
%             hgui.showTextCue=1;
%     end
%     drawnow    

    set(hgui.strh, 'string', getMsgStr(n), 'visible', 'on'); 
    set(hgui.strh, 'FontSize', 30);

    MexIO('init', p);  %SC Inject p to TransShiftMex

%     for i0 = startRep:nReps    %SC Loop for the reps in the phase
%         repString=['rep',num2str(i0)];
%         state.rep=i0;
%         state.params=p;
%         save(fullfile(dirname,'state.mat'),'state');
%         
%         nTrials=length(expt.script.(thisphase).(repString).trialOrder);
% 
%         subsubdirname=fullfile(subdirname,repString);
%         mkdir(subsubdirname);
% 		
% 		% --- Perturbation field ---
% 		p.pertF2=linspace(p.F2Min,p.F2Max,p.pertFieldN);
%         switch (thisphase)
%             case 'ramp'
% 				p.pertAmp=i0/(expt.script.ramp.nReps+1)*subject.shiftRatio*ones(1,p.pertFieldN);
%             case 'stay'
%                 p.pertAmp=subject.shiftRatio*ones(1,p.pertFieldN);				
%             otherwise,
% 				p.pertAmp=zeros(1,p.pertFieldN);	
% 		end
% 		if isequal(subject.shiftDirection,'F1Up')
% 			p.pertPhi=0*ones(1,p.pertFieldN);
% 		elseif isequal(subject.shiftDirection,'F1Down')
% 			p.pertPhi=pi*ones(1,p.pertFieldN);
% 		end
% 		MexIO('init',p);
% 		% --- ~Perturbation field ---

    for k = 1:nTrials
        thisTrial = expt.expt_design.(thisrun).trialTypes(k); % 0: silent; 1: no noise; 2: noise only; 			
        thisWord = expt.expt_design.(thisrun).stimUtters{k};     %SC Retrieve the word from the randomly shuffled list

        hgui.trialType = thisTrial;
        hgui.word = thisWord;

        if (hgui.trialType == 4)	% Speech with masking noise or passively listening to masking noise
            TransShiftMex(3, 'datapb', gainMTB_fb * x_mtb, 0);
        end

        disp('');
        if (ischar(thisWord))
            disp([thisrun,' - ', 'trial #',num2str(k),': trialType = ',num2str(hgui.trialType), ...
                  ' - (', TRIAL_TYPES{hgui.trialType}, ') -  "', thisWord, '"']);
        else
            disp([thisphase,' - ', 'trial #', num2str(k), ': trialType = ',num2str(hgui.trialType),' - Pseudoword -']);
        end

        % Count down    
%         if ~(isequal(thisphase,'start') || isequal(thisphase,'ramp') || isequal(thisphase,'stay') || isequal(thisphase,'end'))
        disp(['Left: ',num2str(expt.expt_design.(thisrun).nTrials - k + 1), '/' , num2str(expt.expt_design.(thisrun).nTrials)]);
%         else 
%             if ~(isequal(thisphase,'ramp') || isequal(thisphase,'stay'))
%                 disp(['Left: ',num2str(expt.script.(thisphase).nTrials-phaseTrialCnt+1),'/',num2str(expt.script.(thisphase).nTrials),...
%                     ', ',num2str((expt.script.(thisphase).nTrials-phaseTrialCnt+1)*hgui.ITI),' sec']);
%             else
%                 disp(['Left: ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials-phaseTrialCnt+1),'/',...
%                     num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials),...
%                     ', ',num2str((expt.script.ramp.nTrials+expt.script.stay.nTrials-phaseTrialCnt+1)*hgui.ITI),' sec']);
%             end
%         end
        % ~Count down

        TransShiftMex(3, 'bpitchshift', 1, 1);
        if (hgui.trialType == 2 || hgui.trialType == 3)   %SC The distinction between train and test words
            if hgui.trialType == 2 % Down shift
                TransShiftMex(3, 'pitchshiftratio', 2 ^ (1 / 12 * expt.expt_config.PITCH_SHIFT_DOWN), 1);
            else
                TransShiftMex(3, 'pitchshiftratio', 2 ^ (1 / 12 * expt.expt_config.PITCH_SHIFT_UP), 1);
            end
        else
            TransShiftMex(3, 'pitchshiftratio', 1, 1);
        end

%         if (thisTrial==5)
%             hgui.skin.facePnt=expt.script.(thisphase).(repString).face(k);
%         end

        UIRecorder('singleTrial', hgui.play, 1, hgui);
        data=get(hgui.UIrecorder,'UserData');           %SC Retrieve the data
        data.timeStamp = clock;
        data.expt = expt;
        data.params.name = thisWord;
        data.params.trialType = thisTrial;

%         if (thisTrial == 1)
%             if ~isempty(data.rms)
%                 switch (thisphase)  %SC Record the RMS peaks in the bout
%                     case 'pre'
%                         rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];
%                     case 'pract1',
%                         rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];                    
%                     case 'pract2',
%                         rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];                    
%                     otherwise,
%                 end
%             end
%         end

%         if (isequal(thisphase,'pract1'))
%             if (thisTrial==1 || thisTrial==2)
%                 if (isfield(data,'vowelLevel') && ~isempty(data.vowelLevel) && ~isnan(data.vowelLevel) && ~isinf(data.vowelLevel))
%                     subjProdLevel=[subjProdLevel,data.vowelLevel];
%                 end
%             end
%         end
    
        dataFN = fullfile(subdirname, ['trial-', num2str(k), '-', num2str(thisTrial)]);
        save(dataFN,'data');
        disp(['Saved data to ', dataFN, '.']);
        disp(' ');

%         phaseTrialCnt=phaseTrialCnt+1;

        % Calculate and show progress
%         if (subject.showProgress)
%             rProgress=calcExpProgress(expt,thisphase,i0,k,rProgress);
%             if (~isnan(rProgress))
%                 progress_mask=zeros(size(progress_meter));
%                 progress_mask(:,1:round(rProgress*100),:)=1;            
%                 set(hgui.progress_imgh,'Cdata',progress_meter.*progress_mask);
%             end
%         end
            
%         end
    end
    
    n = n + 1;
end
set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
set(hgui.msgh,'string',...
	{'Congratulations!';...
	'You have finished the expt.'},'visible','on');
% set(hgui.msgh_imgh,'CData',CDataMessage.finish,'visible','on');
pause(3);
close(hgui.UIrecorder)
pause(2);
% saveExperiment(dirname);

save(fullfile(dirname,'expt.mat'),'expt');
save(fullfile(dirname,'state.mat'),'state');

return

%% 
function phaseScript=genPhaseScript(stage,nReps,trialTypes,trainWords,testWords,pseudoWords,randReps,designNum);
	phaseScript=struct();
	phaseScript.nReps=nReps;
    phaseScript.nTrials=0;
    
	if (designNum==1)

		if mod(nReps,randReps)~=0
			disp('Warning: nReps not divided by randReps. Returning empty phaseScript');
			return
		end

		if (nReps<length(trainWords))
			trainWordsUsed=trainWords(1:nReps);
			testWordsUsed=testWords(1:nReps);
			pseudoWordsUsed1=pseudoWords(1:nReps);
			pseudoWordsUsed2=pseudoWords(1:nReps);
			trainWordsUsed=trainWordsUsed(randperm(nReps));
			testWordsUsed=testWordsUsed(randperm(nReps));
			pseudoWordsUsed1=pseudoWordsUsed1(randperm(nReps));
			pseudoWordsUsed2=pseudoWordsUsed2(randperm(nReps));
		else
			trainWordsUsed=cell(1,0);
			testWordsUsed=cell(1,0);
			pseudoWordsUsed1=[];
			pseudoWordsUsed2=[];
			for n=1:(nReps/length(trainWords))
				trainwu1=trainWords(randperm(length(trainWords)));
				testwu1=testWords(randperm(length(testWords)));
				pseudowu1=pseudoWords(randperm(length(pseudoWords)));
				pseudowu2=pseudoWords(randperm(length(pseudoWords)));
				trainWordsUsed=[trainWordsUsed,trainwu1];
				testWordsUsed=[testWordsUsed,testwu1];
				pseudoWordsUsed1=[pseudoWordsUsed1,pseudowu1];
				pseudoWordsUsed2=[pseudoWordsUsed2,pseudowu2];
			end
		end
		
		if isequal(stage,'pract1') | isequal(stage,'pract2')
			trialTypes=trialTypes(find(trialTypes <= 2));
		end

		repCnt=1;
		for n=1:nReps/randReps
			unitRep=struct;
			tmp_trialType=[];
			for k=1:randReps
				tmp_trialType=[tmp_trialType,trialTypes];
			end		
			unitRep.trialOrder=tmp_trialType(randperm(length(tmp_trialType)));
			unitRep.word=cell(size(unitRep.trialOrder));
			for m=1:length(unitRep.trialOrder)
				if (unitRep.trialOrder(m)==1)		% Speech with auditory feedback
					unitRep.word{m}=trainWordsUsed{repCnt};
				elseif (unitRep.trialOrder(m)==2)	% Speech without auditory feedback
					unitRep.word{m}=testWordsUsed{repCnt};
				elseif (unitRep.trialOrder(m)==3)
					unitRep.word{m}=pseudoWordsUsed1(repCnt);
				elseif (unitRep.trialOrder(m)==4)
					unitRep.word{m}=pseudoWordsUsed2(repCnt);
				end
			end
			for k=1:randReps
				oneRep=struct;
				oneRep.trialOrder=unitRep.trialOrder((k-1)*length(trialTypes)+1:k*length(trialTypes));
				oneRep.word=unitRep.word((k-1)*length(trialTypes)+1:k*length(trialTypes));
				phaseScript.(['rep',num2str(repCnt)])=oneRep;
                phaseScript.nTrials=phaseScript.nTrials+length(oneRep.trialOrder);
				repCnt=repCnt+1;

            end
            
            if (n==nReps/randReps)
                phaseScript.(['rep',num2str(phaseScript.nReps)]).trialOrder=[phaseScript.(['rep',num2str(phaseScript.nReps)]).trialOrder,4];
                phaseScript.(['rep',num2str(phaseScript.nReps)]).word{length(phaseScript.(['rep',num2str(phaseScript.nReps)]).word)+1}=1;
                phaseScript.nTrials=phaseScript.nTrials+1;
            end
		end
	elseif (designNum==2)
		for n=1:nReps
			bt=[zeros(1,length(trainWords)),ones(1,round(length(trainWords)/2))];
            trainWordsUsed=trainWords(randperm(length(trainWords)));
			pseudoWordsUsed=pseudoWords(randperm(length(trainWords)));
%             testWordsUsed2=testWords(randperm(length(testWords)));            
			twCnt=1;
			bt=bt(randperm(length(bt)));
			oneRep=struct;
			oneRep.trialOrder=[];
			oneRep.word=cell(1,0);
            cntTW=1;
			for m=1:length(bt)
				if (bt(m)==0)					
					oneRep.trialOrder=[oneRep.trialOrder,1];
					oneRep.word{length(oneRep.word)+1}=trainWordsUsed{twCnt};
					twCnt=twCnt+1;
				elseif (bt(m)==1)
					oneRep.trialOrder=[oneRep.trialOrder,[5,4,4]];
					oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
					oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW);
					oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
                    cntTW=cntTW+2;
                end
            end
            
            if (isequal(stage,'pract1') || isequal(stage,'pract2'))
                idx=find(oneRep.trialOrder<4);
                oneRep.trialOrder=oneRep.trialOrder(idx);
                oneRep.word=oneRep.word(idx);
            end
            
            if n==nReps
                oneRep.trialOrder=[oneRep.trialOrder,4];    % Dummy trial at the end
                oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(1);
            end
			phaseScript.(['rep',num2str(n)])=oneRep;
            phaseScript.nTrials=phaseScript.nTrials+length(oneRep.trialOrder);
		end
	end
return

function script1=addFaceInfo(script0,dFaces)
    script1=script0;
    
    IDs.male=[];
    IDs.female=[];
    
    for n=1:length(dFaces.d)
        if isequal(dFaces.sex{n},'M')
            if isempty(find(IDs.male==dFaces.subjID(n)))
                IDs.male=[IDs.male,dFaces.subjID(n)];
            end
        elseif isequal(dFaces.sex{n},'F')
            if isempty(find(IDs.female==dFaces.subjID(n)))
                IDs.female=[IDs.female,dFaces.subjID(n)];
            end
        end
    end    
    IDs.male=sort(IDs.male);
    IDs.female=sort(IDs.female);
    
    IDs.male=IDs.male(randperm(length(IDs.male)));
    IDs.female=IDs.female(randperm(length(IDs.female)));    
    
    stages=fields(script1);
    for n=1:length(stages);
        stg=stages{n};
        if ~isempty(find(script1.(stg).rep1.trialOrder==5))
            nPersons=script1.(stg).nReps;
            tPersons=[];
            for k=1:nPersons
                if (mod(k,2)==0)    % female
                    tPersons=[tPersons,IDs.female(1)];
                    IDs.female=IDs.female(2:end);
                else    % male
                    tPersons=[tPersons,IDs.male(1)];
                    IDs.male=IDs.male(2:end);
                end
            end
            idxFaces=[];
            for k=1:length(tPersons)
                idxFaces0=find(dFaces.subjID==tPersons(k));
                idxFaces0=idxFaces0(randperm(length(idxFaces0)));
                idxFaces=[idxFaces,idxFaces0(1:4)];
            end
            idxFaces=idxFaces(randperm(length(idxFaces)));
            cnt=1;
            for k=1:script1.(stg).nReps
                repString=['rep',num2str(k)];
                script1.(stg).(repString).face=zeros(size(script1.(stg).(repString).trialOrder));
                idx=find(script1.(stg).(repString).trialOrder==5)
                script1.(stg).(repString).face(idx)=idxFaces(cnt:cnt+length(idx)-1);
                cnt=cnt+length(idx);
            end
        end
    end
return