%% Transition sensorimotor adaptation
% clear all;
% close all hidden;
% clc;

%% Subject information and settings
% disp('--- Adjust Microphone gain knob to 0 (center) ! ... ---');

subject.nr              = 104;   %SC-Mod(12/04/2007)
subject.name            = 'PS004';
subject.sex             = 'male';  % male /z female
% subject.pertMode        = 'parallel'; %SC parallel / cosine
% subject.shiftDirection  = 'AccelDecel';  %SC inflate / deflate
% subject.shiftRatio      = 0.25;
% subject.date            = clock;

% subject.lvNoise2        = 76;  % dBA. The level of noise for completely masking speech (mode fb = 2).
% subject.lvNoise3        = 0;   % dBA. The level of noise for masking bone conduction (mode fb = 3). 
% 
% subject.microphoneGain  = 0;    % ASSUMED
% subject.pcrKnob         = input('Phone / Ctrl Room knob (-2.5 ~ 0.5) = ');
% subject.m2mDist         = input('Mouth-to-microphone distance (cm): '); % Unit: cm
% subject.closedLoopGain  = calcClosedLoopGain(subject.pcrKnob);
% 
% subject.showProgress    = 1;

%% Other experiment-related parameters
% twoScreens=1;
% break1=60;

% msgBreak={'Please take a break for 2 minutes','',...
%     'Feel free to take off the ear phones and walk around.',...
%     '',...
%     'When you are about to continue, make sure that you ',...
%     'inserted the ear tips as deeply as possible.',...
%     '',...
%     'When you are ready to continue the experiment, press',...
%     'OK to continue.'};

%% Production data directory
% prodDataDir             = ['G:\CS_2004\PROJECTS\TransAdapt\TRIPHF2\S',num2str(subject.nr)];

%% Welcome information
% mymsgbox({'Welcome to the experiment!',...
%     '',...
%     'The experiment consists of three parts, ',...
%     '    1) a hearing test,',...
%     '    2) a production (speaking) test, ',...
% 	'    3) a perceptual discriminability test.'...
%     'The first part will take about 10 minutes; the second part about 50',...
%     'minutes; the third part another 50 minutes ',...
%     '',...
%     'We will now measure your auditory acuity ...'},...
%     'modal','twoScreens',twoScreens,'pause',5);

%% Audiogram measurement #1
wordInfo_dummy.prodDir='';
wordInfo_dummy.stage='';
wordInfo_dummy.rep='';
wordInfo_dummy.word='';

disp('--- Adjust Phones/CtrlRoom knob to -1.5! ... ---');
pause;

cd ../../adapt/perctest
ai_discrim(['SAP-FMRI_',subject.name,'_audiogram1'],'ai',4,pwd,4,wordInfo_dummy,NaN,'twoScreens');
close all hidden;
cd ../../sap-fmri/mcode
