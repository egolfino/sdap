function varargout = UIRecorder(varargin)
% UIRECORDER M-file for UIrecorder.fig
%      UIRECORDER, by itself, creates a new UIRECORDER or raises the existing
%      singleton*.
%
%      H = UIRECORDER returns the handle to a new UIRECORDER or the handle to
%      the existing singleton*.
%
%      UIRECORDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UIRECORDER.M with the given input arguments.
%
%      UIRECORDER('Property','Value',...) creates a new UIRECORDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UIrecorder_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UIrecorder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UIrecorder

% Last Modified by GUIDE v2.5 29-Jul-2012 16:29:08

%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UIrecorder_OpeningFcn, ...
    'gui_OutputFcn',  @UIrecorder_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% SCai: data displaying function
if (~isempty(findStringInCell(varargin,'figIdDat'))) 
    figIdDat=varargin{findStringInCell(varargin,'figIdDat')+1};
end

%% --- Executes just before UIrecorder is made visible.
function UIrecorder_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
if (length(varargin)>=2)
    figIdDat=varargin{2};
else
    figIdDat=[];
end

if (~isempty(findStringInCell(varargin,'trad')))
	isTrad=1;
else 
	isTrad=0;
end

handles.trigByScanner=0;
% handles.TR = 7.5;
% handles.TA = 2.5;
handles.run = '';
handles.trigKey='add';	% To change
handles.debug=0;
handles.vumeterMode=NaN;  % 1: 10 ticks; 2: 3 ticks;

handles.timeCreated=clock;

% handles.msgImgDir='../../adapt/triphcode/uimg';
% handles.utterImgDir='../../adapt/triphcode/utterimg';
if (isTrad)
	handles.msgImgDir=fullfile(handles.msgImgDir,'trad');
	handles.utterImgDir=fullfile(handles.utterImgDir,'trad');
end

set(hObject,'visible','off');
set(handles.UIrecorder,'interruptible','on','busyaction','queue')
%--------------------------------------------------------------------------
%SC Construct the volume/speed indicator template
vumeter  = permute(jet(100),[1,3,2]);
color1=vumeter(20,:,:);
color2=vumeter(53,:,:);
color3=vumeter(90,:,:);
for n=1:100
    if n<=30
        vumeter(n,:,:)=color1;
    elseif n<=70
        vumeter(n,:,:)=color2;
    else
        vumeter(n,:,:)=color3;
    end
end
vumeter2=vumeter;
vumeter(10:10:90,:)=0;    %SC-Commented(12/11/2007)
vumeter0=nan(size(vumeter,2),size(vumeter,1),size(vumeter,3));
vumeter0(:,:,1)=transpose(vumeter(:,:,1));
vumeter0(:,:,2)=transpose(vumeter(:,:,2));
vumeter0(:,:,3)=transpose(vumeter(:,:,3));
vumeter=vumeter0;
% vubounds=[1,30,70,100];%SC The boundaries are at 29 and 69.

vumeter2([30,70],:)=0;
vumeter02=nan(size(vumeter2,2),size(vumeter2,1),size(vumeter2,3));
vumeter02(:,:,1)=transpose(vumeter2(:,:,1));
vumeter02(:,:,2)=transpose(vumeter2(:,:,2));
vumeter02(:,:,3)=transpose(vumeter2(:,:,3));
vumeter2=vumeter02;
%--------------------------------------------------------------------------
%SC Construct the progress indicator template
progressmeter=1*ones(1,100,3);
%SC ~Construct the progress indicator template
%--------------------------------------------------------------------------
if (handles.vumeterMode==1)
    handles.rms_imgh = image(vumeter,'parent',handles.rms_axes);
    handles.speed_imgh=image(vumeter,'parent',handles.speed_axes);
    set(handles.rms_imgh,'CData',zeros(size(vumeter)));
    set(handles.speed_imgh,'CData',zeros(size(vumeter)));
else
    handles.rms_imgh = image(vumeter,'parent',handles.rms_axes);
    handles.speed_imgh=image(vumeter,'parent',handles.speed_axes);
    set(handles.rms_imgh,'CData',zeros(size(vumeter2)));
    set(handles.speed_imgh,'CData',zeros(size(vumeter2))); 
end

if (~isempty(findStringInCell(varargin,'showVuMeter')))    
    set(handles.rms_imgh,'CData',vumeter0);
end



% if (~isempty(findStringInCell(varargin,'showVuMeter')))
%     set(handles.speed_imgh,'CData',vumeter0);
% end

% set(handles.phrase_axes,'Box','off');
% set(handles.axes_msgh,'Box','off');
set(handles.axes_pic,'Box','off');


handles.progress_imgh = image(progressmeter,'parent',handles.progress_axes);

set(handles.rms_label,'string','Volume');
set(handles.speed_label,'string','Speed');

handles.pcrKnob=NaN;
% handles.trialType=4;
% handles.word='Ready...';

% handles.bAuto=1;

handles.time1=[];
handles.time2=[];

% set(handles.auto_btn,'Value',get(handles.auto_btn,'Max'));
skin=struct('pause', imread(fullfile(pwd,'graphics','skin-pause.jpg')),...
    'play', imread(fullfile(pwd,'graphics','skin-play.jpg')),...
    'good', imread(fullfile(pwd,'graphics','choice-yes.gif')),...
    'bad', imread(fullfile(pwd,'graphics','choice-cancel.gif')),...
	'fixation',imread(fullfile(pwd,'graphics','fixation.bmp')),...
	'vumeter',  vumeter,...
    'vumeter2',  vumeter2,...
	'dFaces',getDFaces('./graphics/faces/face*.bmp'),...
	'dPseudowords',dir('./graphics/pseudochars/pseudoword-*.bmp'),...
    'dWords',dir('./graphics/words/word*.bmp'),...
	'faceOrder',[],...
	'facePnt',1);
handles.skin=skin;

handles.pic_imgh=image(handles.skin.fixation,'parent',handles.axes_pic);
set(handles.pic_imgh,'visible','off');

% set(handles.prev,'cdata',skin.prev);
% set(handles.next,'cdata',skin.next);
set(handles.play,'cdata',skin.play);
set(handles.play,'UserData',0);
set(handles.play,'Value',get(handles.play,'Min'));
set(handles.rms_axes,'xtick',[],'ytick',[]);
axis(handles.rms_axes, 'xy');

set(handles.speed_axes,'xtick',[],'ytick',[])
axis(handles.speed_axes, 'xy');

set(handles.axes_pic,'xtick',[],'ytick',[],'box','off','visible','off');

set(handles.progress_axes,'xtick',[],'ytick',[]);
axis(handles.progress_axes, 'xy');

handles.figIdDat=figIdDat;

handles.dataOut=[];
handles.bRmsRepeat=0;
handles.bSpeedRepeat=0;
handles.vocaLen=NaN;    %SC-Mod(2008/01/05) Old value: 300 
handles.lenRange=NaN;   %SC(2008/01/05)

% handles.ITI = 6;			%SC(2009/02/05) Inter-trial interval
handles.TR = 7.5;
handles.TA = 2.5;
handles.stimDelay = 0.5;
handles.recDelay = 0.5;

handles.trialLen = 2.5;
 
handles.showTextCue=0;  %SC(2008/01/06)

handles.dBRange=NaN;
handles.rmsTransTarg_spl=NaN;
% load calibMic;  % gets micGain: wav rms/ Pa rms (Pa^-1)
load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
handles.rmsTransTarg=micRMS_100dBA / (10^((100-handles.rmsTransTarg_spl)/20));

% handles.nextMessage=imread(fullfile(handles.msgImgDir,'message_pre2.bmp'));

handles.exptDir = '';

guidata(hObject, handles);

set(handles.UIrecorder,'keyPressFcn',@key_Callback);
set(handles.strh,'keyPressFcn',@key_Callback);
set(handles.play,'keyPressFcn',@key_Callback);
% set(handles.rec_slider,'keyPressFcn',@key_Callback);
set(handles.msgh,'keyPressFcn',@key_Callback);

set(handles.strh,'string','HELLO','visible','on');

set(hObject,'visible','on');

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes UIrecorder wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% --- Outputs from this function are returned to the command line.
function varargout = UIrecorder_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
varargout{1} = handles;

%% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play

% set(handles.button_next,'visible','off');

% CDataMessageBlank=zeros(750,720,3);
% CDataMessageBlank(:,:,1)=64/255*ones(750,720);
% CDataMessageBlank(:,:,3)=64/255*ones(750,720);

if(get(handles.play,'userdata')==0) % currently in pause mode
    set(handles.play,'cdata',handles.skin.pause,'userdata',1); % now in play mode
    set(handles.msgh,'string','');
	handles.trialType=-1;
	handles.word='Ready...';	

    singleTrial(handles.play,[],handles)   %%SC
else % currently in play mode
    set(handles.play,'cdata',handles.skin.play,'userdata',0); % now in pause mode
    TransShiftMex(2) %%SC stop TransShiftMex
    set(handles.msgh,'string','Press play to continue...');
end

function key_Callback(src, evnt)
hgui = guidata(src);
timeNow = clock;
eTime = etime(timeNow, hgui.timeCreated);

if (isequal(evnt.Key,hgui.trigKey) || isequal(evnt.Key,'a'))
    if ~isempty(hgui.exptDir) && isfile(fullfile(hgui.exptDir, 'lastTrigTime.mat'))
        load(fullfile(hgui.exptDir, 'lastTrigTime.mat')); % gives lastTrigTime
        fprintf('--> Trigger at %.2f s (%.2f s from last trigger) <--\n', ...
                eTime, eTime - lastTrigTime);
    else
        fprintf('--> Trigger at %.2f s <--\n', eTime);
    end
    
    uiresume(hgui.UIrecorder);
else
% 	set(hgui.UIrecorder,'UserData','nogo');
end

lastTrigTime = eTime;
if ~isempty(hgui.exptDir)
    save(fullfile(hgui.exptDir, 'lastTrigTime.mat'), 'lastTrigTime');
end

return

% %% --- Executes on button press in prev.
% function prev_Callback(hObject, eventdata, handles)
% % hObject    handle to prev (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% singleTrial(handles.prev,[],handles);

%% --- Executes on button press in next.
% function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% --- Single trial callback function.
function singleTrial(hObject, eventdata, handles, varargin)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.next,'enable','off','visible','off')
% set(handles.prev,'enable','off','visible','off')
record(handles);

%% SUBFUNCTIONS
%% --------------------------------------------------------------------------


%% --------------------------------------------------------------------------
function startRec(obj,event,handles)
% startRec displays a string and starts a timer object (handles.rect_timer)
% who's TimerFcn Callback (@stopRec) is called after a timeout period, set by
% the Value of the slider (handles.rec_slider)
% the StartFcn / StopFcn of the timer object starts/stops the recording
% (@soundDevice('start')/('stop')
% CDataPhrase=imread('utterimg/phrase.bmp');

fprintf('startRec\n')

handles.dataOut=[];
str=get(handles.strh,'string');
set(handles.rms_imgh,'Cdata',zeros(size(get(handles.rms_imgh,'Cdata'))));
set(handles.speed_imgh,'Cdata',zeros(size(get(handles.speed_imgh,'Cdata'))));

% set(handles.phrase_axes,'CData',CDataPhrase,'visible','on');

% clockStart=clock;
% clockStart(6) = clockStart(6) + get(handles.rec_slider,'Value');
% startat(handles.rec_timer, clockStart);
set(handles.strh,'string',str);

%% --------------------------------------------------------------------------
function stopRec(obj,event,handles)
% this function stops the recording if the timer object (handles.rec_timer)
% is still running. After the recording is stoped, checkData is executed,
% and the next and previous buttons are activated
fprintf('stopRec\n')
% if(strcmp(get(handles.rec_timer,'Running'),'on'))
%     stop(handles.rec_timer)
% end

handles.dataOut=checkData(getData,handles);


guidata(handles.UIrecorder, handles);

% set(handles.next,'enable','on')
% set(handles.prev,'enable','on')
% set(handles.next,'visible','on')
% set(handles.prev,'visible','on')
% if(get(handles.auto_btn,'Value')==get(handles.auto_btn,'Max'))
%     next_Callback(handles.UIrecorder,[],handles)
% end

%% --------------------------------------------------------------------------
function soundDevice(obj,event,handles,action)
% interface to teh external sounddevice
switch(action)
    case 'init'
        TransShiftMex(0)
    case 'start'
        TransShiftMex(1)
    case 'stop'
        TransShiftMex(2)
    otherwise,
end

%% --------------------------------------------------------------------------
function dataOut= getData
% gets the data
dataOut=MexIO('getData');

%% --- Executes on button press in auto_btn.
% function auto_btn_Callback(hObject, eventdata, handles)
% hObject    handle to auto_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_btn

%%
function record(handles)

%     tWords={'BECK','BET','DECK','DEBT','PECK','PEP','PET','TECH'};

    if (handles.trigByScanner==1)
        set(handles.play,'userdata',1);
        uiwait(handles.UIrecorder);
    else
        waitfor(handles.play,'userdata',1);
    end

if (handles.debug==0)
    
    go=get(handles.UIrecorder,'UserData');
    if isequal(go,'nogo')
        return
    end

    handles.dataOut=[];
    guidata(handles.UIrecorder,handles);

    set(handles.strh,'visible','off');
    set(handles.msgh,'visible','off');

    handles.time1=clock;

%     if (handles.trigByScanner==1)
%         if (~isequal(handles.run,'pract1') && ~isequal(handles.run,'pract2'))
    pause(handles.TA);
%         else
%             pause(0.25);
%         end
%     else
%         pause(0.25);
%     end

    set(handles.pic_imgh,'cdata',handles.skin.fixation,'visible','on');
    drawnow;
    pause(handles.stimDelay);
    set(handles.pic_imgh,'visible','off');

    if (isequal(handles.word,'Ready...') || handles.trialType== -1)
        return
    end
    % if (handles.trialType==5)	% 5

    % 	if (rand>0.5)
    % 		set(handles.pic_imgh,'cdata',handles.skin.blueball);
    % 		dataOut.ballColor='blue';
    % 		dataOut.targetKey='leftarrow';
    % 	else
    % 		set(handles.pic_imgh,'cdata',handles.skin.redball);
    % 		dataOut.ballColor='red';
    % 		dataOut.targetKey='rightarrow';
    % 	end
    % 	waitTime=get(handles.rec_slider,'Max')-get(handles.rec_slider,'Value')+2.0+rand*0.25	
    % 	tic;
    % 	uiwait(handles.UIrecorder,waitTime);
    % 	a=toc;
    % 	tKey=get(handles.UIrecorder,'UserData');
    % 	dataOut.respTime=a;
    % 	if isstruct(tKey)
    % 		dataOut.results='timeout';
    % 		if isequal(dataOut.ballColor,'red')
    % 			set(handles.pic_imgh,'cdata',handles.skin.redball_wrong);
    % 		else 
    % 			set(handles.pic_imgh,'cdata',handles.skin.blueball_wrong);
    % 		end
    % 	else
    % 		dataOut.timeOut=0;
    % 		if isequal(tKey,dataOut.targetKey)
    % 			dataOut.result='correct';
    % 			if isequal(dataOut.ballColor,'red')
    % 				set(handles.pic_imgh,'cdata',handles.skin.redball_correct);
    % 			else 
    % 				set(handles.pic_imgh,'cdata',handles.skin.blueball_correct);
    % 			end
    % 		else
    % 			dataOut.result='wrong';
    % 			if isequal(dataOut.ballColor,'red')
    % 				set(handles.pic_imgh,'cdata',handles.skin.redball_wrong);
    % 			else 
    % 				set(handles.pic_imgh,'cdata',handles.skin.blueball_wrong);
    % 			end
    % 		end	
    % 	end
    % 	set(handles.UIrecorder,'UserData',dataOut);
    % 	
    % else

        % Put it here to balance timing !!!

    % 	nums=round(rand(1,3)*(length(handles.skin.dPseudochars)-1))+1;
    % 	imChar{1}=imread(fullfile(pwd,'graphics','pseudochars',handles.skin.dPseudochars(nums(1)).name));
    % 	imChar{2}=imread(fullfile(pwd,'graphics','pseudochars',handles.skin.dPseudochars(nums(2)).name));
    % 	imChar{3}=imread(fullfile(pwd,'graphics','pseudochars',handles.skin.dPseudochars(nums(3)).name));
    % 	imBG(51:130,21:80,:)=imChar{1};
    % 	imBG(51:130,91:150,:)=imChar{2};
    % 	imBG(51:130,161:220,:)=imChar{3};

%         if (handles.trialType >=1 &&  handles.trialType <= 4) % Speech trials
        imBG=imread(fullfile(pwd,'graphics','pseudochars',handles.skin.dPseudowords(ceil(rand*8)).name));	% Balance timing

%             set(handles.pic_imgh, 'Visible', 'off');
        set(handles.strh, 'FontSize', 45);
        set(handles.strh, 'String', handles.word, 'Visible', 'on');
        drawnow;
%             k=findStringInCell(tWords,handles.word);
%             imWord=imread(fullfile(pwd,'graphics','words',handles.skin.dWords(k).name));
%             nFace=round(length(handles.skin.dFaces.d)*rand);    % Balance timing
%             if (nFace==0) nFace=1;  end
%             imFace=imread(['./graphics/faces/',handles.skin.dFaces.d(nFace).name]);
    % 		set(handles.strh,'string',handles.word,'visible','on');
%             set(handles.pic_imgh,'cdata',imWord,'visible','on');
%         else % Baseline trial
%     % 		set(handles.strh,'string',handles.word,'visible','off');
% 
%             imBG=imread(fullfile(pwd,'graphics','pseudochars',handles.skin.dPseudowords(handles.word).name));
% %             imWord=imread(fullfile(pwd,'graphics','words',handles.skin.dWords(ceil(rand*8)).name));
% 
%             if (handles.trialType==5)
%                 nFace=handles.skin.facePnt;
% %                 dt=dir(['./graphics/faces/face',num2str(nFace),'*.bmp']);
%                 dfn=handles.skin.dFaces.d(nFace).name;
%                 imFace=imread(['./graphics/faces/',dfn]);
% 
%                 y=round((size(imBG,1)-size(imFace,1)-1)*0.5);
% %                 y=round((size(imBG,1)-size(imFace,1)-1)*rand);
% %                 if (y==0)
% %                     y=1;
% %                 end
%                 x=round((size(imBG,2)-size(imFace,2)-1)*rand);
%                 if (x==0)
%                     x=1;
%                 end                
%                 imBG([y:y+size(imFace,1)-1],[x:x+size(imFace,2)-1],:)=imFace;
%             else
% %                 nFace=round(length(handles.skin.dFaces.d)*rand);    % Balance timing
% %                 if (nFace==0) nFace=1;  end
% %                 imFace=imread(['./graphics/faces/',handles.skin.dFaces.d(nFace).name]);
%             end
% 
%     % 		set(handles.strh,'string',handles.word,'visible','off');
%             set(handles.pic_imgh,'cdata',imBG,'visible','on');
%         end

        if (handles.trialType == 4)
            TransShiftMex(3, 'fb', 2); % Noise masking
        else
            TransShiftMex(3, 'fb', 1); % Normal or perturbed auditory feedback, no added noise        
        end
        
        TransShiftMex(1);
        pause(handles.trialLen);  % Changed 2008/06/18 to make the pause longer +1.5 --> +2.0
%         if get(handles.play,'userdata') == 0 % in pause mode
%             record(handles); % re-do recording
%         end
        TransShiftMex(2);

        [dataOut, bRmsGood, bSpeedGood]= checkData(getData, handles);

        bRmsRepeat = handles.bRmsRepeat;
        bSpeedRepeat = handles.bSpeedRepeat;
%         if (handles.trialType==3 | handles.trialType==4)
%             if (bRmsRepeat==1)
%                 bRmsRepeat=0;
%             end
%             if (bSpeedRepeat==1)
%                 bSpeedRepeat=0;
%             end
%         end

        if ((~bRmsGood && bRmsRepeat) || (~bSpeedGood && bSpeedRepeat))
            record(handles)    
        else
        % data is saved as UserData in the fig handle (wicht is the signal for
        % the host function to launch the next single trial
        set(handles.UIrecorder, 'UserData', dataOut)
        end
    % end
    if (handles.trigByScanner==0)
        pause(0.25);
    end
    set(handles.strh,'visible','off');
    % set(handles.pic_imgh,'cdata',handles.skin.fixation);
    set(handles.pic_imgh, 'visible', 'on');
    drawnow;
else
    dataOut=struct;
    dataOut.signalIn=[]; dataOut.signalOut=[];
    dataOut.rms=[];
    set(handles.UIrecorder,'UserData',dataOut);
end
return

% --- Executes on button press in button_next.
% function button_next_Callback(hObject, eventdata, handles)
% % hObject    handle to button_next (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if (isfield(handles,'nextMessage'))
%     set(handles.msgh_imgh,'CData',handles.nextMessage,'visible','on');
% end
% set(handles.button_next,'visible','off');
% 
% set(handles.play,'visible','on');

function dFaces=getDFaces(fileMask)
    dFaces=struct;
    dFaces.d=dir(fileMask);
    dFaces.sex=cell(1,length(dFaces.d));
    dFaces.subjID=nan(1,length(dFaces.d));
    
    for n=1:length(dFaces.d)
        fn=dFaces.d(n).name;
        idx=strfind(fn,'.bmp');
        dFaces.sex{n}=dFaces.d(n).name(idx-1);
        idx1=strfind(fn,'-s');
        idx2=strfind(fn,['-',dFaces.sex{n}]);
        dFaces.subjID(n)=str2num(dFaces.d(n).name(idx1+2:idx2-1));
    end
return
