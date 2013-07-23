clear all;
addpath 'C:\Users\egolfino\Documents\MATLAB\Cogent2000v1.32\Toolbox'

% Subject Directory
dataDIR = 'C:\SDAP';

%Setup for display
full_screen = 1;                 % window mode ( 0=partial window, 1=full screen )
keyboard_timing = 'exclusive'; % 'nonexclusive'= timing is poor; 'exclusive' timing is accurate
screen_resolution = 3;           % 1024 x 768
background_color = [0, 0, 0];    % black 
font_color = [1, 1, 1];          % white
font_size = 100;				     % font size
font_type = 'Arial'; 			 % font type
number_buffers = 5;              % number of offscreen buffers
number_bitspixel = 0;            % number of bits/pixel 0 is maximum

%Setup for trial 
trialDur = 14444;   			 % length of a trial in ms.
totrunDur=260000;                % length of a run in ms.
totalxDur = 11940;               % total duration xhair
number_runs = 3;                 % total possible number of runs
repetitionTime = 2500;           % Image acquisition time      
begnumdummies = 3;               % Number of dummy scans 

% Setup for data storage
day=date;
time = fix(clock);
hour = num2str(time(4));
min = num2str(time(5));

%Gets user's input on subject number and run converts these into
%usable variables 
subNumStr = input('Please enter a subject number: ', 's');
subNum = str2double(subNumStr);

%These numbers are very important ensure that the investigator knows what
%she's entering!

while isempty(subNum) 
    warning('That is not a number');
    subNumStr = input('Please re-enter a subject number: ', 's');
    subNum = str2double(subNumStr);
end

runNumStr = input('Please enter a run number: ', 's');
runNum = str2num(runNumStr);

while isempty(runNum) 
    warning('That is not a number');
    runNumStr = input('Please re-enter a run number: ', 's');
    runNum = str2double(runNumStr);
end

while runNum > number_runs
    warning('File not found');
    runNumStr = input('Please renter a run number: ', 's');
    runNum = str2double(runNumStr);
end

datafile=strcat('run',runNumStr,'.dat');
     
if subNum < 10
    eval(strcat('mkdir C:\SDAP_DATA\SDAP_PIL0',subNumStr,'_LOC080713\'));   
    logfile=strcat('C:\SDAP_DATA\SDAP_PIL0',subNumStr,'_LOC080713\RUN0',runNumStr,'_',day,'-',hour,'-',min,'.log');
else
    eval(strcat('mkdir C:\SDAP_DATA\SDAP_PIL',subNumStr,'_LOC080713\')); 
    logfile=strcat('C:\SDAP_DATA\SDAP_PIL',subNumStr,'_LOC080713\RUN0',runNumStr,'_',day,'-',hour,'-',min,'.log');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           END OF SETUP                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cgloadlib; % loads Cogent Graphics commands - without which cg* commands don't run
config_display( full_screen, screen_resolution, background_color, font_color, font_type, font_size, number_buffers, number_bitspixel);
config_keyboard(20,5,keyboard_timing);
config_data( datafile );
config_log( logfile );
start_cogent; %matlab session is given a higher priority than all these other tasks

numstim = countdatarows;
totxdur(1)=0;
for i = 1:numstim
    stim{i} = getdata( i, 1 ); %gets stimulus string
    stimdur(i) = getdata(i,2); %loads length of time stimulus is presented
    totxdur(i+1) = getdata(i,3); %loads length of time + is presented
end;

% Log information about experiment and stimuli
logstring (['Subject: ' subNumStr]);
logstring (['Run: ' runNumStr]);

% Presents instructions to the subject
cgpencol(1,1,1);
cgfont(font_type,50);
cgtext('The run is about to begin. Please begin the task',0,100);
cgtext('as soon as you see the + and stop the task',0,50);
cgtext('when the + disappears. Please try not to make',0,0);
cgtext('unnecessary movements.',0,-50);
cgflip(0,0,0);
clearpict( 0 );

% Prepare crosshair in buffer window 2
loadpict( 'fixation_blackBG.bmp', 2 ); % Draw fixation point in display buffer 2

keys = getkeymap; % getkeymap command outputs a list of codes assigned to each key by cogent
clearkeys;
readkeys;
logkeys; %logs key presses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DUMMY SCAN                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%drawpict( 2 );  
%dummylength=repetitionTime*begnumdummies;
%wait(dummylength);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           STIMULUS PRESENT                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for a trigger signal from MRI scanner ('+') or press 'q' key
[keyPressed,scanTime]=waitkeydown(inf, [keys.PadAdd,keys.Escape,keys.Equals,keys.A,keys.Return]);
if keyPressed == keys.Escape     % if escape is pressed, abort experiment
    logstring('...aborted');  % print out error message and write to log
    stop_cogent;              % exits Cogent mode
    clear all;
	return;                    % exit from this function
end

for i = 1:numstim    
    % Forces cogent to stop and clears all variables
    % way to abort the localizer in case subject panics
    readkeys;
    n=getkeydown;
    idKeys=find(n==[keys.Escape] | n==[keys.Q]);
    if ~isempty(idKeys)     % if escape is pressed, abort experiment
        logstring('...aborted');  % print out error message and write to log
        stop_cogent;              % exits Cogent mode
        clear all;
        return;                    % exit from this function
    end
    settextstyle(font_type,font_size);
    preparestring( stim{i}, 1 ); %loads stimulus string in buffer window 1
    stimTime=drawpict( 1 ); %draws stimulus in buffer window 1
    waituntil(scanTime+totxdur(i)+stimdur(i));
    % logs stimulus and time
    logstring([stim{i},': ' num2str(stimTime)]);  
    drawpict( 2 ); %draws + in buffer window 2
    %waituntil(scanTime+dummylength+totxdur(i)); %waituntil is more precise!
    waituntil(scanTime+totxdur(i+1));
    clearpict(1);
end

xhairTime=drawpict( 2 );  
%waituntil(scanTime+dummylength+totrunDur);
waituntil(3000);
logstring(['End of run : ' num2str(xhairTime)]);

clearpict(0,background_color);

cgpencol(1,1,1);
cgfont(font_type,50);
cgtext('Thank You!',0,100);
cgtext('Please try not to make unnecessary movements',0,50); 
cgtext('while we prepare for the next run.',0,0);
cgflip(1,1,1);

[k,t,n] = waitkeydown(inf,[keys.Escape,keys.Q,keys.Return]); % Press Q or Esc to quit

%Forces cogent to stop and clears all variables
stop_cogent; %priority returned to normal
clear all;

return

