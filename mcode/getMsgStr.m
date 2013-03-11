function msg = getMsgStr(thisphase)
if isnumeric(thisphase)
    msg={'The run is about to begin. Please read', ...
         'the sentence out loud as soon as it', ...
         'appears. Do NOT read anything after the', ...
         '+ appears. Try not to move too much.'};
else
    switch(thisphase)
        case 'pre'
            msg={'Welcome to the experiment!',...
                '',...
                'When you are ready, please press the "play" button to start.',...
                };
        case 'pract1'         % pract1
            msg={'Now, you will see a volume meter below the word.',...
                '',...
                'It will tell you how loud you''ve spoken. Please carefully',...
                'adjust you loudness to hit the middle green area. ',...
                '',...
                'A trial will be repeated until you successfully land in this area.',...
                '',...
                'Press "play" to continue.',...
                };
        case 'pract2'		% pract2
            msg={'Good job!',...
                '',...
                'Now, in addition to the volume meter, a speed meter is',...
                'displayed. It will tell you how fast you''ve spoken ',...
                'Please carefully adjust your word length to hit the middle green area',...
                '',...
                'A trial will be repeated until you simultaneously lands in the green ',...
                'areas of both the volume and speed meters. Press "play" to continue.'};
        case 'start'	% start
            msg={'Good job!',...
                '',...
                'Now, the trials will not repeat even if you fall out of the green areas.',...
                'However, please try to keep the volume and speed you''ve learned ',...
                'throughout the rest of the experiment.',...
                '',...
                'Press "play" to continue'};
        case 'ramp'
            msg={''};
        otherwise,
            msg={''};
end
end
    
