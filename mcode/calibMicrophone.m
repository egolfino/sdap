function calibMicrophone
%%
calibWavFN = 'calib_vowel.wav';
if ~isfile(calibWavFN)
    error('Cannot find the wav file for calibration: %s\n', calibWavFN);
end

%%
TransShiftMex(0);
TransShiftMex(3, 'srate', 16000);
TransShiftMex(3, 'framelen', 32);

[w, fs] = wavread(calibWavFN);

foo=input('Press any enter to continue...','s');

% soundsc(w, fs);

foo=input('Wait for SLM reading to stabilize. Press any enter to start TransShiftMex when ready...','s');

TransShiftMex(1);
pause(0.05);
TransShiftMex(2);

TransShiftMex(1);
pause(2);
TransShiftMex(2);
sig=TransShiftMex(4);
sig=sig(:,1);

clear functions;

figure;
sr=12e3;
taxis=0:(1/sr):((length(sig)-1)/sr);
plot(taxis,sig); 
hold on;



title('Set tStart!');
coord1=ginput(1);
tStart=coord1(1);
ys=get(gca,'YLim');
plot(repmat(tStart,1,2),ys,'k--');

title('Set tEnd');
coord2=ginput(1);
tEnd=coord2(1);
plot(repmat(tEnd,1,2),ys,'k-');

sig_sel=sig(find(taxis>=tStart & taxis<=tEnd));

level_SLM=input('level_SLM = ');

micRMS_100dBA=10^((100-level_SLM)/20)*rms(sig_sel);

fprintf('rms(sig_sel) = %.5f\n', rms(sig_sel));
fprintf(1, 'micRMS_100dBA target is 0.1 [0.09 - 0.11]\n');
fprintf('micRMS_100dBA = %.5f\n', micRMS_100dBA);