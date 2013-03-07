function calibMicrophone_EMMA
TransShiftMex(2);
MexIO('reset');

foo=input('Press any enter to continue','s');
MexIO('reset');
TransShiftMex(1);
pause(18);
TransShiftMex(2);

sig=TransShiftMex(4);
sig=sig(:,1);
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

fprintf('rms(sig_sel) = %.5f\n',rms(sig_sel))
fprintf('micRMS_100dBA = %.5f\n',micRMS_100dBA);


return