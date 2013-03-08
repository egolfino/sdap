function testpb
srate = 48000 / 4;
framelen = 96 / 4;

TransShiftMex(3, 'srate', srate);
TransShiftMex(3, 'framelen', framelen);
TransShiftMex(6);

load('E:\STUT_DATA\PFS_M04\APSTV2_STUT_EH\main\rep4\trial-9-6.mat');
x = resample(data.signalIn, 48e3, 12e3);

TransShiftMex(3, 'datapb', x);

TransShiftMex(12);
pause(2.5);
TransShiftMex(2);
return