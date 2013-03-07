function test_ds3
addpath e:\speechres\sdap\Audapter-DAF\BIN\debug

DOWNSAMP_FACT = 3;
TransShiftMex(3, 'srate', 48000 / DOWNSAMP_FACT);
TransShiftMex(3, 'framelen', 96 / DOWNSAMP_FACT);
TransShiftMex(3, 'bpitchshift', 1);
TransShiftMex(3, 'pvocframlen', 256);
TransShiftMex(3, 'pvochop', 64);

TransShiftMex(6);
TransShiftMex(3, 'pitchshiftratio', 1);
TransShiftMex(1);
pause(8);
TransShiftMex(2);

TransShiftMex(6);
TransShiftMex(3, 'pitchshiftratio', 1.1);
TransShiftMex(1);
pause(8);
TransShiftMex(2);

TransShiftMex(6);
TransShiftMex(3, 'pitchshiftratio', 1);
TransShiftMex(1);
pause(8);
TransShiftMex(2);


sig = TransShiftMex(4);
return