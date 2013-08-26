function testFeedback()
%%
configFN = '../expt/expt_config.txt';

%%

expt_config = read_parse_expt_config(configFN);

expt_config.MOUTH_MIC_DIST = input('Please input mouth-mic distance (cm): ');

p = getTSMDefaultParams(expt_config.SUBJECT_GENDER, ...
                            'DOWNSAMP_FACT', expt_config.DOWNSAMP_FACT, ...
                            'FRAME_SIZE', expt_config.FRAME_SIZE / expt_config.DOWNSAMP_FACT, ...
                            'closedLoopGain', expt_config.CLOSED_LOOP_GAIN, ...
                            'mouthMicDist', expt_config.MOUTH_MIC_DIST);
MexIO('init', p);

%%
TransShiftMex(3, 'fb', 1);
MexIO('reset');

TransShiftMex(1);
foo = input('Hit Enter to stop ...');
TransShiftMex(2);

%% Load the multi-talker babble noise
[x_mtb, fs_mtb]=wavread('mtbabble48k.wav');

x_mtb = x_mtb - mean(x_mtb);
mb_rms = rms(x_mtb);
x_mtb = x_mtb / mb_rms;

TransShiftMex(3, 'datapb', x_mtb, 0);
TransShiftMex(3, 'rmsff_fb', ...
              [expt_config.SMN_KERNEL_0, expt_config.SMN_KERNEL_1, ...
               expt_config.SMN_ON_RAMP, expt_config.SMN_OFF_RAMP]);
TransShiftMex(3, 'fb4gaindb', expt_config.SMN_GAIN);
% TransShiftMex(3, 'fb4gaindb', 50);
TransShiftMex(3, 'fb', 4);

MexIO('reset');
foo = input('Hit Enter to start the SMN mode ...');

TransShiftMex(1);
foo = input('Hit Enter to stop ...');
TransShiftMex(2);

sig = TransShiftMex(4);
figure;
plot(sig(:, 1));
return