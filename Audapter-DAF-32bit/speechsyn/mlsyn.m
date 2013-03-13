%MLSYN  - matlab interface to Klatt synthesizer
%
%	usage:  [s, sRate] = mlsyn(defPars, varPars, score)
%
% This procedure is an interface to the cascade/parallel formant synthesizer of 
% Dennis Klatt, based on SYN version 2.1 (98/07/17)
%
% DEFPARS is a struct with fieldnames matching KLSYN parameter names whose values
% overwrite defaults; for a list of parameter defaults enter "mlsyn defaults"
%
% VARPARS is a cellstr vector of names [1 x nParams] specifying the time-varying parameters
%
% SCORE is a [nSteps x nParams+1] matrix of which the first column gives the time offset
% in ms at which the specified changes are to become active, and the remaining columns
% (in VARPARS order) give the values to replace (all values persist until replaced)
%
% returns the vector of synthesized samples S [nSamps x 1] (int16 format)
% optionally returns sampling rate SRATE
%
% see also PARSEKLSYNDOC, which loads parameters from KLSYN DOC files
%
%Examples
%
%>> mlsyn defaults                           % list parameter defaults
%>> soundsc(double(mlsyn),10000);            % generate and play the default utterance
%
%>> defPars = struct('DU',750, 'SR',8000);   % change duration and sampling rate
%>> varPars = {'F0','AV'};  % vary F0 (100 and 200 Hz) and voicing amplitude
%>> score = [0 0 0 ;100 1000 50; 300 2000 50; 500 0 0];   % voicing from 100-500 ms
%>> s = mlsyn(defPars, varPars, score);      % synthesize
%
%>> [defPars, varPars, score] = ParseKLSYNDoc('FOO');     % load params from FOO.DOC
%>> [s,sRate] = mlsyn(defPars, varPars, score);           % synthesize FOO
%

% mkt 03/07
