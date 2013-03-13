function [s,sRate] = slsyn(fileName,varargin);
%
% By Steven M. Lulich
% 4/23/2008
%
% For Speech Communication (HST.710) Lab 3

% Get the *.xls file
% [fileName, pathName] = uigetfile('*.xls','Choose *.xls synthesis file...');
% fileName = [pathName '/' fileName];

% Read the parameter data
[defParsT dummy1 dummy2] = xlsread(fileName,1,'D14:D76');
[score dummy1 dummy2] = xlsread(fileName,1,'A83:L10000');

clear dummy1 dummy2;

% Convert defParsT parameter data to proper format
defPars.DU = defParsT(1);
defPars.UI = defParsT(2);
defPars.SR = defParsT(3);
defPars.NF = defParsT(4);
defPars.F0 = defParsT(14);
defPars.TL = defParsT(18);
defPars.F1 = defParsT(23);
defPars.F2 = defParsT(27);
defPars.F3 = defParsT(29);
defPars.F4 = defParsT(31);
defPars.F5 = defParsT(33);
defPars.F6 = defParsT(35);
defPars.FTP = defParsT(41);
defPars.FTZ = defParsT(43);
defPars.A2F = defParsT(45);
defPars.A3F = defParsT(46);
defPars.A4F = defParsT(47);
defPars.A5F = defParsT(48);
defPars.A6F = defParsT(49);
defPars.B2F = defParsT(51);
defPars.B3F = defParsT(52);
defPars.B4F = defParsT(53);
defPars.B5F = defParsT(54);
defPars.B6F = defParsT(55);

clear defPars2;

% Create varPars variable
varPars = {'F0','AV','AH','AF','F1','F2','F3','F4','F5','F6','AI'};

% Synthesize the utterance
[s,sRate] = mlsyn(defPars,varPars,score);
% Output variables:  s, sRate
