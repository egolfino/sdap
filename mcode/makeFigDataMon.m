function figIdDat=makeFigDataMon2()
fid=figure('Position',[300,200,900,640],'Name','TransAdapt: Data monitor');
axes1=subplot('Position',[0.05,0.575,0.275,0.4]);   % Input waveform
axes2=subplot('Position',[0.05,0.1,0.275,0.4]);   % Output waveform
axes3=subplot('Position',[0.40, 0.575, 0.525, 0.4]);
axes4=subplot('Position',[0.40, 0.1, 0.525,0.4]);
% axes5=subplot('Position',[0.7,0.575,0.275,0.4]);
% axes6=subplot('Position',[0.7,0.1,0.275,0.4]);
% figIdDat=[fid,axes1,axes2,axes3,axes4,axes5,axes6];
figIdDat=[fid,axes1,axes2,axes3,axes4];
return