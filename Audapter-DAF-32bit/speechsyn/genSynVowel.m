function varargout=genSynVowel(varargin)
%% Add paths
addpath('./');
addpath('../mcode/commonmcode');
addpath('../mcode/commonmcode/graph');

%%
    vF0=[240,240];
    vF1=[600,400];
    vF2=[650*3,2200];
    
    baseFreq=650;
    if ~isempty(findStringInCell(varargin,'baseFreq'))
        baseFreq=varargin{findStringInCell(varargin,'baseFreq')+1};
    end
    
    vF3=[baseFreq*5,baseFreq*5];
    vF4=[baseFreq*7,baseFreq*7];
    vF5=[baseFreq*9,baseFreq*9];
    vF6=[baseFreq*11,baseFreq*11];
    fileName='./vowel';
    
    if ~isempty(findStringInCell(varargin,'vF0'))
        vF0=varargin{findStringInCell(varargin,'vF0')+1};
    end    
    if ~isempty(findStringInCell(varargin,'vF1'))
        vF1=varargin{findStringInCell(varargin,'vF1')+1};    
    end
    if ~isempty(findStringInCell(varargin,'vF2'))
        vF2=varargin{findStringInCell(varargin,'vF2')+1};    
    end
    
    if ~isempty(findStringInCell(varargin,'fileName'))
        fileName=varargin{findStringInCell(varargin,'fileName')+1};
    end
    
    toPlay=1;
    if ~isempty(findStringInCell(varargin,'toPlay'))
        toPlay=varargin{findStringInCell(varargin,'toPlay')+1};
    end
%%    
	fs=48000;
% 	
	sInt=2;		% msec
	dur1=400;	% msec
	dur0=50;
	dur=dur1+2*dur0;

    tAxis2=0:sInt:dur;

	load score_iau;		% gives score
    % Column 1: time;       % Column 2: F0 (0.1*Hz)
    % Column 3: AV          % Column 4: AH
    % Column 5: AF          % Column 6: F1 (Hz)
    % Column 7: F2 (Hz)     % Column 8: F3 (Hz)
    % Column 9: F4 (Hz)     % Column 10: F5 (Hz)
    % Column 11: F6 (Hz)    % Column 12: AI
    
	load defParsT_iau;	% gives defParsT
    tab=nan(0,size(score,2));
    for n=1:length(tAxis2)
        tab=[tab;score(end,:)];
        tab(end,1)=tAxis2(n);
        tab(end,2)=10*(vF0(1)+(vF0(2)-vF0(1))*(n-1)/(length(tAxis2)-1));
        tab(end,6)=vF1(1)+(vF1(2)-vF1(1))*(n-1)/(length(tAxis2)-1);
        tab(end,7)=vF2(1)+(vF2(2)-vF2(1))*(n-1)/(length(tAxis2)-1);
        tab(end,8)=vF3(1)+(vF3(2)-vF3(1))*(n-1)/(length(tAxis2)-1);
        tab(end,9)=vF4(1)+(vF4(2)-vF4(1))*(n-1)/(length(tAxis2)-1);
        tab(end,10)=vF5(1)+(vF5(2)-vF5(1))*(n-1)/(length(tAxis2)-1);
        tab(end,11)=vF6(1)+(vF6(2)-vF6(1))*(n-1)/(length(tAxis2)-1);
	end	
	
%     xlswrite('iau1.xls',sInt,'iau','D15');
%     xlswrite('iau1.xls',tab,'iau',['A83:','L',num2str(83+size(tab,1)-1)]);	
	varPars = {'F0','AV','AH','AF','F1','F2','F3','F4','F5','F6','AI'};		
	defPars.DU = defParsT(1);	defPars.UI = defParsT(2);
	defPars.SR = defParsT(3);	defPars.NF = defParsT(4);
	defPars.F0 = defParsT(14);	defPars.TL = defParsT(18);
	defPars.F1 = defParsT(23);	defPars.F2 = defParsT(27);
	defPars.F3 = defParsT(29);	defPars.F4 = defParsT(31);
	defPars.F5 = defParsT(33);	defPars.F6 = defParsT(35);
	defPars.FTP = defParsT(41);	defPars.FTZ = defParsT(43);
	defPars.A2F = defParsT(45);	defPars.A3F = defParsT(46);
	defPars.A4F = defParsT(47);	defPars.A5F = defParsT(48);
	defPars.A6F = defParsT(49);	defPars.B2F = defParsT(51);
	defPars.B3F = defParsT(52);	defPars.B4F = defParsT(53);
	defPars.B5F = defParsT(54);	defPars.B6F = defParsT(55);
	
	defPars.DU=dur;
	
	[w,fs0]=mlsyn(defPars,varPars,tab);
    w=double(w);
    w=w/max(abs(w))*0.5;
%     figure();
%     plot(w);

	if (size(w,1)>size(w,2))	w=w';	end

	if (~isempty(findStringInCell(varargin,'level')))
		t1=round(tAxis2(Nub)*fs0/1e3);
		t2=round(tAxis2(Nlb)*fs0/1e3)-1;
        tRMSOut=calcAWeightedRMS(w(t1:t2),fs0);
        soundLv1=20*log10(tRMSOut/(dBSPL2WaveAmp(0,1000,pcrKnob)/sqrt(2)));  % dBA SPL        
        gain=10^((level-soundLv1)/20);
    else
        gain=1;
	end
	
	w=w*gain;
	w=resample(w,fs,fs0);
	
%%
    if toPlay
        wavplay(w,fs);
    end
    wavwrite(w,fs,[fileName,'.wav']);
    t=tAxis2 * 1e-3;    % s
    f0=tab(:,2)/10;
    f1=tab(:,6);
    f2=tab(:,7);
    save([fileName,'.mat'],'t','f0','f1','f2');
return
