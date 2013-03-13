function varargout=genSynTriph(r,dir,varargin)
%%
    if (~isempty(findStringInCell(varargin,'level')))
        level=varargin{findStringInCell(varargin,'level')+1};
        pcrKnob=varargin{findStringInCell(varargin,'level')+2};
	end
	
%%
	fs=48000;
	
	sInt=2;		% msec
	dur1=400;	% msec
	dur0=50;
	dur=dur1+2*dur0;
	f2min=1000;	f2max=2000;	% Hz
	f2ub=1900;	f2lb=1100;	% Hz
	f1min=500;	f1max=750;	% Hz
	
	fieldWidth=f2ub-f2lb;
	maxPert=0.25*fieldWidth;
	
	N1=round(dur1/sInt)+1;
	N0=round(dur0/sInt);
	N=N1+2*N0;
	
	f2s=linspace(f2max,f2min,N1);
	p=polyfit([1,round(N1/2)+1,N1],[500,750,500],2);
	f1s=polyval(p,1:N1);
	
	f1s=[f1s(1)*ones(1,N0),f1s,f1s(end)*ones(1,N0)];
	f2s=[f2s(1)*ones(1,N0),f2s,f2s(end)*ones(1,N0)];
    
    f1s_o=f1s;
    f2s_o=f2s;
	
    Nslope=[min(find(diff(f2s)<0)):1:max(find(diff(f2s)<0))+1];    
    
	Nub=interp1(f2s(N0+1:N0+N1),N0+1:N0+N1,f2ub);
	Nlb=interp1(f2s(N0+1:N0+N1),N0+1:N0+N1,f2lb);
	Nfield=[Nub,ceil(Nub):1:floor(Nlb),Nlb];       
	if (Nfield(1)==Nfield(2))
		Nfield=Nfield(2:end);
	end
	if (Nfield(end)==Nfield(end-1))
		Nfield=Nfield(1:end-1);
	end
	
	% The perturbation field
	p2=polyfit([Nub,mean([Nub,Nlb]),Nlb],[0,maxPert*r,0],2);
	pv=polyval(p2,Nfield);
	        
	if (isequal(dir,'f2_up'))
		f2s(Nfield)=f2s(Nfield)+pv;
    elseif (isequal(dir,'f2_down'))
		f2s(Nfield)=f2s(Nfield)-pv;
    elseif (isequal(dir,'AccelDecel') | isequal(dir,'DecelAccel'))
        for n=1:length(pv)
            if (isequal(dir,'AccelDecel'))
                f2s(Nfield(n))=f2s(Nfield(n))-pv(n);
            else
                f2s(Nfield(n))=f2s(Nfield(n))+pv(n);
            end
            if (f2s(Nfield(n))>max(f2s_o(Nslope)))
                f1s(Nfield(n))=interp1(f2s_o(Nslope),f1s_o(Nslope),max(f2s_o(Nslope)));
            elseif (f2s(Nfield(n))<min(f2s_o(Nslope)))
                f1s(Nfield(n))=interp1(f2s_o(Nslope),f1s_o(Nslope),min(f2s_o(Nslope)));
            else
                f1s(Nfield(n))=interp1(f2s_o(Nslope),f1s_o(Nslope),f2s(Nfield(n)));
            end

        end
    end

    tAxis2=0:sInt:dur;
    
    if (~isempty(findStringInCell(varargin,'plot')))
%         plot(tAxis2,f1s,'r');  hold on;    plot(tAxis2,f2s,'r');
        figure; set(gca,'FontSize',12);
        plot(tAxis2,f1s_o,'b','LineWidth',1.5); hold on;
        plot(tAxis2,f2s_o,'b','LineWidth',1.5);
        plot([tAxis2(1),tAxis2(end)],[f2ub,f2ub],'k--');
        plot([tAxis2(1),tAxis2(end)],[f2lb,f2lb],'k--');
        set(gca,'YLim',[0,2500]);
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
%         plot(f1s,f2s,'r');  hold on;
%         plot(f1s_o,f2s_o,'k');        
    end
% 	clf;
% 	plot(f1s,f2s);
%%
	

% 	[N1,T]=xlsread('iau1.xls');
	load score_iau;		% gives score
	load defParsT_iau;	% gives defParsT
    tab=nan(0,size(score,2));
    for n=1:length(tAxis2)
        tab=[tab;score(end,:)];
        tab(end,1)=tAxis2(n);
        tab(end,6)=f1s(n);
        tab(end,7)=f2s(n);
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
	t0=1/fs;
	t1=0:t0:dur0/1e3;
	env=sin(pi/2*t1/(dur0/1e3)).^2;
	env=[env,ones(1,length(w)-2*length(env)),fliplr(env)];
	w=w.*env;

%% Output	
	if (nargout==1)
		varargout{1}=w;
	elseif (nargout==2)
		varargout{1}=w;
		varargout{2}=48000;	% Hz
	end
	
%%
	if (~isempty(findStringInCell(varargin,'play')))
		if (length(w)<120000)
			w=[w,zeros(1,120000-length(w))];
		end	
		TransShiftMex(3,'datapb',w);
		TransShiftMex(12);
		pause(2.5);
		TransShiftMex(2);
	end
return
