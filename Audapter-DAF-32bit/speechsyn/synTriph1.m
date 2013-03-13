function varargout=synTriph1()
    load('avgTraj_CL1');
    
    af1=fliplr(af1);    af2=fliplr(af2);
    
    dur=300;        % triphthong duration, msec
    sInt=2;         % msec
    
    tAxis1=linspace(0,dur,length(af1));
    tAxis2=0:sInt:dur;
    
    f1=interp1(tAxis1,af1,tAxis2);
    f2=interp1(tAxis1,af2,tAxis2);   
    
        
	load score_iau;		% gives score
	load defParsT_iau;	% gives defParsT
    
    tab=nan(0,size(score,2));
    for n=1:length(tAxis2)
        tab=[tab;N(end,:)];
        tab(end,1)=tAxis2(n);
        tab(end,6)=f1(n);
        tab(end,7)=f2(n);
    end
    
%     xlswrite('iau1.xls',sInt,'iau','D15');
%     xlswrite('iau1.xls',tab,'iau',['A83:','L',num2str(83+size(tab,1)-1)]);
    
    [w,fs]=slsyn('iau1.xls');
    w=double(w);
    w=w/max(abs(w))*0.5;
    figure();
    plot(w);
    wavplay(w,fs);
return