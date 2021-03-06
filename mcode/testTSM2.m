function varargout = testTSM2(data,p)    
%% 
    fs = data.params.sr;
    
    sigIn = data.signalIn;

    sigIn = resample(sigIn, 48000, fs);     
    sigInCell = makecell(sigIn, 64);    
    
    TransShiftMex(6);   % Reset;\
	
    p.rmsClipThresh=0.01;
	p.bRMSClip=1;
	
    MexIO('init',p);

    for n = 1 : length(sigInCell)
        TransShiftMex(5, sigInCell{n});
	end
        
	data=MexIO('getData');
    
    [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
    [k1,k2]=detectVowel(data.fmts(:,1),data.fmts(:,2),iv1,iv2,'eh','rms',data.rms(:,1));

%% Output
    if (nargout==1)
        varargout{1}=data.fmts(:,1:2);
    end
return