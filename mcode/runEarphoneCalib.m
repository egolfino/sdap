function runEarphoneCalib
    TransShiftMex(2);
    
    calibFreqs=[500,1000,2000,3000,4000,5000,6000];
    
    for n=1:length(calibFreqs)
        disp(['f = ',num2str(calibFreqs(n)),'Hz']);
        TransShiftMex(3,'wgfreq',calibFreqs(n));
        TransShiftMex(11);
        
        disp('Press Enter to proceed to the next frequency');
        pause
        TransShiftMex(2);
    end
    
    TransShiftMex(2);
return