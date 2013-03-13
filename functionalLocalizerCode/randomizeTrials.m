function [randomized, separation]=randomizeTrials

    numconditions = 4;
    repconditions = 6;

    totconditions = numconditions*repconditions;
    matrixconditions = repmat([1:numconditions],1,repconditions);

    randomid=randperm(totconditions);

    randomized=matrixconditions(randomid);
    
    interval=1000/(totconditions+1);
    steps=interval.*[1:totconditions];
    
    separation=2500+steps;
    separation=separation(randomid);

end