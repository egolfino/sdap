function pairs=randomizeTrials

conditionID = [1:8]';
sentenceID = [1:3]';

nTrials= 48;

nCons = numel(conditionID); 
nSentences = numel(sentenceID); 

randomizedCons=mod(randperm(nTrials),nCons)+1;
frandomizedCons=conditionID(randomizedCons);
randomizedSentences=zeros(size(randomizedCons));
for iCond=1:nCons-2
    idxCond=find(frandomizedCons==conditionID(iCond));
    randomizedSentences(idxCond)=mod([1:length(idxCond)],nSentences)+1;
end
pairs =[frandomizedCons randomizedSentences'];
