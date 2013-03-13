function pairs=randomizeTrials

conditionID = [1:8]';
sentenceID = [1:5]';

nTrials= 40;

nCons = numel(conditionID); 
nSentences = numel(sentenceID); 

randomizedCons=mod(randperm(nTrials),nCons)+1;
frandomizedCons=conditionID(randomizedCons);
randomizedSentences=zeros(size(randomizedCons));
for iCond=1:nCons-2
    idxCond=find(frandomizedCons==conditionID(iCond));
    randomizedSentences(idxCond)=randperm(nSentences);
end
pairs =[frandomizedCons randomizedSentences'];
