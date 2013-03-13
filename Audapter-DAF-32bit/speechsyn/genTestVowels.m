function genVowelsFemale()
%% Add paths
addpath('./');
addpath('../mcode/commonmcode');
addpath('../mcode/commonmcode/graph');

%% Vowel information
vowelName=  {'IY','IH','EH','AE','AH','AA','AO','UW','UH','ER','AY','EY','AW','OW'};

vowelF1(1,:)=[288, 446, 650, 856, 795, 880, 570, 311, 480, 490, 1000,630, 1000,570];
vowelF1(2,:)=[288, 446, 650, 856, 795, 880, 570, 311, 480, 490, 500, 420, 420, 311];
vowelF2(1,:)=[2700,2500,2400,1900,1450,1250,800, 820, 1200,1600,1300,2000,1200,800];
vowelF2(2,:)=[2700,2500,2400,1900,1450,1250,800, 820, 1200,1600,2700,2700,820, 780];

F0Mids.female=160:20:300;
F0Mids.male=90:10:160;

toPlay=0;

baseFreq.female=650;
baseFreq.male=550;

%TODO: Interpolate between the point vowels

%% 
for sexes={'male','female'}
    sex=sexes{1};
    for k=1:length(F0Mids.(sex))
        F0Mid=F0Mids.(sex)(k);
        F0Low=F0Mid/1.2;
        F0High=F0Mid*1.2;
        for n=1:length(vowelF1)
            fn=['../testVowels/',sex(1),'_f0-',num2str(F0Mid),'-',num2str(F0Mid),...
                '_f1-',num2str(vowelF1(1,n)),'-',num2str(vowelF1(2,n)),...
                '_f2-',num2str(vowelF2(1,n)),'-',num2str(vowelF2(2,n))];
            genSynVowel('baseFreq',baseFreq.(sex),'vF1',vowelF1(:,n),'vF2',vowelF2(:,n),...
                'vF0',repmat(F0Mid,1,2),'fileName',fn,'toPlay',toPlay);

            fn=['../testVowels/',sex(1),'_f0-',num2str(F0Mid),'-',num2str(F0Low),...
                '_f1-',num2str(vowelF1(1,n)),'-',num2str(vowelF1(2,n)),...
                '_f2-',num2str(vowelF2(1,n)),'-',num2str(vowelF2(2,n))];
            genSynVowel('baseFreq',baseFreq.(sex),'vF1',vowelF1(:,n),'vF2',vowelF2(:,n),...
                'vF0',[F0Mid,F0Low],'fileName',fn,'toPlay',toPlay);

            fn=['../testVowels/',sex(1),'_f0-',num2str(F0Mid),'-',num2str(F0High),...
                '_f1-',num2str(vowelF1(1,n)),'-',num2str(vowelF1(2,n)),...
                '_f2-',num2str(vowelF2(1,n)),'-',num2str(vowelF2(2,n))];
            genSynVowel('baseFreq',baseFreq.(sex),'vF1',vowelF1(:,n),'vF2',vowelF2(:,n),...
                'vF0',[F0Mid,F0High],'fileName',fn,'toPlay',toPlay);
        end
    end
end
return