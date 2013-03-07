function processNimStimFaces(varargin)
rawImageDir='C:\Users\scai\nimStimFaces';
procImagedir='C:\speechres\sap-fmri\mcode\graphics\faces';

d{1}=dir(fullfile(rawImageDir,['*_ca_*.bmp']));
d{2}=dir(fullfile(rawImageDir,['*_ne_*.bmp']));
d{3}=dir(fullfile(rawImageDir,['*_ha_c*.bmp']));

faceCnt=1;

recSubjIDs=cell(1,0);
recSubjSex=cell(1,0);
cntSubjID=0;

for m=1:length(d)
	for n=1:length(d{m})
        subjID=d{m}(n).name(1:3);
        subjID=upper(subjID);
        if (isempty(findStringInCell(recSubjIDs,subjID)))
            recSubjIDs{length(recSubjIDs)+1}=subjID;
            recSubjSex{length(recSubjSex)+1}=subjID(end);
            cntSubjID=length(recSubjIDs);
            subjSex=recSubjSex{end};
        else
            cntSubjID=findStringInCell(recSubjIDs,subjID);
            subjSex=recSubjSex{cntSubjID};
        end
        
        
		im=imread(fullfile(rawImageDir,d{m}(n).name));
        im1=imresize(im,[120,100]);
        fn=['face',num2str(faceCnt),'-s',num2str(cntSubjID),'-',subjSex,'.bmp'];
		imwrite(im1,fullfile(procImagedir,fn));
		disp([fullfile(procImagedir,fn),' written.']);
		faceCnt=faceCnt+1;
	end
end
return