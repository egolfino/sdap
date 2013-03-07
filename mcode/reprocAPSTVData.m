function dataOut=reprocAPSTVData(data,varargin)
p=data.params;
sigIn=data.signalIn;

if ~isempty(varargin)
    for i1=1:2:length(varargin)
        p.(varargin{i1})=varargin{i1+1};
    end
end
    
MexIO('reset');
MexIO('init',p);

sigIn=resample(sigIn,48000,data.params.sr);
sigInCell=makecell(sigIn,64);
for n = 1 : length(sigInCell)
    tic;
    TransShiftMex(5,sigInCell{n});
end

dataOut=MexIO('getData');

return