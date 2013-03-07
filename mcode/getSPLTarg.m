function splt=getSPLTarg(varargin)
    if (isempty(varargin)) % Transition prod
        splt=76;        
    else
        mouthMicDist=varargin{1};
        splt=76+20*log10(10/mouthMicDist);% assume the distance is x cm. it should be 75+20*log10(10/x)
    end
return