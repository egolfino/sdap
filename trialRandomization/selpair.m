function d=selpair(b,s)
%The input b is a (s,2) matrix
%The input s is the number of tones
%The output d is the matrix b but without consecutive rows with 
%the same freq or amp, other detections are included
%just uncomment the line with the while you want and comment the others
%The returned value might not have the same size as the input 
%(less lines), this is because the function can't satisfy the condition,   
%it happens in the end of the data.

e=0;
for c=1:s
  d(c,:)=b(1+e,:);
  b(1+e,:)=[];e=0;

    %detect same freq
    %while((isempty(b)) | (d(c,1)==b(1+e,1)))

    %detect same amp
    %while((isempty(b)) | (d(c,2)==b(1+e,2)))

    %detect same freq and amp
    %while((isempty(b)) | (d(c,:)==b(1+e,:)))

    %detect same freq or amp
    while((isempty(b)) | ((d(c,1)==b(1+e,1))|(d(c,2)==b(1+e,2))))
        e=e+1;
        if ((e>=(numel(b(:,1))-1)) | (numel(b(:,1))<2)) % * error fixed *

            %Sometimes the function cannot satisfy the condition very near
            %the end of the data, some lines are ignore and the output size
            %is smaller than the input, if you wish the size to be always the
            %same just uncomment the next line (=no condition for that data)
            d=[d;b];
            idxBase=find(d==0);
            baseOrder=repmat([1:3],1,4)';
            d(idxBase)=baseOrder(randperm(12));
        return
        end
    end
end