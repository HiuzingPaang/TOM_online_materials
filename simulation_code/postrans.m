function vectors = postrans(eventmat,ndimension,nscaler)

if nargin<2
    ndimension = 50;
    nscaler = 10000;  % User defined scalar, default to 10000 by the authors of Attention is all u need
end

if nargin<3
    nscaler = 10000;  % User defined scalar, default to 10000 by the authors of Attention is all u need
end
    
    %% set parameters for model
    nwholeitem = length(eventmat);
    idx_e = [find(eventmat==1),nwholeitem+1];
    npos = max([length(idx_e)-1,idx_e(2:end)-idx_e(1:end-1)]);
    
    %% model part
    pos_only = nan(npos,ndimension);
    
    for k = 1:npos
        for i =1:floor(ndimension/2)
            denominator = nscaler^(2*i/ndimension);
            pos_only(k,2*i-1) = sin(k/denominator);
            pos_only(k,2*i) = cos(k/denominator);
        end
        if mod(ndimension,2)==1
            denominator = nscaler^(2*(i+1)/ndimension);
            pos_only(k,end) = sin(k/denominator);
        end
    end
   
    onlylocal = nan(nwholeitem,ndimension);
    for eg = 1:length(idx_e)-1
        e0 = idx_e(eg);
        e1 = idx_e(eg+1);
        onlylocal(e0:e1-1,:) = pos_only(1:e1-e0,:);
    end
    
    vectors = onlylocal;
end