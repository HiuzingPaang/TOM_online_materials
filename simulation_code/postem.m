function vectors = postem(eventmat,beta,ndimension)

if nargin<3
    ndimension = 50;
end
    
    %% set parameters for model
    l_prob01 = beta;
    l_prob10 = beta;
    nwholeitem = length(eventmat);
    idx_e = [find(eventmat==1),nwholeitem+1];
    nlocal = max(idx_e(2:end)-idx_e(1:end-1));

    %% model part
    context_only_local = nan(nlocal,ndimension);
    p_only_local = rand(nlocal,ndimension);
    context_only_local(1,:) = round(rand(1,ndimension));
    
    for t = 2:nlocal
        context_only_local(t,:) = context_only_local(t-1,:);
        idx0 = find(context_only_local(t-1,:)==0);
        idx1 = find(context_only_local(t-1,:)==1);
        context_only_local(t,idx0(p_only_local(t,idx0)<=l_prob01)) = 1;
        context_only_local(t,idx1(p_only_local(t,idx1)<=l_prob10)) = 0;
    end

    onlylocal = nan(nwholeitem,ndimension);
    for eg = 1:length(idx_e)-1
        e0 = idx_e(eg);
        e1 = idx_e(eg+1);
        onlylocal(e0:e1-1,:) = context_only_local(1:e1-e0,:);
    end
    
    vectors = onlylocal;
end