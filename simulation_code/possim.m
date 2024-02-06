function vectors = possim(eventmat,beta,ndimension)

if nargin<3
    ndimension = 50;
end

    
    %% set parameters for model
    l1 = beta;
    l = sqrt(1-l1^2);
    nwholeitem = length(eventmat);
    idx_e = [find(eventmat==1),nwholeitem+1];
    nlocal = max(idx_e(2:end)-idx_e(1:end-1));

	%% model part
    l_only = nan(ndimension,ndimension);
    
    for i = 1:ndimension
        l_only(i,1) = l^(i-1);
        if i == 1
            l_only(1,2:end) = 0;
        else
            for j = 1:i-1 
                l_only(i,j+1) = l^(i-j-1)*l1;
            end
        end
        if i > 1 && i < ndimension
            l_only(i,i+1:end) = 0;
        end
    end
    l_only_vector = l_only(ndimension-nlocal+1:ndimension,:);
    
    onlylocal = nan(nwholeitem,ndimension);
    for eg = 1:length(idx_e)-1
        e0 = idx_e(eg);
        e1 = idx_e(eg+1);
        onlylocal(e1-1:-1:e0,:) = l_only_vector(end-(e1-e0)+1:end,:);
    end
    
    vectors = onlylocal;
end