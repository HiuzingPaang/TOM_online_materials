function vectors = temreset(eventmat,rho,lambda,ndimension)

if nargin<4
    ndimension = 100;
end

    %% set parameters for model
    prob01 = rho;
    prob10 = rho;
    reset_lamda = lambda;
    nwholeitem = length(eventmat);
    idx_e = find(eventmat==1);
    if ismember(length(idx_e),1)
        with_boundary = 0;
    else
        with_boundary = 1;
        probtime = idx_e(2:end);
    end
    
    %% model part
    context = nan(nwholeitem,ndimension);
    p = rand(nwholeitem,ndimension);
    context(1,:) = round(rand(1,ndimension));
    
    for t = 2:nwholeitem
        idx0 = find(context(t-1,:)==0);
        idx1 = find(context(t-1,:)==1);
        if with_boundary == 0
            context(t,:) = context(t-1,:);    
            context(t,idx0(p(t,idx0)<=prob01)) = 1;
            context(t,idx1(p(t,idx1)<=prob10)) = 0;
        else
            tmp_context = context(t-1,:);     
            tmp_context(1,idx0(p(t,idx0)<=prob01)) = 1;
            tmp_context(1,idx1(p(t,idx1)<=prob10)) = 0;
            if any(t==probtime)
                tmp_reset_p = rand(1,ndimension);
                context(t,tmp_reset_p<=reset_lamda) = context(1,tmp_reset_p<=reset_lamda);
                context(t,tmp_reset_p>reset_lamda) = tmp_context(1,tmp_reset_p>reset_lamda);                    
            else         
                context(t,:) = tmp_context;
            end 
        end
    end

    vectors = context;
end