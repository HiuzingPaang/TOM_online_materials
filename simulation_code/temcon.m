function vectors = temcon(eventmat,rho,ndimension)

if nargin<3
    ndimension = 50;
end
    
    %% set parameters for model
    prob01 = rho;
    prob10 = rho;
    nwholeitem = length(eventmat);
    
    %% model part
    context = nan(nwholeitem,ndimension);
    p = rand(nwholeitem,ndimension);
    context(1,:) = round(rand(1,ndimension));
    
    for t = 2:nwholeitem
        context(t,:) = context(t-1,:);
        idx0 = find(context(t-1,:)==0);
        idx1 = find(context(t-1,:)==1);
        context(t,idx0(p(t,idx0)<=prob01)) = 1;
        context(t,idx1(p(t,idx1)<=prob10)) = 0;
    end
    
    vectors = context;
end