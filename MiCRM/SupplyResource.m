function R0 = SupplyResource(nMediator,Ra0)
    R0 = zeros(nMediator,1);
    r = randi(nMediator);
    s = randi(nMediator-1);
    R0(r) = Ra0;
    if r>=s 
        R0(s) = Ra0;
    elseif r < s
        R0(s+1) = Ra0;
    end
return