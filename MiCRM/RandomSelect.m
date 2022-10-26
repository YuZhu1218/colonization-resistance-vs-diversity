function r = RandomSelect(Nm,p)
    r = zeros(1,Nm);
    random = randperm(Nm);
    r(random>p*Nm) = 1;
end
