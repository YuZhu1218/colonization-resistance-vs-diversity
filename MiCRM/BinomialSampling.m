function C = BinomialSampling(Nc,Nm,q)

C = zeros(Nc,Nm);
rndc = rand(Nc,Nm);
C(rndc <= q) = 1;

return;