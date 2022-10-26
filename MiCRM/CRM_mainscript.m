  % MiCRM model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R
% CC: keeping track of resident community composition after attempted invasion
        
clear
rndseed0 = 6130; 

nSample = 10000; % # of samples being screened
rndseed = round(100*nSample*rand(1,nSample));  

InvFrac = 1e-3;
nGen = 200;
nInitialCell = 1e4; % total initial cells
dilTh = 1e7; % coculture dilution threshold
GenPerRound = log(dilTh/nInitialCell)/log(2);
nRound = round(nGen/GenPerRound); % number of rounds of propagation
extTh = 0.1; % population extinction threshold

nCellType = 200; % # of cell types in the initial pool
nMediator = 20; % # of mediators
fss = 0.1; % fraction of the supply resource
R0 = 1000; % total amount of supply resource
Ra0 = R0 / (fss * nMediator); % fractional amount of supply resource

tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
del = 0; % avg. decay rate (per hour)
taust = 1; % timecale for supply of external resource
m = 1; % threshold consumption level for growth
l = 0.5; % fraction of resources secreted as by-product

% parameters for cia
qa = 0.5; % the probability of that resources A is preferrably consumed by family A
qs = 0.5; % the probability of that resources S is preferrably consumed by family S
muc = 10; % mean of cia
sigmac = 3; % variance of cia
pasn = 0.5; % the probability that species belong to family A
pasr = 0.5; % the probability that resources belong to family A

% parameters for Dab
s = 0.3; 

R0T = zeros(nMediator,nSample);
dabT = zeros(nMediator,nMediator,nSample);
ciaT = zeros(nCellType,nMediator,nSample);

CmpDT = zeros(nCellType,nSample);
CmpDCT = zeros(nCellType, nSample);

NE0D = zeros(1,nSample);
NE0DC = zeros(1,nSample);

CmpFDIF = zeros(nCellType+1,nSample);
InvEffD1 = zeros(nSample);

for ns = 1:nSample
    disp(ns)
    tic
    
    rng(rndseed(ns),'twister');
    
    % assign species and resources to A and S
    P = repmat(BinomialSampling(nCellType,1,pasn),1,nMediator);
    R = repmat(BinomialSampling(1,nMediator,pasr),nCellType,1);
    
    %%Parameters
    % cia
    cia_m = FamilyEncounter(P, R, qa, qs).* muc / nMediator;
    cia_v = FamilyEncounter(P, R, qa, qs).* sigmac / nMediator;
    sigma = cia_m./cia_v;
    alpha = (cia_m.^2)./cia_v;
    cia = gamrnd(alpha, sigma,nCellType,nMediator);
    
    % Dab
    conc = zeros(1,nMediator);
    conc(:,:) = 1/(s*nMediator);
    dab = DirichletSampling(conc,nMediator)';
    
    % R0
    R0 = SupplyResource(nMediator,Ra0);
    
    CmpB = 1 / nCellType * ones(1,nCellType);
    [Ne, CmpF, Ne0, CmpF0] = WellmixedInteraction(nRound,CmpB,nInitialCell,cia,l,m,R0,dab,taust,tauf,dtau,dilTh,extTh);

    V0 = zeros(1,nCellType);
    V0D(Ne0) = 1;
    
    V0DC = zeros(1,nCellType);
    V0DC(Ne) = 1;
    
    NE0D(:,ns) = sum(V0D,2);
    NE0DC(:,ns) = sum(V0DC,2);
    
    Cmp0D = zeros(1,nCellType);
    Cmp0D(Ne0) = CmpF0;

    Cmp0DC = zeros(1,nCellType);
    Cmp0DC(Ne) = CmpF;
    
    CmpDT(:,ns) = Cmp0D;
    CmpDCT(:,ns) = Cmp0DC;
    
    R0T(:,ns) = R0;
    dabT(:,:,ns) = dab;
    ciaT(:,:,ns) = cia;
    
    % Properties of the invader
    RI = BinomialSampling(1,nMediator,pasr);
    PI = repmat(BinomialSampling(1,1,pasn),1,nMediator);
    cia_mI = FamilyEncounter(PI, RI, qa, qs).* muc / nMediator;
    cia_vI = FamilyEncounter(PI, RI, qa, qs).* sigmac^2 / nMediator;
    sigmaI = cia_mI./cia_vI;
    alphaI = (cia_mI.^2)./cia_vI;
    ciaI = [cia;gamrnd(alphaI, sigmaI,1,nMediator)];
    
    % Compostion when invading
    Cmp0DI = [(1-InvFrac)*Cmp0D, InvFrac];
    [NeI, CmpFI, Ne0I, CmpF0I] = WellmixedInteraction(nRound,Cmp0DI,nInitialCell,ciaI,l,m,R0,dab,taust,tauf,dtau,dilTh,extTh);
    
    CmpFDIF(NeI,ns) = CmpFI;
    InvEffD1(ns) = CmpFDIF(nCellType+1,ns)/InvFrac;
end

save(strcat('CRM_V1_fss',num2str(round(100*fss)),'_taust',num2str(round(100*taust)),'_CSD',num2str(nInitialCell,2),'_DilTh',num2str(dilTh,2),'_ExtTh',num2str(extTh),'_l',num2str(l),'_m',num2str(m),'_Ra0',num2str(Ra0),'_qa',num2str(round(100*qa)),'_qs',num2str(round(100*qs)),'_muc',num2str(muc),'_sigmac',num2str(sigmac),'_pasn',num2str(pasn),'_pasr',num2str(pasr),'_rndseed',num2str(rndseed0),'.mat'))

toc