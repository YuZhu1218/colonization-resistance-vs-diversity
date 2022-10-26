function [Ne, Cmp, Ne0, Cmp0] = WellmixedInteraction(nRound,cellRatioArray,nInitialCell,cia,l,m,R0Array,dab,taust,tauf,dtau,DilTh,ExtTh)

    tau0 = 0;
    [nCellType, nMediator] = size(cia);
    
    %% Initial state 
    cMedVector = R0Array; % concentrations of interaction mediators
    
    %% Cell-growth time-course
    tauScaleArray = tau0:dtau:tauf;
    nTauScale = length(tauScaleArray);

    tc = zeros(1,nRound*nTauScale);
    Xc = zeros(nCellType,nRound*nTauScale);
    Cc = zeros(nMediator,nRound*nTauScale);

    cct = 0;
    nCellVector = nInitialCell * cellRatioArray'; % initial number of each cell type
    
    for iRound = 1 : nRound
        %%%%%%%%%
        cMedVector = nInitialCell / sum(nCellVector) * cMedVector; 
        %%%%%%%%%
        nCellVector = nInitialCell * cellRatioArray'; % initial number of each cell type

        tau0 = 0; % in hours
        tau = tau0;

        nCellOnEachScale = zeros(nCellType,nTauScale);
        cMedOnEachScale = zeros(nMediator,nTauScale);
        rzs = zeros(nCellType,nTauScale);
        rzm = zeros(nMediator,nTauScale);
        count = 0;
        while (tau<=tauf-dtau) && (sum(nCellOnEachScale(:,max(count,1)))<DilTh)

            count = count+1;
            tau = tauScaleArray(count);

            rIntPerCellVector = cia * cMedVector - m;
            
            Re = cMedVector * ones(1,nCellType);
            drdt = (R0Array - cMedVector)./taust - (nCellVector'*(Re'.*cia))' + ((nCellVector'*cia)*dab*l)';
 
            nCellVector = nCellVector + dtau * (rIntPerCellVector .* nCellVector);
            nCellVector(nCellVector < ExtTh) = 0;  

            cMedVector = cMedVector + dtau*drdt;
            cMedVector(cMedVector<0) = 0;

            nCellOnEachScale(:,count) = nCellVector;
            cMedOnEachScale(:,count) = cMedVector;
            rzs(:,count) = rIntPerCellVector;
            rzm(:,count) = drdt;

        end
        cellRatioArray = 1/sum(nCellOnEachScale(:,count))*nCellOnEachScale(:,count)';

        if cct ==0
            tc(cct+1:cct+count) = tauScaleArray(1:count);
        else
            tc(cct+1:cct+count) = tc(cct) + tauScaleArray(1:count);
        end
        Xc(:,cct+1:cct+count) = nCellOnEachScale(:,1:count);
        Cc(:,cct+1:cct+count) = cMedOnEachScale(:,1:count);
        cct = cct+count;

    end
    indx = 1:nCellType;
    Ne0 = indx(nCellOnEachScale(:,count)>0);
    Cmp0 = cellRatioArray(Ne0);
    % get Cmp as percentage each cell type contributes to the total community
    if sum(Cmp0) > 0
        Cmp_sum = zeros(1,size(Cmp0,2));
        Cmp_sum(1,:) = sum(Cmp0);

        Cmp0 = Cmp0./Cmp_sum;
    end

    nGen = log(DilTh/nInitialCell)/log(2);
    r = (nCellOnEachScale(:,count)./nCellOnEachScale(:,1)).^(20/nGen);
    stp = (r > abs(0.9*max(r)));
    Ne = indx(stp);
    Cmp = cellRatioArray(Ne);
    % get Cmp as percentage each cell type contributes to the total community
    if sum(Cmp) > 0
        Cmp_sum = zeros(1,size(Cmp,2));
        Cmp_sum(1,:) = sum(Cmp);

        Cmp = Cmp./Cmp_sum;
    end

return;