function r = FamilyEncounter(familyMat, resourceMat, qA, qS)
    sz = size(familyMat);
    r = zeros(sz);
    szdim1 = size(familyMat,1);
    szdim2 = size(familyMat,2);
    for i = 1:szdim1
        for j = 1:szdim2
            if (familyMat(i,j) == 1) && (resourceMat(i,j) == 1)
                r(i,j) = 1+qA;
            elseif (familyMat(i,j) == 0) && (resourceMat(i,j) == 1)
                r(i,j) = 1-qA;
            elseif (familyMat(i,j) == 0) && (resourceMat(i,j) == 0)
                r(i,j) = 1+qS;
            elseif (familyMat(i,j) == 1) && (resourceMat(i,j) == 0)
                r(i,j) = 1-qS;
            end
        end
    end
return;