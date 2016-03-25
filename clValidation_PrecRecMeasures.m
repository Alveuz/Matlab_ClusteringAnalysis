function [prec, rec] = clValidation_PrecRecMeasures( vecAssign, vecObjective )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

vecObjctSize = size(vecObjective,1);

if(vecObjctSize==1)

    vecX1Idxes = find(vecAssign);
    vecX2Idxes = find(vecObjective);

    vecIntersection = intersect(vecX1Idxes, vecX2Idxes);
    vecInterSize    = length(vecIntersection);

    vecX1Size = sum(vecAssign);
    vecX2Size  = sum(vecObjective);

    prec = vecInterSize/vecX1Size;
    rec  = vecInterSize/vecX2Size;
else
    vecX1Idxes  = find(vecAssign);
    vecX1Size   = sum(vecAssign);
    
    prec        = zeros(1,vecObjctSize);
    rec         = zeros(1,vecObjctSize);
    
    for i=1:vecObjctSize
        
        vecX2Idxes      = find(vecObjective(i,:));
        vecX2Size       = sum(vecObjective(i,:));
        
        vecIntersection = intersect(vecX1Idxes, vecX2Idxes);
        vecInterSize    = length(vecIntersection);
        
        prec(1,i)   = vecInterSize/vecX1Size;
        rec(1,i)    = vecInterSize/vecX2Size;
    end
end

end

