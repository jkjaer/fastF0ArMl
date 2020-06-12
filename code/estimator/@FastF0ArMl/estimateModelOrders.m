function [estArOrder, estPitchOrder] = estimateModelOrders(Obj)
    Obj.bayesFactor = pitchArBicModelComparison(Obj.nData, ...
        Obj.objValNoPitch, Obj.objValPitch, Obj.validPitchOrders);
%     Obj.bayesFactor = pitchArModelComparison(Obj.nData, Obj.objValNoPitch, ...
%         Obj.objValPitch, Obj.dftRange, Obj.nDft, 3, Obj.validPitchOrders);
    [idxPitch, idxAr] = ...
        find(Obj.bayesFactor == max(max(Obj.bayesFactor)));
    estPitchOrder = idxPitch-1;
    estArOrder = idxAr-1;
end