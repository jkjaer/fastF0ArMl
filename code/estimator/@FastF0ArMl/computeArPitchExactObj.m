 function [objValNoPitch, objValPitch] = ...
        computeArPitchExactObj(Obj, dataVector)
    for ii = 0:Obj.maxPitchOrder
        if ii == 0
            iiPitchGrid = nan;
            objValNoPitch = nan(Obj.maxArOrder+1,1);
        else
            % the +1 is to compensate for MATLAB's indexing
            iiPitchGrid = ...
                Obj.fullPitchGrid((Obj.dftRange(ii,1):Obj.dftRange(ii,2))+1);
            nPitches = length(iiPitchGrid);
            if ii == 1
                objValPitch = ...
                    nan(Obj.maxArOrder+1, nPitches, Obj.maxPitchOrder);
            end
        end
        
        for jj = 0:Obj.maxArOrder
            if ii == 0
                objValNoPitch(jj+1) = ...
                    computeObjectiveNaively(dataVector, iiPitchGrid, ...
                    ii, jj, Obj.startIndex, Obj.nData);
            else
                objValPitch(jj+1,1:nPitches,ii) = ...
                    computeObjectiveNaively(dataVector, iiPitchGrid, ...
                    ii, jj, Obj.startIndex, Obj.nData).';
            end
        end
    end
end
