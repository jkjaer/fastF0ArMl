function [estPitch, estSinusAmps, estSinusPhases, estArParams, ...
        estExVar] = estimatePitchArParameters(Obj, dataVector, ...
        pitchOrder, arOrder, refinementTol)
    if pitchOrder == 0
        estPitch = nan;
    else
        [~, pitchIdx] = ...
            min(Obj.objValPitch(arOrder+1,:, pitchOrder));
        coarsePitchEst = ...
            Obj.fullPitchGrid(Obj.dftRange(pitchOrder,1)+pitchIdx);
        % refine the pitch estimate if needed
        if refinementTol < Obj.samplingFreq/Obj.nDft
            % add refinement here if neceessary!!!
            pitchLimits = coarsePitchEst/Obj.samplingFreq+[-1, 1]/Obj.nDft;
            costFunction = @(f0) computeObjectiveNaively(dataVector, ...
                f0, pitchOrder, arOrder, Obj.startIndex, Obj.nData);
            [f0l, f0u] = fibonacciSearch(costFunction, ...
                pitchLimits(1), pitchLimits(2), refinementTol);
            estPitch = (f0l+f0u)*Obj.samplingFreq/2;
        else
            estPitch = coarsePitchEst;
        end
    end
    [estExVar, estLinParams] = ...
        computeObjectiveNaively(dataVector, estPitch/Obj.samplingFreq, ...
        pitchOrder, arOrder, Obj.startIndex, Obj.nData);
    if arOrder == 0
        estArParams = zeros(0);
    else
        estArParams = estLinParams(1:arOrder);
    end
    if pitchOrder == 0
        estSinusAmps = zeros(0);
        estSinusPhases = zeros(0);
    else
        [estSinusAmps, estSinusPhases] = ...
            Obj.estimatePitchAmpsAndPhases(estArParams, ...
            estLinParams(arOrder+1:end), estPitch);
    end
end