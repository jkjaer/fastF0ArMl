function [objValNoPitch, objValPitch] = computeObjOnGrid(Obj, dataVector)
    if strcmp(Obj.alg, 'E')
        [objValNoPitch, objValPitch] = ...
            Obj.computeArPitchFastExactObj(dataVector);
        if Obj.benchmarkAlg
            f = @() Obj.computeArPitchFastExactObj(dataVector);
            Obj.compTime = timeit(f);
        end
    elseif strcmp(Obj.alg, 'A')
        [objValNoPitch, objValPitch] = ...
            Obj.computeArPitchFastApproxObj(dataVector);
        if Obj.benchmarkAlg
            f = @() Obj.computeArPitchFastApproxObj(dataVector);
            Obj.compTime = timeit(f);
        end
    elseif strcmp(Obj.alg, 'A1')
        [objValNoPitch, objValPitch] = ...
            Obj.computeArPitchFastApprox1Obj(dataVector);
        if Obj.benchmarkAlg
            f = @() Obj.computeArPitchFastApprox1Obj(dataVector);
            Obj.compTime = timeit(f);
        end
    elseif strcmp(Obj.alg, 'A2')
        [objValNoPitch, objValPitch] = ...
            Obj.computeArPitchFastApprox2Obj(dataVector);
        if Obj.benchmarkAlg
            f = @() Obj.computeArPitchFastApprox2Obj(dataVector);
            Obj.compTime = timeit(f);
        end
    elseif strcmp(Obj.alg, 'N')
        [objValNoPitch, objValPitch] = ...
            Obj.computeArPitchExactObj(dataVector);
        if Obj.benchmarkAlg
            f = @() Obj.computeArPitchExactObj(dataVector);
            Obj.compTime = timeit(f);
        end
    end
end