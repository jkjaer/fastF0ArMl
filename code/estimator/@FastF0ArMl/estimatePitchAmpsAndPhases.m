function [estSinusAmps, estSinusPhases] = ...
        estimatePitchAmpsAndPhases(Obj, estArParams, ...
        shapedPitchLinearParameters, estPitch)
    arOrder = length(estArParams);
    pitchOrder = length(shapedPitchLinearParameters)/2;
    shapedPitchComplexLinearParameters = reshape(...
        shapedPitchLinearParameters', 2, pitchOrder)'*[1;-1i];
    % remove the contribution of the AR-parameters to the pitch
    % parameters
    if arOrder == 0
        pitchComplexLinearParameters = ...
            shapedPitchComplexLinearParameters;
    else
        pitchComplexLinearParameters = ...
            shapedPitchComplexLinearParameters./(1-...
            exp(-1i*(2*pi*estPitch/Obj.samplingFreq)*...
            (1:pitchOrder)'*(1:arOrder))*estArParams);
    end
    estSinusPhases = angle(pitchComplexLinearParameters);
    estSinusAmps = abs(pitchComplexLinearParameters);
end