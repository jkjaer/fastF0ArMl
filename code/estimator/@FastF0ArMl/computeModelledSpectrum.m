function [modelledSpectrum, arPsd, pitchPsd, freqVector] = ...
        computeModelledSpectrum(Obj, nDft)
    if isempty(Obj.objValNoPitch)
        error('Please run the estimator (Obj.estimate) first!'); 
    end
    if nargin < 2
        nDft = Obj.nDft;
    end
    % AR spectrum
    arPsd = Obj.estExVar./(Obj.samplingFreq*...
        abs(fft([1;-Obj.estArParams], nDft,1)).^2);
    % pitch spectrum
    pitchPsd = zeros(nDft,1);
    pitchDftIdx = ...
        round(Obj.estPitch*nDft/Obj.samplingFreq);
    pitchPsd(1+pitchDftIdx*(1:Obj.estPitchOrder)) = ...
        length(Obj.dataVector)*Obj.estSinusAmps.^2/...
        (Obj.samplingFreq);
    pitchPsd(nDft-1-pitchDftIdx*(1:Obj.estPitchOrder)) = ...
        length(Obj.dataVector)*Obj.estSinusAmps.^2/...
        (Obj.samplingFreq);
    pitchPsd = pitchPsd/4;
    % total spectrum
    modelledSpectrum = pitchPsd+arPsd;
    freqVector = Obj.samplingFreq*(0:nDft-1)'/nDft;
end
