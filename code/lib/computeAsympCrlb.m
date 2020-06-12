function aCrlb = computeAsympCrlb(nData, freqs, sinusAmps, arParams, ...
        exVar, samplingFreq)
    if nargin < 6
        samplingFreq = 1;
    end
    pitchOrder = length(sinusAmps);
    arOrder = length(arParams);
    arPsdVals = ...
        exVar./abs(1-exp(-1i*2*pi*freqs(:)*(1:arOrder)/samplingFreq)*...
        arParams(:)).^2;
    aCrlb = samplingFreq^2*24/((2*pi)^2*nData^3*sum(((1:pitchOrder)').^2.*...
        sinusAmps(:).^2./arPsdVals));
end