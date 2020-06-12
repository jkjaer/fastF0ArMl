function [preWhitenedSpectrum, freqVector] = ...
        computePrewhitenedSpectrum(Obj, nDft)
    if nargin < 2
        nDft = Obj.nDft;
    end
    freqVector = Obj.samplingFreq*(0:nDft-1)'/nDft;
    prewhitenedData = ifft(fft(Obj.dataVector, nDft, 1).*...
        fft([1;-Obj.estArParams], nDft, 1));
    preWhitenedSpectrum = abs(fft(prewhitenedData, nDft)).^2/...
        (Obj.samplingFreq*length(Obj.dataVector));
end