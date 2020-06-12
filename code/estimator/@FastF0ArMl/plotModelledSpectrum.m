function plotModelledSpectrum(Obj, plotRangeDb)
    [modelledSpectrum, arPsd, ~, freqVector] = ...
        Obj.computeModelledSpectrum;
    perEst = abs(fft(Obj.dataVector, Obj.nDft)).^2/...
        (Obj.samplingFreq*length(Obj.dataVector));
    psdsInDb = 10*log10([perEst, modelledSpectrum, arPsd]);
    plot(freqVector, psdsInDb, 'lineWidth', 2);
    xlim([0,Obj.samplingFreq/2])
    if nargin > 1
        ylim([-plotRangeDb, 0]+max(max(psdsInDb)));
    end
    xlabel('Frequency [Hz]');
    ylabel('PSD [dB/Hz]')
    legend('Periodogram', 'Harmonic+AR spectrum', 'AR PSD');
end