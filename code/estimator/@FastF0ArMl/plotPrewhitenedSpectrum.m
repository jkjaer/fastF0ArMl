function plotPrewhitenedSpectrum(Obj, plotRangeDb)
    [preWhitenedSpectrum, freqVector] = ...
            Obj.computePrewhitenedSpectrum;
    perEst = abs(fft(Obj.dataVector, Obj.nDft)).^2/...
        (Obj.samplingFreq*length(Obj.dataVector));
    psdsInDb = 10*log10([perEst, preWhitenedSpectrum]);
    plot(freqVector, psdsInDb, 'lineWidth', 2);
    xlim([0,Obj.samplingFreq/2])
    if nargin > 1
        ylim([-plotRangeDb, 0]+max(max(psdsInDb)));
    end
    xlabel('Pitch [Hz]');
    ylabel('PSD [dB/Hz]')
    legend('Periodogram', 'Prewhitened periodogram');
end
