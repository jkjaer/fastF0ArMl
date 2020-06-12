function [tachoSignal, time] = ...
        extractTachoSignalFromWav(fileName, extension)
    [rawData, samplingFreq] = audioread([fileName, '.', extension]);
    info = audioinfo([fileName, '.', extension]);
    nBits = info.BitsPerSample;
    [nData, nChannels] = size(rawData);
    % the tacho signal is the LSB of the audio signal in each channel
    time = (0:nData-1)'/samplingFreq;
    tachoSignal = nan(nData,nChannels);
    for iChannel = 1:nChannels
       tachoSignal(:,iChannel) = mod(round(2^(nBits-1)*rawData(:,iChannel)),2);
    end
end

