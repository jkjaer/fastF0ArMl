function [fullPitchGrid, dftRange] = ...
        computePitchGrid(nDft, pitchBounds, maxPitchOrder, samplingFreq)
    fullPitchGrid = samplingFreq*(0:nDft-1)'/nDft;
    minDftIdx = max(0, round(pitchBounds(1)*nDft/samplingFreq));
    maxDftIdx = round(min(nDft/2-1, pitchBounds(2)*nDft/samplingFreq));
    dftRange = nan(maxPitchOrder, 2);
    dftRange(1,:) = [minDftIdx,maxDftIdx];
    for ii = 2:maxPitchOrder
        % the max pitch must be smaller than fs/(2*modelOrder) to avoid 
        % problems with aliasing
        maxDftIdx = round(min(nDft/(2*ii)-1, ...
            pitchBounds(2)*nDft/samplingFreq));
        dftRange(ii,:) = [minDftIdx, maxDftIdx];
    end
end
