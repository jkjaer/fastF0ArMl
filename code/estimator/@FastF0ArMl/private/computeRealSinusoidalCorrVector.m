function crossCorrVectors = computeRealSinusoidalCorrVector(nData, ...
        pitchGrid, pitchOrder)
    nPitches = length(pitchGrid);
    orderPitchMtx = (1:2*pitchOrder)'*(pi*pitchGrid(:)');
    crossCorrVectors = ...
        [nData*ones(1, nPitches);...
        sin(nData*orderPitchMtx)./sin(orderPitchMtx)]/2;
end
