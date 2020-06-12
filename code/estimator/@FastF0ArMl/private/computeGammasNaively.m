function [gammaTpH, gammaTmH] = computeGammasNaively(nData, pitchGrid, ...
        pitchOrder)
    nPitches = length(pitchGrid);
    revUnitVectorTmH = [zeros(pitchOrder-1,1); 1];
    gammaTmH = nan(pitchOrder,nPitches);
    revUnitVectorTpH = revUnitVectorTmH;
    gammaTpH = gammaTmH;
    timeIndices = (0:nData-1)'-(nData-1)/2;
    for ii = 1:nPitches
        tmp = exp(1i*timeIndices*2*pi*pitchGrid(ii)*(1:pitchOrder));
        sinusMtx = [real(tmp), imag(tmp)];
        sinusMtxTpH = sinusMtx(:,1:pitchOrder);
        sinusMtxTmH = sinusMtx(:,pitchOrder+1:end);
        gammaTpH(:,ii) = (sinusMtxTpH'*sinusMtxTpH)\revUnitVectorTpH;
        gammaTmH(:,ii) = (sinusMtxTmH'*sinusMtxTmH)\revUnitVectorTmH;
        gammaTpH(:,ii) = gammaTpH(:,ii)/sqrt(gammaTpH(end,ii));
        gammaTmH(:,ii) = gammaTmH(:,ii)/sqrt(gammaTmH(end,ii));
    end
end