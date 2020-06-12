function [objectiveGrid, estLinParams] = ...
        computeObjectiveNaively(dataVector, pitchGrid, pitchOrder, ...
        arOrder, startIndex, nData)
    dataPower = dataVector'*dataVector/nData;
    maxArOrder = length(dataVector)-nData;
    timeIndices = (0:nData+maxArOrder-1)'+startIndex;
    if arOrder == 0
        arDataMatrix = zeros(nData+arOrder,0);
    else
        arDataMatrix = toeplitz([0; dataVector(1:end-1)], zeros(1,arOrder));
    end
    if any(isnan(pitchGrid)) && size(arDataMatrix,2) == 0
        % no pitch and AR part
        objectiveGrid = dataPower;
        estLinParams = nan;
    elseif any(isnan(pitchGrid)) % no pitch
        estLinParams = arDataMatrix\dataVector;
        objectiveGrid = dataPower-...
            dataVector'*arDataMatrix*estLinParams/nData;
    else
        nPitches = length(pitchGrid);
        objectiveGrid = nan(nPitches,1);
        estLinParams = nan(arOrder+2*pitchOrder,nPitches);
        for ii = 1:nPitches
            cplxMtx = exp(1i*2*pi*timeIndices*pitchGrid(ii)*(1:pitchOrder));
            sinusoidalMatrix = ...
                kron(real(cplxMtx),[1,0])+kron(imag(cplxMtx),[0,1]);
            % compute the linear parameters
            systemMtx = [arDataMatrix, sinusoidalMatrix];
            estLinParams(:,ii) = systemMtx\dataVector;
            % compute the objective
            objectiveGrid(ii) = dataPower-dataVector'*systemMtx*...
                estLinParams(:,ii)/nData;
        end
    end
end
