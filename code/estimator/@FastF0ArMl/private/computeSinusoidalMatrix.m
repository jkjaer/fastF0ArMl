function sinusoidalMatrix = computeSinusoidalMatrix(nData, startIndex, ...
        freqList, dataAreRealValued)
    if nargin < 4
        dataAreRealValued = false;
    end
    timeIndices = (0:nData-1)'+startIndex;
    sinusoidalMatrix = exp(1i*timeIndices*freqList);
    if dataAreRealValued
        sinusoidalMatrix = kron(real(sinusoidalMatrix),[1,0]) +...
            kron(imag(sinusoidalMatrix),[0,1]);
    end
end

