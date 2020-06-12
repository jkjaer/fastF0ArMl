function [dataVector, truePitch, sinusAmps, arParams, exVar] = ...
        generateData(nData, snrDb, pitchOrder, arOrder, pitchBounds)
    % gererate periodic signal
%     minPitch = 2.5/nData;
    minPitch = pitchBounds(1);
    maxPitch = pitchBounds(2)/pitchOrder;
    truePitch = rand(1)*(maxPitch-minPitch)+minPitch;
    sinusAmps = ones(pitchOrder,1).*2.^(-(0:pitchOrder-1)');
    sinusPhases = 2*pi*rand(pitchOrder,1);
    cisMtx = exp(1i*2*pi*truePitch*(0:nData-1)'*(1:pitchOrder));
    perSignal = real(cisMtx*(sinusAmps.*exp(1i*sinusPhases)));
    % generate AR-parameters
    rootAbsMin = 0.5;
    rootAbsMax = 0.90;
    arParams = genRndArParams(arOrder, rootAbsMin, rootAbsMax);
    arVect = [1; -arParams];
    % pre-filter the harmonic signal
    perSignal = filter(arVect, 1, perSignal);
    % compute the value of the excitation variance
    perPower = perSignal'*perSignal/nData;
    unitWgn = randn(nData,1);
    exVar = 10^(-snrDb/10)*perPower/(unitWgn'*unitWgn/nData);
    noisySignal = perSignal+sqrt(exVar)*unitWgn;
    % post-filter the signal to get the data vector with AR-noise
    dataVector = filter(1,arVect,noisySignal);
end