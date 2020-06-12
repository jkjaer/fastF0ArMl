function tests = iterArF0Test()
    addpath('../../estimator');
    tests = functiontests(localfunctions);
end


function testComputedEstimate(testCase)
    rng(1)
    % setup
    nData = 199;
    pitchOrder = 5;
    fTrue = rand(1)/(2*pitchOrder);
    maxArOrder = 3;
    startIndex = -nData/2;
    pitchBounds = [1.1/(nData), 0.5]; % problems with the inversion if the pitch is too low
    % generate some data (periodic signal in AR(1) noise)
    dataVector = genAr1HarmSignal(fTrue, pitchOrder, 0.8, 0.1, ...
        nData, startIndex);

    % estimate
    EstUnderTest =IterArF0(nData, pitchOrder, maxArOrder, ...
                pitchBounds, 1, 'E');
    % compute the actual objective
    estPitch = EstUnderTest.estimate(dataVector, 100);
    %% really hard to test
end

function [sinMtx, cSinMtx] = genSinMtx(om, order, timeIndices)
    cSinMtx = exp(1i*om*timeIndices*(1:order));
    sinMtx = [real(cSinMtx), imag(cSinMtx)];
end

function arMtx = compArMtx(dataVector, arOrder)
    if arOrder > 0
        arMtx = toeplitz([0; dataVector(1:end-1)],zeros(1,arOrder));
    else
        arMtx = zeros(length(dataVector),0);
    end
end

function signal = genAr1HarmSignal(fTrue, pitchOrder, arParam, exVar, ...
        nData, startIndex)
    timeIndices = (0:nData-1)'+startIndex;
    [~, cSinMtx] = genSinMtx(2*pi*fTrue, pitchOrder, timeIndices);
    periodicSignal = real(cSinMtx*...
        (ones(pitchOrder,1).*exp(1i*2*pi*rand(pitchOrder,1))));
    arNoise = sqrt(exVar)*filter(1,[1,-arParam],randn(nData,1));
    signal = periodicSignal+arNoise;
end
