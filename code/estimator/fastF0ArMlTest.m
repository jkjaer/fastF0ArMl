function tests = fastF0ArMlTest()
    tests = functiontests(localfunctions);
end


function testComputedObjectiveFastExact(testCase)
    rng(1)
    % setup
    nData = 199;
    pitchOrder = 5;
    fTrue = rand(1)/(2*pitchOrder);
    arOrderList = 0:4;
    nArOrders = length(arOrderList);
    for ii = 1:nArOrders
        arOrder = arOrderList(ii);
        startIndex = -(nData+arOrder-1)/2;
        timeIndices = (0:nData+arOrder-1)'+startIndex;
        pitchBounds = [1.1/nData, 0.5]; % problems with the inversion if the pitch is too low
        % generate some data (periodic signal in AR(1) noise)
        dataVector = genAr1HarmSignal(fTrue, pitchOrder, 0.8, 0.1, ...
            arOrder, nData, startIndex);

        % estimate
        EstUnderTest = FastF0ArMl(nData, pitchOrder, arOrderList(ii), ...
                    pitchBounds, 1, 'E');
        % compute the actual objective
        EstUnderTest.estimate(dataVector(1:nData));
        % compute the expected objective
        arDataMatrix = compArMtx(dataVector, arOrder);
        % for q=0 and p=0
        expObj = dataVector'*dataVector/nData;
        testCase.assertEqual(EstUnderTest.objValNoPitch(1), expObj, 'absTol', ...
            1e-12);
        % for q=0 and p>0
        for p = 1:arOrder
            expObj = (dataVector'*dataVector-dataVector'*arDataMatrix(:,1:p)*...
                (arDataMatrix(:,1:p)\dataVector))/nData;
            testCase.assertEqual(EstUnderTest.objValNoPitch(p+1), expObj, 'absTol', ...
                1e-12);
        end

        % for q>0 and any p
        for q = 1:pitchOrder
            pitchGrid = EstUnderTest.fullPitchGrid(...
                (EstUnderTest.dftRange(q,1):EstUnderTest.dftRange(q,2))+1);
            nPitches = length(pitchGrid);
            for p = 0:arOrder
                expObjectiveGrid = nan(1,nPitches);
                for kk = 1:nPitches
                    sinusoidalMatrix = ...
                        genSinMtx(2*pi*pitchGrid(kk), q, timeIndices);
                    systemMatrix = [arDataMatrix(:,1:p), sinusoidalMatrix];
                    linearParameters = systemMatrix\dataVector;
                    % real is only applied for numerical reasons
                    expObjectiveGrid(kk) = (dataVector'*dataVector-...
                        dataVector'*systemMatrix*linearParameters)/nData;
                end
                testCase.assertEqual(EstUnderTest.objValPitch(p+1,1:nPitches,q), ...
                    expObjectiveGrid, 'absTol', 1e-12);
            end
        end
        % check estimated linear parameters for estimated orders and pitch
        sinusoidalMatrix = genSinMtx(2*pi*EstUnderTest.estPitch, ...
            EstUnderTest.estPitchOrder, timeIndices);
        systemMatrix = [arDataMatrix(:,1:EstUnderTest.estArOrder), ...
            sinusoidalMatrix];
        linearParameters = systemMatrix\dataVector;
        if EstUnderTest.estArOrder > 0
            expArParams = linearParameters(1:EstUnderTest.estArOrder);
        else
            expArParams = [];
        end
        cplxSinusParams = reshape(linearParameters(...
            EstUnderTest.estArOrder+1:end),EstUnderTest.estPitchOrder ,2)*[1;-1i];
        if arOrder > 0
            cplxSinusParams = cplxSinusParams./(1-exp(-1i*2*pi*EstUnderTest.estPitch*...
                (1:EstUnderTest.estPitchOrder)'*(1:EstUnderTest.estArOrder))*expArParams);
        end
        expSinusAmps = abs(cplxSinusParams);
        expSinusPhases = angle(cplxSinusParams);
        testCase.assertEqual(EstUnderTest.estArParams, expArParams, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusAmps, expSinusAmps, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusPhases, expSinusPhases, ...
            'absTol', 1e-12);
    end
end

function testComputedObjectiveFastApprox(testCase)
    rng(1)
    % setup
    nData = 199;
    pitchOrder = 5;
    fTrue = rand(1)/(2*pitchOrder);
    arOrderList = 0:3;
    nArOrders = length(arOrderList);
    for ii = 1:nArOrders
        arOrder = arOrderList(ii);
        startIndex = -(nData+arOrder-1)/2;
        timeIndices = (0:nData+arOrder-1)'+startIndex;
        pitchBounds = [1.1/nData, 0.5]; % problems with the inversion if the pitch is too low
        % generate some data (periodic signal in AR(1) noise)
        dataVector = genAr1HarmSignal(fTrue, pitchOrder, 0.8, 0.1, ...
            arOrder, nData, startIndex);

        % estimate
        EstUnderTest = FastF0ArMl(nData, pitchOrder, arOrder, ...
                    pitchBounds, 1, 'A');
        % compute the actual objective
        EstUnderTest.estimate(dataVector(1:nData));
        % compute the expected objective
        arDataMatrix = compArMtx(dataVector, arOrder);
        % for q=0 and p=0
        expObj = dataVector'*dataVector/nData;
        testCase.assertEqual(EstUnderTest.objValNoPitch(1), expObj, 'absTol', ...
            1e-12);
        % for q=0 and p>0
        for p = 1:arOrder
            expObj = (dataVector'*dataVector-dataVector'*arDataMatrix(:,1:p)*...
                (arDataMatrix(:,1:p)\dataVector))/nData;
            testCase.assertEqual(EstUnderTest.objValNoPitch(p+1), expObj, 'absTol', ...
                1e-12);
        end

        % for q>0 and any p
        for q = 1:pitchOrder
            pitchGrid = EstUnderTest.fullPitchGrid(...
                (EstUnderTest.dftRange(q,1):EstUnderTest.dftRange(q,2))+1);
            nPitches = length(pitchGrid);
            for p = 0:arOrder
                expObjectiveGrid = nan(1,nPitches);
                for ii = 1:nPitches
                    sinusoidalMatrix = ...
                        genSinMtx(2*pi*pitchGrid(ii), q, timeIndices);
                    sinLinParams = sinusoidalMatrix\dataVector;
                    residualSin = dataVector-sinusoidalMatrix*sinLinParams;
                    if p > 0
                        ccVector = [dataVector';arDataMatrix(:,1:p)']*...
                            residualSin;
                        ZZ = toeplitz(ccVector(1:end-1));
                        arParams = inv(ZZ)*ccVector(2:end);
                        signalVarAr = ccVector(2:end)'*arParams;
                    else
                        arParams = [];
                        signalVarAr = 0;
                    end
                    % real is only applied for numerical reasons
                    expObjectiveGrid(ii) = (residualSin'*residualSin-...
                        signalVarAr)/nData;
                end
                testCase.assertEqual(EstUnderTest.objValPitch(p+1,1:nPitches,q), ...
                    expObjectiveGrid, 'absTol', 1e-12);
            end
        end
        % check estimated linear parameters for estimated orders and pitch
        sinusoidalMatrix = genSinMtx(2*pi*EstUnderTest.estPitch, ...
            EstUnderTest.estPitchOrder, timeIndices);
        systemMatrix = [arDataMatrix(:,1:EstUnderTest.estArOrder), ...
            sinusoidalMatrix];
        linearParameters = systemMatrix\dataVector;
        if EstUnderTest.estArOrder > 0
            expArParams = linearParameters(1:EstUnderTest.estArOrder);
        else
            expArParams = [];
        end
        cplxSinusParams = reshape(linearParameters(...
            EstUnderTest.estArOrder+1:end),EstUnderTest.estPitchOrder ,2)*[1;-1i];
        if arOrder > 0
            cplxSinusParams = cplxSinusParams./(1-exp(-1i*2*pi*EstUnderTest.estPitch*...
                (1:EstUnderTest.estPitchOrder)'*(1:EstUnderTest.estArOrder))*expArParams);
        end
        expSinusAmps = abs(cplxSinusParams);
        expSinusPhases = angle(cplxSinusParams);
        testCase.assertEqual(EstUnderTest.estArParams, expArParams, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusAmps, expSinusAmps, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusPhases, expSinusPhases, ...
            'absTol', 1e-12);
    end
end

function testComputedObjectiveFastApprox1(testCase)
    rng(4)
    % setup
    nData = 199;
    pitchOrder = 5;
    fTrue = rand(1)/(2*pitchOrder);
    arOrderList = 0:3;
    nArOrders = length(arOrderList);
    for ii = 1:nArOrders
        arOrder = arOrderList(ii);
        startIndex = -(nData+arOrder-1)/2;
        timeIndices = (0:nData+arOrder-1)'+startIndex;
        pitchBounds = [1.1/nData, 1]; % problems with the inversion if the pitch is too low
        % generate some data (periodic signal in AR(1) noise)

        dataVector = genAr1HarmSignal(fTrue, pitchOrder, 0.8, 1, ...
            arOrder, nData, startIndex);

        % estimate
        EstUnderTest = FastF0ArMl(nData, pitchOrder, arOrder, ...
                    pitchBounds, 1, 'A1');
        % compute the actual objective
        EstUnderTest.estimate(dataVector(1:nData));
        % compute the expected objective
        arDataMatrix = compArMtx(dataVector, arOrder);
        % for q=0 and p=0
        expObj = dataVector'*dataVector/nData;
        testCase.assertEqual(EstUnderTest.objValNoPitch(1), expObj, 'absTol', ...
            1e-12);
        % for q=0 and p>0
        for p = 1:arOrder
            expObj = (dataVector'*dataVector-dataVector'*arDataMatrix(:,1:p)*...
                (arDataMatrix(:,1:p)\dataVector))/nData;
            testCase.assertEqual(EstUnderTest.objValNoPitch(p+1), expObj, 'absTol', ...
                1e-12);
        end

        % for q>0 and any p
        for q = 1:pitchOrder
            pitchGrid = EstUnderTest.fullPitchGrid(...
                (EstUnderTest.dftRange(q,1):EstUnderTest.dftRange(q,2))+1);
            nPitches = length(pitchGrid);
            for p = 0:arOrder
                expObjectiveGrid = nan(1,nPitches);
                for ii = 1:nPitches
                    sinusoidalMatrix = ...
                        genSinMtx(2*pi*pitchGrid(ii), q, timeIndices);
                    systemMatrix = [arDataMatrix(:,1:p)'*arDataMatrix(:,1:p), ...
                        arDataMatrix(:,1:p)'*sinusoidalMatrix;...
                        sinusoidalMatrix'*arDataMatrix(:,1:p), ...
                        (nData+arOrder)*eye(2*q)/2];
                    linearParameters = inv(systemMatrix)*...
                        [arDataMatrix(:,1:p)'; sinusoidalMatrix']*dataVector;
                    % real is only applied for numerical reasons
                    expObjectiveGrid(ii) = (dataVector'*dataVector-...
                        dataVector'*[arDataMatrix(:,1:p), sinusoidalMatrix]*...
                        linearParameters)/nData;
                end
                testCase.assertEqual(EstUnderTest.objValPitch(p+1,1:nPitches,q), ...
                    expObjectiveGrid, 'absTol', 1e-12);
            end
        end
        % check estimated linear parameters for estimated orders and pitch
        sinusoidalMatrix = genSinMtx(2*pi*EstUnderTest.estPitch, ...
            EstUnderTest.estPitchOrder, timeIndices);
        systemMatrix = [arDataMatrix(:,1:EstUnderTest.estArOrder), ...
            sinusoidalMatrix];
        linearParameters = systemMatrix\dataVector;
        if EstUnderTest.estArOrder > 0
            expArParams = linearParameters(1:EstUnderTest.estArOrder);
        else
            expArParams = [];
        end
        cplxSinusParams = reshape(linearParameters(...
            EstUnderTest.estArOrder+1:end),EstUnderTest.estPitchOrder ,2)*[1;-1i];
        if arOrder > 0
            cplxSinusParams = cplxSinusParams./(1-exp(-1i*2*pi*EstUnderTest.estPitch*...
                (1:EstUnderTest.estPitchOrder)'*(1:EstUnderTest.estArOrder))*expArParams);
        end
        expSinusAmps = abs(cplxSinusParams);
        expSinusPhases = angle(cplxSinusParams);
        testCase.assertEqual(EstUnderTest.estArParams, expArParams, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusAmps, expSinusAmps, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusPhases, expSinusPhases, ...
            'absTol', 1e-12);
    end
end

function testComputedObjectiveFastApprox2(testCase)
    rng(4)
    % setup
    nData = 199;
    pitchOrder = 5;
    fTrue = rand(1)/(2*pitchOrder);
    arOrderList = 0:3;
    nArOrders = length(arOrderList);
    for ii = 1:nArOrders
        arOrder = arOrderList(ii);
        startIndex = -(nData+arOrder-1)/2;
        timeIndices = (0:nData+arOrder-1)'+startIndex;
        pitchBounds = [1.1/nData, 1]; % problems with the inversion if the pitch is too low
        % generate some data (periodic signal in AR(1) noise)

        dataVector = genAr1HarmSignal(fTrue, pitchOrder, 0.8, 1, ...
            arOrder, nData, startIndex);

        % estimate
        EstUnderTest = FastF0ArMl(nData, pitchOrder, arOrder, ...
                    pitchBounds, 1, 'A2');
        % compute the actual objective
        EstUnderTest.estimate(dataVector(1:nData));
        % compute the expected objective
        arDataMatrix = compArMtx(dataVector, arOrder);
        % for q=0 and p=0
        expObj = dataVector'*dataVector/nData;
        testCase.assertEqual(EstUnderTest.objValNoPitch(1), expObj, 'absTol', ...
            1e-12);
        % for q=0 and p>0
        for p = 1:arOrder
            expObj = (dataVector'*dataVector-dataVector'*arDataMatrix(:,1:p)*...
                (arDataMatrix(:,1:p)\dataVector))/nData;
            testCase.assertEqual(EstUnderTest.objValNoPitch(p+1), expObj, 'absTol', ...
                1e-12);
        end

        % for q>0 and any p
        for q = 1:pitchOrder
            pitchGrid = EstUnderTest.fullPitchGrid(...
                (EstUnderTest.dftRange(q,1):EstUnderTest.dftRange(q,2))+1);
            nPitches = length(pitchGrid);
            for p = 0:arOrder
                expObjectiveGrid = nan(1,nPitches);
                for ii = 1:nPitches
                    sinusoidalMatrix = ...
                        genSinMtx(2*pi*pitchGrid(ii), q, timeIndices);
                    systemMatrix = [arDataMatrix(:,1:p)'*arDataMatrix(:,1:p), ...
                        arDataMatrix(:,1:p)'*sinusoidalMatrix;...
                        sinusoidalMatrix'*arDataMatrix(:,1:p), ...
                        (nData+arOrder)*eye(2*q)/2];
                    linearParameters = inv(systemMatrix)*...
                        [arDataMatrix(:,1:p)'; sinusoidalMatrix']*dataVector;
                    % real is only applied for numerical reasons
                    expObjectiveGrid(ii) = (dataVector'*dataVector-...
                        dataVector'*[arDataMatrix(:,1:p), sinusoidalMatrix]*...
                        linearParameters)/nData;
                end
                testCase.assertEqual(EstUnderTest.objValPitch(p+1,1:nPitches,q), ...
                    expObjectiveGrid, 'absTol', 1e-12);
            end
        end
        % check estimated linear parameters for estimated orders and pitch
        sinusoidalMatrix = genSinMtx(2*pi*EstUnderTest.estPitch, ...
            EstUnderTest.estPitchOrder, timeIndices);
        systemMatrix = [arDataMatrix(:,1:EstUnderTest.estArOrder), ...
            sinusoidalMatrix];
        linearParameters = systemMatrix\dataVector;
        if EstUnderTest.estArOrder > 0
            expArParams = linearParameters(1:EstUnderTest.estArOrder);
        else
            expArParams = [];
        end
        cplxSinusParams = reshape(linearParameters(...
            EstUnderTest.estArOrder+1:end),EstUnderTest.estPitchOrder ,2)*[1;-1i];
        if arOrder > 0
            cplxSinusParams = cplxSinusParams./(1-exp(-1i*2*pi*EstUnderTest.estPitch*...
                (1:EstUnderTest.estPitchOrder)'*(1:EstUnderTest.estArOrder))*expArParams);
        end
        expSinusAmps = abs(cplxSinusParams);
        expSinusPhases = angle(cplxSinusParams);
        testCase.assertEqual(EstUnderTest.estArParams, expArParams, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusAmps, expSinusAmps, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusPhases, expSinusPhases, ...
            'absTol', 1e-12);
    end
end

function testComputedObjectiveExactNaive(testCase)
    rng(1)
    % setup
    nData = 199;
    pitchOrder = 5;
    fTrue = rand(1)/(2*pitchOrder);
    arOrderList = 0:3;
    nArOrders = length(arOrderList);
    for ii = 1:nArOrders
        arOrder = arOrderList(ii);
        startIndex = -(nData+arOrder-1)/2;
        timeIndices = (0:nData+arOrder-1)'+startIndex;
        pitchBounds = [1.1/nData, 0.5]; % problems with the inversion if the pitch is too low
        % generate some data (periodic signal in AR(1) noise)
        dataVector = genAr1HarmSignal(fTrue, pitchOrder, 0.8, 0.1, ...
            arOrder, nData, startIndex);

        % estimate
        EstUnderTest = FastF0ArMl(nData, pitchOrder, arOrder, ...
                    pitchBounds, 1, 'N');
        % compute the actual objective
        EstUnderTest.estimate(dataVector(1:nData));
        % compute the expected objective
        arDataMatrix = compArMtx(dataVector, arOrder);
        % for q=0 and p=0
        expObj = dataVector'*dataVector/nData;
        testCase.assertEqual(EstUnderTest.objValNoPitch(1), expObj, 'absTol', ...
            1e-12);
        % for q=0 and p>0
        for p = 1:arOrder
            expObj = (dataVector'*dataVector-dataVector'*arDataMatrix(:,1:p)*...
                (arDataMatrix(:,1:p)\dataVector))/nData;
            testCase.assertEqual(EstUnderTest.objValNoPitch(p+1), expObj, 'absTol', ...
                1e-12);
        end

        % for q>0 and any p
        for q = 1:pitchOrder
            pitchGrid = EstUnderTest.fullPitchGrid(...
                (EstUnderTest.dftRange(q,1):EstUnderTest.dftRange(q,2))+1);
            nPitches = length(pitchGrid);
            for p = 0:arOrder
                expObjectiveGrid = nan(1,nPitches);
                for ii = 1:nPitches
                    sinusoidalMatrix = ...
                        genSinMtx(2*pi*pitchGrid(ii), q, timeIndices);
                    systemMatrix = [arDataMatrix(:,1:p), sinusoidalMatrix];
                    linearParameters = systemMatrix\dataVector;
                    % real is only applied for numerical reasons
                    expObjectiveGrid(ii) = (dataVector'*dataVector-...
                        dataVector'*systemMatrix*linearParameters)/nData;
                end
                testCase.assertEqual(EstUnderTest.objValPitch(p+1,1:nPitches,q), ...
                    expObjectiveGrid, 'absTol', 1e-12);
            end
        end
        % check estimated linear parameters for estimated orders and pitch
        sinusoidalMatrix = genSinMtx(2*pi*EstUnderTest.estPitch, ...
            EstUnderTest.estPitchOrder, timeIndices);
        systemMatrix = [arDataMatrix(:,1:EstUnderTest.estArOrder), ...
            sinusoidalMatrix];
        linearParameters = systemMatrix\dataVector;
        if EstUnderTest.estArOrder > 0
            expArParams = linearParameters(1:EstUnderTest.estArOrder);
        else
            expArParams = [];
        end
        cplxSinusParams = reshape(linearParameters(...
            EstUnderTest.estArOrder+1:end),EstUnderTest.estPitchOrder ,2)*[1;-1i];
        if arOrder > 0
            cplxSinusParams = cplxSinusParams./(1-exp(-1i*2*pi*EstUnderTest.estPitch*...
                (1:EstUnderTest.estPitchOrder)'*(1:EstUnderTest.estArOrder))*expArParams);
        end
        expSinusAmps = abs(cplxSinusParams);
        expSinusPhases = angle(cplxSinusParams);
        testCase.assertEqual(EstUnderTest.estArParams, expArParams, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusAmps, expSinusAmps, ...
            'absTol', 1e-12);
        testCase.assertEqual(EstUnderTest.estSinusPhases, expSinusPhases, ...
            'absTol', 1e-12);
    end
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
        arOrder, nData, startIndex)
    timeIndices = (0:nData-1)'+startIndex;
    [~, cSinMtx] = genSinMtx(2*pi*fTrue, pitchOrder, timeIndices);
    periodicSignal = real(cSinMtx*...
        (ones(pitchOrder,1).*exp(1i*2*pi*rand(pitchOrder,1))));
    arNoise = sqrt(exVar)*filter(1,[1,-arParam],randn(nData,1));
    signal = [periodicSignal+arNoise; zeros(arOrder,1)];
end
