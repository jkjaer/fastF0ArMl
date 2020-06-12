function tests = computeObjectiveNaivelyTest()
    tests = functiontests(localfunctions);
end

function testComputedObjectiveReal(testCase)
    % setup
    nData = 199;
    arOrder = 3;
    startIndex = -(nData+arOrder-1)/2;
    pitchOrder = 5;
    nDft = 5*(nData+arOrder)*pitchOrder;
    pitchBounds = [0.5/nData, 1]; % problems with the inversion if the pitch is too low
    arMethod = 'circulant';
    dataVector = [randn(nData,1); zeros(arOrder,1)];
    % compute the expected objective
    [fullPitchGrid, dftRange] = computePitchGrid(nDft, pitchBounds, ...
        pitchOrder, 1);
    pitchGrid = fullPitchGrid(dftRange(pitchOrder,1):dftRange(pitchOrder,2));
    nPitches = length(pitchGrid);
    arDataMatrix = computeArDataMatrix(dataVector, arOrder, arMethod);
    expObjectiveGrid = nan(nPitches,1);
    for ii = 1:nPitches
        sinusoidalMatrix = computeSinusoidalMatrix(nData+arOrder, startIndex, ...
            2*pi*pitchGrid(ii)*(1:pitchOrder), true);
        systemMatrix = [arDataMatrix, sinusoidalMatrix];
        linearParameters = systemMatrix\dataVector;
        % real is only applied for numerical reasons
        expObjectiveGrid(ii) = (dataVector'*dataVector-...
            dataVector'*systemMatrix*linearParameters)/nData;
    end
    % compute the actual objective
    actObjectiveGrid = computeObjectiveNaively(dataVector, ...
        pitchGrid, pitchOrder, arOrder, startIndex, nData);
    testCase.assertEqual(actObjectiveGrid, expObjectiveGrid, 'absTol', ...
        1e-12);
end
