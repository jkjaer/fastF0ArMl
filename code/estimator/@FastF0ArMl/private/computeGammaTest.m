function tests = computeGammaTest()
    tests = functiontests(localfunctions);
end

function testComputedGamma(testCase)
    % setup
    rng(10)
    nData = 100;
    maxArOrder = 3;
    maxPitchOrder = 5;
    nDft = 5*maxPitchOrder*nData;
    pitchBounds = [1, nData/2-1]/nData;
    [fullPitchGrid, dftRange] = ...
    	computePitchGrid(nDft, pitchBounds, maxPitchOrder, 1);
    % compute the expected objective for a pitch
    expGammaTpH = cell(maxPitchOrder,1);
    expGammaTmH = cell(maxPitchOrder,1);
    for q = 1:maxPitchOrder
        [fullGammaTpH, fullGammaTmH] = ...
            computeGammasNaively(nData+maxArOrder, ...
            fullPitchGrid((dftRange(q,1):dftRange(q,2))+1), q);
        expGammaTpH{q} = fullGammaTpH;
        expGammaTmH{q} = fullGammaTmH;
    end
%   compute the actual objectives
    sinusCrossCorrVector = ...
        computeRealSinusoidalCorrVector(nData+maxArOrder, ...
        fullPitchGrid((dftRange(1,1):dftRange(1,2))+1), maxPitchOrder);
    [actGammaTpH, actGammaTmH] = computeGamma(maxPitchOrder, ...
        dftRange, sinusCrossCorrVector);
%   test produces results
    testCase.assertEqual(actGammaTpH, expGammaTpH, 'absTol', 1e-12);
    testCase.assertEqual(actGammaTmH, expGammaTmH, 'absTol', 1e-12);
end