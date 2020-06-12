function tests = orderRecursiveForwardSubstitutionTest()
    tests = functiontests(localfunctions);
end

% Validate the produced indices
function testAlg(testCase)
    rng(1);
    nData = 32;
    nRows = 5;
    nCases = 4;
    cholVct = nan(nRows*(nRows+1)/2,nCases);
    rightHandSides = randn(nRows, nCases);
    expResVar = nan(nRows,nCases);
    for ii = 1:nCases
        tmp = randn(nRows,nData);
        psMtxPrior = tmp*tmp';
        cholFact = chol(psMtxPrior)';
        cholIdx = 1:nRows;
        for jj = 1:nRows
            cholVct(cholIdx,ii) = cholFact(jj:end,jj);
            cholIdx = cholIdx(1:end-1) + (nRows-jj+1);
            expResVar(jj,ii) = rightHandSides(1:jj,ii)'*...
                (psMtxPrior(1:jj,1:jj)\rightHandSides(1:jj,ii));
        end
    end
    actResVar = orderRecursiveForwardSubstitution(cholVct, rightHandSides);
    testCase.assertEqual(actResVar, expResVar, 'absTol', 1e-15);
end