function tests = genSchurToeplitzAlgTest()
    tests = functiontests(localfunctions);
end

% Validate the produced indices
function testAlg(testCase)
    nData = 1000;
    maxArOrder = 25;
    dataVector = [randn(nData,1);zeros(maxArOrder,1)];
    arMtx = toeplitz([0;dataVector(1:end-1)],zeros(1,maxArOrder));
    covVect = [dataVector, arMtx]'*dataVector/(nData+maxArOrder);
    TplzMtx = toeplitz(covVect(1:end-1));
    cholFact = chol(TplzMtx)';
    expCholVect = nan(maxArOrder*(maxArOrder+1)/2,1);
    expPredErrVar = nan(maxArOrder+1,1);
    expPredErrVar(1) = covVect(1);
    cholIdx = 1:maxArOrder;
    for p = 1:maxArOrder
        expCholVect(cholIdx) = cholFact(p:end,p);
        [~, expPredErrVar(p+1)] = levinson(covVect(1:p+1));
        cholIdx = cholIdx(1:end-1) + (maxArOrder-p+1);
    end
    [actPredErrVar, actCholVect] = genSchurToeplitzAlg(covVect);
    testCase.assertEqual(actPredErrVar, expPredErrVar, 'absTol', 1e-15);
    testCase.assertEqual(actCholVect, expCholVect, 'absTol', 1e-15);
end