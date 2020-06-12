function tests = cholDowndatingTest()
    tests = functiontests(localfunctions);
end

% Validate the produced indices
function testAlg(testCase)
    rng(1);
    nData = 32;
    nRows = 5;
    nCases = 4;
    newVect = randn(nRows,nCases);
    cholVectInit = nan(nRows*(nRows+1)/2,nCases);
    expCholVectPost = cholVectInit;
    for ii = 1:nCases
        tmp = randn(nRows,nData);
        psMtxPrior = tmp*tmp';
        cholFactPrior = chol(psMtxPrior)';
        psMtxPost = psMtxPrior-newVect(:,ii)*newVect(:,ii)';
        cholFactPost = chol(psMtxPost)';
        cholIdx = 1:nRows;
        for jj = 1:nRows
            cholVectInit(cholIdx,ii) = cholFactPrior(jj:end,jj);
            expCholVectPost(cholIdx,ii) = cholFactPost(jj:end,jj);
            cholIdx = cholIdx(1:end-1) + (nRows-jj+1);
        end
    end
    actCholVectPost = cholDowndating(cholVectInit, newVect);
    testCase.assertEqual(actCholVectPost, expCholVectPost, ...
        'absTol', 1e-14);
end