function tests = fibonacciSearchTest
    tests = functiontests(localfunctions);
end

% Validate that the returned interval of uncertainty satisfies requirements
function testReturnedIntervalOfUncertainty(testCase)
    f = @(x) x.^2-2*x-3;
    fMin = 1;
    requiredUncertaintyArray = 10.^(-(0:7));
    nRequiredUncertaintyArray = length(requiredUncertaintyArray);
    lowerBound = -2;
    upperBound = 5;
    for iRequiredUncertainty = 1:nRequiredUncertaintyArray
        [lowerFinalBound, upperFinalBound] = fibonacciSearch(f,lowerBound,...
            upperBound,requiredUncertaintyArray(iRequiredUncertainty));
        fMinIsInFinalInterval = lowerFinalBound<fMin && fMin<upperFinalBound;
        testCase.assertTrue(fMinIsInFinalInterval);
        finalIntervalIsSmallEnough = (upperFinalBound-lowerFinalBound)<=...
            requiredUncertaintyArray(iRequiredUncertainty);
        testCase.assertTrue(finalIntervalIsSmallEnough);
    end
end

% Validate number of inputs to the contructer
function testConstructorInputs(testCase)
    f = @(x) x.^2-2*x-3;
    lowerBound = -2;
    upperBound = 5;
    requiredUncertainty = 1e-3;
    inputArray = {f,lowerBound,upperBound,requiredUncertainty};
    nInputArray = length(inputArray);
    validInputNoArray = [4];
    for iInput = 1:nInputArray+1
        nInputs = iInput-1;
        if ismember(nInputs,validInputNoArray)
            [lowerFinalBound, upperFinalBound] = fibonacciSearch(...
                inputArray{1:nInputs});
            testCase.assertGreaterThan(upperFinalBound-lowerFinalBound,0);
        else
            testCase.assertError(@()fibonacciSearch(inputArray{1:nInputs}),...
                'fibonacciSearch:argChk');
        end
    end
end

% Input validation for setFunction
function testInputsInSetFunction(testCase)
    validFuncArray = {@(x) 2*x-3, @(h) ones(1,10)*h};
    nValidFuncArray = length(validFuncArray);
    invalidFuncArray = {nan, 1i, 47, inf, -inf,''};
    funcArray = [validFuncArray,invalidFuncArray];
    nFuncArray = length(funcArray);
    lowerBound = -2;
    upperBound = 5;
    requiredUncertainty = 1e-3;
    for iFunc = 1:nFuncArray
        if iFunc <= nValidFuncArray
            [lowerFinalBound, upperFinalBound] = fibonacciSearch(...
                funcArray{iFunc},lowerBound,upperBound,requiredUncertainty);
            testCase.assertGreaterThan(upperFinalBound-lowerFinalBound,0);
        else
            testCase.assertError(@()fibonacciSearch(...
                funcArray{iFunc},lowerBound,upperBound,requiredUncertainty),...
                'fibonacciSearch:argChk');
        end
    end
end

% Input validation for setInitialBounds
function testInputsInSetInitialBounds(testCase)
    validBoundsArray = {[0,1]};
    nValidBoundsArray = length(validBoundsArray);
    inValidBoundsArray = {[0,1;1,2],[1,0],[0,1]*1i,[-inf,inf],nan(1,2)};
    boundsArray = [validBoundsArray,inValidBoundsArray];
    nBoundsArray = length(boundsArray);
    f = @(x) x.^2-2*x-3;
    requiredUncertainty = 1e-3;
    for iBounds = 1:nBoundsArray
        if iBounds <= nValidBoundsArray
            [lowerFinalBound, upperFinalBound] = fibonacciSearch(f,...
                boundsArray{iBounds}(:,1),boundsArray{iBounds}(:,2),...
                requiredUncertainty);
            testCase.assertGreaterThan(upperFinalBound-lowerFinalBound,0);
        else
            testCase.assertError(@()fibonacciSearch(f,...
                boundsArray{iBounds}(:,1),boundsArray{iBounds}(:,2),...
                requiredUncertainty),...
                'fibonacciSearch:argChk');
        end
    end
end

% Input validation for setRequiredInterval
function testInputsInSetRequiredInterval(testCase)
    f = @(x) x.^2-2*x-3;
    lowerBound = 0;
    upperBound = 1;
    validRequiredIntervalArray = {1e-3,1e-6,0.1,0.50};
    nValidRequiredIntervalArray = length(validRequiredIntervalArray);
    invalidRequiredIntervalArray = {[0.2,0.3],nan,inf,0,2};
    requiredIntervalArray = [validRequiredIntervalArray,invalidRequiredIntervalArray];
    nRequiredIntervalArray = length(requiredIntervalArray);
    for iInterval = 1:nRequiredIntervalArray
        if iInterval <= nValidRequiredIntervalArray
            [lowerFinalBound, upperFinalBound] = fibonacciSearch(...
                f,lowerBound,upperBound,requiredIntervalArray{iInterval});
            testCase.assertGreaterThan(upperFinalBound-lowerFinalBound,0);
        else
            testCase.assertError(@()fibonacciSearch(...
                f,lowerBound,upperBound,requiredIntervalArray{iInterval}),...
                'fibonacciSearch:argChk');
        end
    end
end