%% fibonacciSearch
% Single variable bounded minisation using a Fibonacci search
%
%% Syntax:
% [lowerFinalBound, upperFinalBound] = fibonacciSearch(objectiveFunction,...
%     lowerInitialBound, upperInitialBound, requiredInterval);
%
%% Description:
% Narrows down the minimum of a function from an initial interval to smaller
% interval using af Fibonacci search.
% The implementation is based on Sec. 4.3 in 
% A. Antoniou and W.-S. Lu, Practical Optimization: Algorithms and Engineering
% Applications. Springer, Mar. 2007.
% Note that the index for the Fibonacci number used in the book is 1
% smaller than the index used below.
% * objectiveFunction: A function handle of the function to be minimised
% * lowerInitialBound: The lower initial bound of the minimiser
% * upperInitialBound: The upper initial bound of the minimiser
% * requiredInterval: The required final interval of uncertainty
%
%% Examples:
% f = @(x) x.^2-2*x-3;
% lowerBound = -2;
% upperBound = 5;
% requiredUncertainty = 1e-3;
% [lowerBound, upperBound] = fibonacciSearch(f,lowerBound,upperBound,...
%     requiredUncertainty);
% estimatedMinimiser = (upperBound+lowerBound)/2;
%
function [lowerFinalBound, upperFinalBound] = fibonacciSearch(objectiveFunction,...
        lowerInitialBound, upperInitialBound, requiredInterval)
    % Check the input arguments
    if nargin<4
        error('fibonacciSearch:argChk','Four input arguments are required.');
    end
    validateObjectiveFunction(objectiveFunction);
    validateBound(lowerInitialBound);
    validateBound(upperInitialBound);
    initialInterval = upperInitialBound-lowerInitialBound;
    if initialInterval<=0
        error('fibonacciSearch:argChk',...
            'The lower bound must be smaller than the upper bound.');
    end
    validateInterval(requiredInterval);
    % Compute the list of Fibonacci numbers
    fibonacciNoList = computeFibonacciNumbers(initialInterval,requiredInterval);
    % Perform the Fibonacci search
    [lowerFinalBound, upperFinalBound] = narrowBounds(objectiveFunction,...
        lowerInitialBound, upperInitialBound, fibonacciNoList);
end

% Return an error if the objective function is invalid
function validateObjectiveFunction(objectiveFunction)
    if ~isa(objectiveFunction, 'function_handle')
        error('fibonacciSearch:argChk',...
                'The input argument must be a function handle');
    end
end

% Return an error if the bound is invalid
function validateBound(bound)
    if ~isscalar(bound) || ~isfinite(bound)
        error('fibonacciSearch:argChk',...
            'The bounds must be real-valued scalars.');
    end
end

% Return an error if the interval is invalid
function validateInterval(interval)
    if ~isscalar(interval) || ~isfinite(interval) ||...
            interval<=0
        error('fibonacciSearch:argChk',...
            'The required interval must be a positive real-valued scalar.');
    end
end

% Compute the list of Fibonacci numbers so that the initial interval is narrowed
% down to an interval not greater than the required interval
function fibonacciNoList = computeFibonacciNumbers(initialInterval,...
        requiredInterval)
    % The factor of two ensures that the Fibonacci search produces a final
    % interval smaller than requiredInterval instead of 2*requiredInterval.
    intervalRatio = initialInterval/requiredInterval;
    if intervalRatio < 2
        error('fibonacciSearch:argChk',...
        ['The final interval must be smaller than half the ',...
        'initial interval.']);
    end
    nFibonacci = fibonacciIndex(intervalRatio);
    if nFibonacci > 70
        error('fibonacciSearch:argChk',...
        ['The required final interval requires a Fibonacci number ',...
        'of index greater than 70. This cannot be computed reliably.']);
    end
    fibonacciNoList = fibonacci(2:nFibonacci);
end

% Run the Fibonacci search algorithm to narrow down the lower and upper bound 
% for the minimiser to within a tolerance of at most +/-requiredInterval/2. 
% The algorithm is based on algorithm 4.1 in A. Antoniou and W.-S. Lu, 
% Practical Optimization: Algorithms and Engineering Applications. 
% Springer, Mar. 2007.
function [lowerBound, upperBound] = narrowBounds(objectiveFunction,...
        lowerBound, upperBound, fibonacciNoList)
    nFibonacci = length(fibonacciNoList);
    startInterval = upperBound-lowerBound;
    iInterval = startInterval*fibonacciNoList(nFibonacci-1)/...
        fibonacciNoList(nFibonacci);
    variableLowerVal = upperBound-iInterval;
    funcLowerVal = objectiveFunction(variableLowerVal);
    variableUpperVal = lowerBound+iInterval;
    funcUpperVal = objectiveFunction(variableUpperVal);
    nMaxIterations = nFibonacci-2;
    % Run the first nFibonacci-2 iterations
    for iIteration = 1:nMaxIterations
        if iIteration<nMaxIterations
            iInterval = iInterval*fibonacciNoList(nFibonacci-iIteration-1)/...
                fibonacciNoList(nFibonacci-iIteration);
        else
            % To avoid that the interior points are the same after the last
            % iteration, we multiply by slightly more than 0.5.
            iInterval = 1.01*iInterval/2;
        end
        if funcLowerVal > funcUpperVal
            % The minimum is in the interval [variableLowerVal;upperBound]
            lowerBound = variableLowerVal;
            variableLowerVal = variableUpperVal;
            variableUpperVal = lowerBound+iInterval;
            funcLowerVal = funcUpperVal;
            funcUpperVal = objectiveFunction(variableUpperVal);
        else
            % The minimum is in the interval [lowerBound;variableUpperVal]
            upperBound = variableUpperVal;
            variableUpperVal = variableLowerVal;
            variableLowerVal = upperBound-iInterval;
            funcUpperVal = funcLowerVal;
            funcLowerVal = objectiveFunction(variableLowerVal);
        end
        % If the final search tolerance is within the precision of
        % the computer, then stop
        if variableLowerVal > variableUpperVal
            break;
        end
    end
    % In the last iteration, only the final bounds are used so we
    % do not waste computational resources on updating the other
    % parameters
    if funcLowerVal > funcUpperVal
        lowerBound = variableLowerVal;
    else
        upperBound = variableUpperVal;
    end
end
