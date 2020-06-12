%% fibonacciIndex
% Compute the index of the smallest Fibonacci number upper bounding a
% real-valued number.
%% Syntax:
%# index = fibonacciIndex(positiveNumber)
%
%% Description:
% Compute the index $n$ of the smallest Fibonacci number $F(n)$ upper
% bounding a real-valued number $2<=x<=F(70)$. That is, if the number $x$
% satisfies
%
%   $$F(n-1) < x <= F(n)$$
%
% the index $n$ is returned. The function also accept a vector or matrix of
% numbers as inputs.
% * positiveNumber: scalar, vector, or matrix of postive numbers
%
%% Examples:
%# index = fibonacciIndex(77)
%# indexVector = fibonacciIndex([2,47])
function indexVector = fibonacciIndex(positiveNumber)
    fibonacci70 = 190392490709135; % Fibonacci number for index 70
    % If any index elements are +/- infinite, complex-valued, too small, or 
    % too large, return an error
    fibonacciNoIsTooSmall = (positiveNumber<2);
    fibonacciNoIsTooBig = (positiveNumber>fibonacci70);
    fibonacciNoIsNonInteger = ~isreal(positiveNumber) | isinf(positiveNumber);
    fibonacciNoIsInvalid = fibonacciNoIsTooSmall | fibonacciNoIsTooBig |...
        fibonacciNoIsNonInteger;
    if any(fibonacciNoIsInvalid)
        error('fibonacciIndex:argChk',...
            ['A infinite, complex-valued, too small (< 2) ',...
            'Fibonacci number has been given.']);
    else
        goldenRatio = 0.5+sqrt(5)/2;
        indexVector = floor(log(positiveNumber*sqrt(5)+0.5)/log(goldenRatio));
        % If the input number is not a Fibonacci number and larger than the
        % Fibonacci number with the computed index, add one to the index.
        fibonacciNo = fibonacci(indexVector);
        numberIsFibonacciNo = (abs(fibonacciNo-positiveNumber)<1e-16) | ...
            isnan(positiveNumber);
        if any(~numberIsFibonacciNo)
            indexVector(~numberIsFibonacciNo) = indexVector(~numberIsFibonacciNo)+...
                double(fibonacciNo(~numberIsFibonacciNo)<...
                positiveNumber(~numberIsFibonacciNo));
        end
    end
end