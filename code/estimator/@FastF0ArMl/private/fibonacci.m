%% fibonacci 
% Compute Fibonacci numbers
%% Syntax:
%# fibonacciNo = fibonacci(index)
%
%% Description:
% Compute the Fibonacci numbers corresponding the indices in the
% index vector. For index n, the Fibonacci number is
%
%   $$F(n) = F(n-1)+F(n-2)$$
%
% with $$F(0) = 0$$ and $$F(1) = 1$$. For negative indices, the recursion can be
% rearraged as
%
%   $$F(n-2) = F(n)-F(n-1)$$
%
% The indices must be real-valued integers with a magnitude
% smaller than 71.
% * index: scalar, vector, og matrix of indices
%
%% Examples:
%# fibonacciNo = fibonacci(10)
%# fibonacciNos = fibonacci([1:5,15,-3])
function fibonacciNo = fibonacci(index)
    % If any index elements are +/- infinite, complex-valued, not an
    % integer, or larger than abs(70), issue an error
    indexElementIsNonInteger = ~isreal(index) | ...
        abs(index-round(index))>1e-16;
    indexElementIsTooLarge = abs(index)>70;
    indexIsInvalid = indexElementIsNonInteger | indexElementIsTooLarge;
    if any(indexIsInvalid)
        error('fibonacci:argChk',...
            ['An infinite, complex-valued, non integer, or too large ',...
            '(index > 70) index has been given.']);
    else
        goldenRatio = 0.5+sqrt(5)/2;
        fibonacciNo = round(...
            (goldenRatio.^index-(1-goldenRatio).^index)/sqrt(5));
    end
end