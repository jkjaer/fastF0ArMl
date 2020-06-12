function arParameters = genRndArParams(order, rootAbsMin, rootAbsMax)
    nRealRoots = mod(order,2);
    nComplexRootPairs = (order-nRealRoots)/2;
    % generate the real roots
    realRoots = sign(randn(nRealRoots,1)).*(...
        (rootAbsMax-rootAbsMin)*rand(nRealRoots,1)+rootAbsMin);
    % generate the complex roots
    complexRoots = ...
        ((rootAbsMax-rootAbsMin)*rand(nComplexRootPairs,1)+rootAbsMin).*...
        exp(1i*pi*rand(nComplexRootPairs,1));
    % compute AR coefficients
    b = poly([realRoots;complexRoots;conj(complexRoots)])';
    arParameters = -b(2:end);
end
