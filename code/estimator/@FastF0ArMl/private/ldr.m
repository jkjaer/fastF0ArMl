function [predErrVar, rho, nu] = ldr(dataVector, arOrder, nData)
    % buffers
    rho = nan(arOrder,1);
    predErrVar = nan(arOrder+1,1);
    nu = nan(arOrder,1);
    % p == 0
    predErrVar(1) = (dataVector'*dataVector)/nData;
    % p == 1
    if arOrder > 0
        rho(1) = (dataVector(2:end)'*dataVector(1:end-1))/nData;
        nu(1) = rho(1)/predErrVar(1);
        predErrVar(2) = predErrVar(1)*(1-nu(1)^2);
        % we run the recursion in the reversed version beta to avoid one flipud
        betaRev = nu(1);
        for p = 2:arOrder % p >= 1
            rho(p) = (dataVector(p+1:end)'*dataVector(1:end-p))/nData;
            nu(p) = (rho(p)-betaRev'*rho(1:p-1))/predErrVar(p);
            predErrVar(p+1) = predErrVar(p)*(1-nu(p)^2);
            if p < arOrder
                betaRev = [0;betaRev]+nu(p)*[1; -flipud(betaRev)];
            end
        end
    end
end