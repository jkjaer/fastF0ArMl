function [predErrVar, arParams] = myLpc(dataVector, arOrder)
    dataVector = [dataVector;zeros(arOrder,1)];
    nData = length(dataVector);
    % buffers
    predErrVar = nan(arOrder+1,1);
    rho = nan(arOrder,1);
    arParams = nan(arOrder,arOrder);
    % p == 0
    predErrVar(1) = (dataVector'*dataVector)/nData;
    % p == 1
    if arOrder > 0
        rho(1) = (dataVector(2:end)'*dataVector(1:end-1))/nData;
        nu = rho(1)/predErrVar(1);
        predErrVar(2) = predErrVar(1)*(1-nu^2);
        arParams(1,1) = nu;
        for p = 2:arOrder % p >= 1
            rho(p) = (dataVector(p+1:end)'*dataVector(1:end-p))/nData;
            nu = (rho(p)-flipud(arParams(1:p-1,p-1))'*rho(1:p-1))/predErrVar(p);
            predErrVar(p+1) = predErrVar(p)*(1-nu^2);
            arParams(1:p,p) = [arParams(1:p-1,p-1);0]+...
                nu*[-flipud(arParams(1:p-1,p-1));1];
        end
    end
end
