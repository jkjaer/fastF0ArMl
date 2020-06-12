function [predErrVar, cholVct] = genSchurToeplitzAlg(rho)
    arOrder = length(rho)-1;
    predErrVar = rho;
    cholVct = nan(arOrder*(arOrder+1)/2,1);
    u = rho(1:end)/sqrt(rho(1));
    v = u(2:end);
    cholIdx = 1:arOrder;
    for p = 1:arOrder
       cholVct(cholIdx) = u(1:end-1);
       nu = v(1)/u(1);
       det = (1-nu^2);
       predErrVar(p+1) = predErrVar(p)*det;
       if p < arOrder
           sqrtDet = sqrt(det);
           u = (u(1:end-1)-nu*v)/sqrtDet;
           v = sqrtDet*v(2:end)-nu*u(2:end);
           cholIdx = cholIdx(1:end-1)+(arOrder-p+1);
       end
    end
end

% function [predErrVar, cholFact] = genSchurToeplitzAlg(rho)
%     arOrder = length(rho)-1;
%     predErrVar = rho;
%     cholFact = cell(arOrder,1);
%     u = rho(1:end)/sqrt(rho(1));
%     v = u(2:end);
%     for p = 1:arOrder
%        cholFact{p} = u(1:end-1);
%        nu = v(1)/u(1);
%        det = (1-nu^2);
%        predErrVar(p+1) = predErrVar(p)*det;
%        if p < arOrder
%            sqrtDet = sqrt(det);
%            u = (u(1:end-1)-nu*v)/sqrtDet;
%            v = sqrtDet*v(2:end)-nu*u(2:end);
%        end
%     end
% end
