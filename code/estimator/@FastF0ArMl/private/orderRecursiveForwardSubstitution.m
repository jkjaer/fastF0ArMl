function resVar = orderRecursiveForwardSubstitution(cholVct, qRho)
    [maxArOrder, nFreqs] = size(qRho);
    resVar = nan(maxArOrder,nFreqs);
    cholIdx = nan(maxArOrder);
    linParam = nan(maxArOrder,nFreqs);
    for p = 1:maxArOrder
        if p == 1
            cholIdx(:,1) = 1:maxArOrder;
            linParam(1,:) = qRho(1,:)./cholVct(cholIdx(1,1),:);
            resVar(1,:) = linParam(1,:).^2;
        else
            cholIdx(p:end,p) = cholIdx(end,p-1)+(1:maxArOrder-p+1);
            linParam(p,:) = (qRho(p,:)-sum(...
                cholVct(cholIdx(p,1:p-1),:).*linParam(1:p-1,:),1))./...
                cholVct(cholIdx(p,p),:);
            resVar(p,:) = resVar(p-1,:)+linParam(p,:).^2;
        end
    end
end