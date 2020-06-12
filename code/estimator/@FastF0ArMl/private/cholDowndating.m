function cholVct = cholDowndating(cholVct, v)
    nRows = size(v,1);
    cholIdx = 1:nRows;
    for p = 1:nRows
        nu = ones(nRows-p+1,1)*(v(1,:)./cholVct(cholIdx(1),:));
        sqrtDet = sqrt(1-nu.^2);
        cholVct(cholIdx,:) = (cholVct(cholIdx,:)-nu.*v)./sqrtDet;
        if p < nRows
            v = sqrtDet(2:end,:).*v(2:end,:)-...
                nu(2:end,:).*cholVct(cholIdx(2:end),:);
            cholIdx = cholIdx(1:end-1)+(nRows-p+1);
        end
    end
end