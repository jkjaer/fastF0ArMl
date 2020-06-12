function covVct = computeCovVct(dataVector, maxLag, nData)
    covVct = nan(maxLag+1,1);
    for ii = 0:maxLag
        covVct(ii+1) = (dataVector(ii+1:end)'*dataVector(1:end-ii))/nData;
    end
end