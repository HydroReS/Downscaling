function[dsdata]=CDFmatchTemp(cmdata,coarsedist,refdist)

dsdata = NaN(size(cmdata));

if ~isnan(cmdata)
    
    y = cdf(coarsedist,cmdata);
    
    dsdata = icdf(refdist,y);
    
end
end