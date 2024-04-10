%% Combined CDF

function mycdf = customcdf(p0,dist1,x)

mycdf = NaN(size(x));

mycdf(x == 0) = p0;
mycdf(x ~= 0) = p0+(1-p0)*cdf(dist1,x(x~=0));

end