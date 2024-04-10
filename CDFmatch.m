function[dsdata]=CDFmatch(cmdata,coarsedist,refdist,norain,p0,p0R)


dsdata = NaN(size(cmdata));

if ~isnan(cmdata)

    if p0 > p0R

        zero_sample = zeros(round((p0R/p0)*10000,0),1);
        a = 0.01;
        b = icdf(refdist,(p0-p0R)/(1-p0R)); % Finding corresponding
        nozero_sample = (b-a).*rand(round((1-p0R/p0)*10000,0),1)+a;
        rand_sample = cat(1,zero_sample,nozero_sample);
        ind_rand = randi([1 numel(rand_sample)],1,numel(cmdata(cmdata<=norain)));
        dsdata(cmdata<=norain) = rand_sample(ind_rand);

        y = (cdf(coarsedist,cmdata(cmdata>norain))*(1-p0)+p0-p0R)/(1-p0R);

        dsdata(cmdata>norain) = icdf(refdist,y);

    else
        dsdata(cmdata <= norain) = 0;

        y = cdf(coarsedist,cmdata(cmdata>norain));

        dsdata(cmdata>norain) = icdf(refdist,y);
    end
end
end