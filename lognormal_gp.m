function [best_fit,threshold,qnt_threshold,p0,A_sq_best, ind_min_A_sq] = ...
    lognormal_gp(tserie,pup,norain)


tserie(tserie<=norain) = 0;
p0 = numel(find(tserie==0))/numel(tserie);

sort_precip = sort(unique(tserie));

N = numel(unique(tserie));
m = (1:N).';

A_sq = NaN(length(pup),1);
pd_ln_gp = cell(length(pup),1);
% Fitting Distributions precipitation
for i = 1:length(pup)
    
    pd_ln_gp{i} = paretotails(tserie(tserie>0),0,pup(i),@mylognormal);
      
    
    temp = ((2.*m-1)/N).*(log(customcdf(p0,pd_ln_gp{i},sort_precip)) +...
            log(1-customcdf(p0,pd_ln_gp{i},flip(sort_precip))));
        
    
    A_sq(i) = -N - sum(temp);
    
end

[~, ind_min_A_sq] = min(A_sq);

A_sq_best = A_sq(ind_min_A_sq);

best_fit = pd_ln_gp{ind_min_A_sq};

qnt_threshold = pup(ind_min_A_sq)*(1-p0) + p0;

threshold = quantile(tserie(tserie>0),pup(ind_min_A_sq));

end