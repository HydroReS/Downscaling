function [best_fit,qnt_threshold,threshold,A_sq_best] = ...
    gaussian_gp(tserie,pup,pl)

sort_precip = sort(unique(tserie));

N = numel(unique(tserie));
m = (1:N).';

A_sq = NaN(length(pup),length(pl));
pd_ln_gp = cell(length(pup),length(pl));
% Fitting Distributions precipitation
for i = 1:length(pup)
    for j = 1:length(pl)
        
        pd_ln_gp{i,j} = paretotails(tserie,pl(j),pup(i),@mygaussian);        
        
        temp = ((2.*m-1)/N).*(log(cdf(pd_ln_gp{i,j},sort_precip)) +...
            log(1-cdf(pd_ln_gp{i,j},flip(sort_precip))));
        
        A_sq(i,j) = -N - sum(temp);
    end
end

[ind_min_A_sq_i,ind_min_A_sq_j] = find(A_sq == min(A_sq,[],'all'));

best_fit = pd_ln_gp{ind_min_A_sq_i,ind_min_A_sq_j};

A_sq_best = A_sq(ind_min_A_sq_i,ind_min_A_sq_j);

qnt_threshold(1) = pup(ind_min_A_sq_i);

qnt_threshold(2) = pl(ind_min_A_sq_j);

threshold(1) = quantile(tserie,pup(ind_min_A_sq_i));
threshold(2) = quantile(tserie,pl(ind_min_A_sq_j));

end