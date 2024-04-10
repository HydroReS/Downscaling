%% Downscale fields of variable x (e.g. precipitation) using mixed CDF matching
%
%Input: 1) Ref_field: 3D array (lat,lon,time) of reference field for
%estimation of empirical CDF
%       2) Ref_lat: 2D array of lat info for reference field
%       3) Ref_lon: 2D array of lon info for reference field
%       4) Coarse_field: 3D array(lat,lon,time) of coarse field for
%       estimation of empirical CDF. Note that in principle, Ref_field and
%       Coarse_field should overlap in time (i.e refer to the same period)
%       5) Coarse_lat: 2D array of lat info for coarse field
%       6) Coarse_lon: 2D array of lat info for coarse field
%       7) CS_to_DL_field: 3D array of coarse field to be
%       downscaled...MUST have same orientation with Coarse field!!!!
%
%IMPORTANT!!!!! 3rd dimension of the Ref,Coarse,CS_to_DL field must be time<------------
%
%Output: 1) DL_field:Downscaled field from coarse to reference resolution/spatial extent
%
%NOTE: a) The coarse field's extent must be greater or equal to reference
%field's extent.
%       b)Coarse_field and CS_to_DL_field must have equal extent
%       c) units of reference, coarse and CS_to_DL fields MUST have the
%       same units
%       d) missing values in field are denoted by NaN
%%
function [DL_field, Ref_best_fit, Ref_p0, Coarse_best_fit, Coarse_p0 ] = mixed_CDF_downscale_field(Ref_field,Ref_lat,Ref_lon,Coarse_field,Coarse_lat,Coarse_lon,CS_to_DL_field,pup,norain)

CS_to_DL_field = round(CS_to_DL_field,3);

tic

if(size(Coarse_field,1)~=size(CS_to_DL_field,1) || size(Coarse_field,2)~=size(CS_to_DL_field,2))
    error('dimensions of coarser fields (for calibration) and fields for downscaling should be consistent')
end

% Prelocating space
DL_field = single(NaN(size(Ref_field,1),size(Ref_field,2),size(CS_to_DL_field,3)));
Ref_best_fit = cell(size(Ref_field,1),size(Ref_field,2));
Coarse_best_fit = cell(size(Ref_field,1),size(Ref_field,2));
Ref_p0 = single(NaN(size(Ref_field,1),size(Ref_field,2)));
Coarse_p0 = single(NaN(size(Ref_field,1),size(Ref_field,2)));
toc
%for all rows and columns
for i=1:size(Ref_field,1)
    for j=1:size(Ref_field,2)
        fprintf('i=%d, j=%d\n',i,j);
        %check if at least 50% of the reference series is not NaN
        %and proceed with empirical CDF estimation and downscaling
        %otherwise assign NaN
         
        clear Ref_ecdf Coarse_ecdf
            
        % match locations between coarse and reference pixels         
        %Calculate distance between reference point Ref_field(i,j) and points in
        %Coarse field and identify the nearest point to match CDFs.
        x=Coarse_lon(:);%create column vector of all coordinates
        y=Coarse_lat(:);
        d=ipdm([x,y],[Ref_lon(i,j),Ref_lat(i,j)]);%calculate distance
        k=find(d==min(d));%find index of smallest distance
       [ii, jj]=ind2sub(size(Coarse_lat),k(1));%get corresponding row,col index for nearest point
            
        if(length(find(isnan(squeeze(Ref_field(i,j,:)))))<0.5*length(Ref_field(i,j,:))) && (length(find(isnan(squeeze(CS_to_DL_field(ii(1),jj(1),:)))))<0.5*length(CS_to_DL_field(ii(1),jj(1),:)))
        
	    toc
           
	    Ref = squeeze(Ref_field(i,j,:));

	    Coarse = squeeze(Coarse_field(ii(1),jj(1),:));

	    p0Ri = numel(find(Ref(Ref<=norain)))/numel(Ref);

	    [f,x] = ecdf(Coarse);

	    ind = find(f>=p0Ri);

	    p0Ci = f(ind(1));

	    norain_coarse = x(ind(1));
 
            [Ref_best_fit{i,j},~,~,Ref_p0(i,j),~,~] = lognormal_gp_mod(squeeze(Ref_field(i,j,:)),pup,norain,0.5);
            [Coarse_best_fit{i,j},~,~,Coarse_p0(i,j),~,~] = lognormal_gp_mod(squeeze(Coarse_field(ii(1),jj(1),:)),pup,norain_coarse,0.5);
            
            toc
             
            DL_field(i,j,:) =CDFmatch(squeeze(CS_to_DL_field(ii(1),jj(1),:)),Coarse_best_fit{i,j},Ref_best_fit{i,j},norain_coarse,Coarse_p0(i,j),Ref_p0(i,j));
            
            toc
        else
            DL_field(i,j,:) = NaN*ones(size(CS_to_DL_field,3),1);
	    toc	
        end
    end
end


end
