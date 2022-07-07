function [isr,psr]=run_fingerprinting(data,n_scans,ref_scan,target_scan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run fingerprinting to report identification success and perfect separation rates
% by Stephanie Noble and Jean Ye 07/07/2022
%
% input:
%   data = n x n x m, n=# nodes, m=# subs x # scans
%       ordered like: sub1_scan1 sub1_scan2 sub2_scan1 sub2_scan2
%
% output:
%   isr = identification success rate (for each subject, there is a perfect match)
%       currently built for only a single reference and target
%       TODO: generalize to multiple targets
%   psr = perfect separation rate (for each subject, all match
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n_nodes1,n_nodes2,n_subs_by_scans]=size(data);
if n_nodes1~=n_nodes2
    error('Data must be n_nodes by n_nodes by n_scans.')
end
n_nodes=n_nodes1;

n_subs=n_subs_by_scans/n_scans;

data_vec=zeros((n_nodes-1)*n_nodes/2,n_subs_by_scans);
trilmask_data=logical(tril(ones(n_nodes),-1));
for i=1:n_subs_by_scans
    d=data(:,:,i);
    data_vec(:,i)=d(trilmask_data);
end

corr_all_mats=corr(data_vec);
trilmask=logical(tril(ones(n_scans),-1));


for this_sub=1:n_subs

    % track perfect separations
    idx_w_start=(this_sub-1)*n_scans+1;
    idx_w_end=(this_sub)*n_scans;

    corr_within=corr_all_mats(idx_w_start:idx_w_end,idx_w_start:idx_w_end);
    corr_within=corr_within(trilmask);

    corr_between=corr_all_mats(idx_w_start:idx_w_end,:);
    corr_between(:,idx_w_start:idx_w_end)=[];

    ps(this_sub)=min(corr_within(:))>max(corr_between(:));

    % track identification successes
    idx_ref=(this_sub-1)*n_scans+ref_scan;
    idx_target=[0:n_subs-1]*n_scans+target_scan;

    corr_all__ref_scan=corr_all_mats(idx_ref,:);
    corr_all__ref_scan=corr_all__ref_scan(idx_target);
    max_id=find(corr_all__ref_scan==max(corr_all__ref_scan));

    id(this_sub)=max_id==this_sub;

    % springboard for permutation test
%     corr_all__ref_scan=corr_all__ref_scan(randperm(length(corr_all__ref_scan)));
%     max_id__perm=find(corr_all__ref_scan==max(corr_all__ref_scan));

end

psr=sum(ps)/n_subs*100; % perfect separation rate
isr=sum(id)/n_subs*100; % identification success rate
