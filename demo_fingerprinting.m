% fingerprinting demo
% generates data and runs fingerprinting script


% set params

n_subs=4;
n_scans=5;
n_nodes=268;
amt_noise=1; % higher is less - 0.001 for 1 scan to get id % 1 for multiple; 0.001 for multiple scans can give id=0.25 and psr=0

% generate data

data=zeros(n_nodes,n_nodes,n_subs*n_scans);

for i=1:n_subs
    t=rand(n_nodes);
    for j=1:n_scans
        t_noise=(rand(n_nodes,n_nodes)-0.5)/amt_noise;
        data(:,:,(i-1)*n_scans+j)=t + t_noise;
    end
end

% visualize everything

[n_nodes,~,n_subs_by_scans]=size(data);
n_subs=n_subs_by_scans/n_scans;
data_vec=reshape(data,n_nodes*n_nodes,n_subs_by_scans);
imagesc(data_vec)

corr_all_mats=corr(data_vec);
figure; imagesc(corr_all_mats)

ref_scan=1; target_scan=2;
[isr,psr]=run_fingerprinting(data,n_scans,ref_scan,target_scan);

fprintf('Perfect Separation Rate: %0.1f\nIdentification Success Rate=%0.1f\n',psr,isr)

