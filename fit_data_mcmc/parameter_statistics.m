function [parametertable,matertable] = parameter_statistics(pfull_final,logP,func_struct,burnin,confidence_level,names)

Nwalkers = size(pfull_final,2);

% Select walkers past the burn-in
p = pfull_final(:,:,round(burnin*size(pfull_final,3)));

% Reshape to 2D array
p = p(:,:);

% Reorder parameters in each row lowest to highest
sortp = sort(p,2);

% Get confidence region
alpha = 1 - confidence_level;
leftCR = sortp(:,floor((alpha/2)*size(sortp,2)));
rightCR = sortp(:,ceil((1-alpha/2)*size(sortp,2)));

% Get mean parameter values
p_mean = mean(p,2);

% Get median parameter values
p_median = median(p,2);

% Get MLE parameters
[~,MLEind] = max(logP(2,:,end));
p_MLE = pfull_final(:,MLEind,end);

% Combine into a table
parametertable = array2table([leftCR,p_median,p_mean,p_MLE,rightCR],'VariableNames',{'LeftCR','Median','Mean','MLE','RightCR'},...
    'RowNames',names);
display(parametertable)
checkoutsiderange = (p_MLE > rightCR) | (p_MLE < leftCR);
if sum(checkoutsiderange) > 0
    warning('MLE is outside confidence region. Either not run long enough to be sampling a stationary distribution or a parameter is pushing up against a bound. Examine burn-in plot and corner plot and run again')
end

%% Get atomic ratios

atom_perc_mat = zeros(2,Nwalkers);
[lig_cov_vec,equiv_d_vec,core_ratio_vec] = deal(zeros(1,Nwalkers));

for i = 1:Nwalkers
    xwalk = pfull_final(:,i,end);
    [atomic_perc,lig_cov,equiv_d,core_ratio] = calc_atomic_ratios(xwalk,func_struct);
    atom_perc_mat(:,i) = atomic_perc;
    lig_cov_vec(i) = lig_cov;
    equiv_d_vec(i) = equiv_d;
    core_ratio_vec(i) = core_ratio;
end

combined_mater = [atom_perc_mat;core_ratio_vec;equiv_d_vec;lig_cov_vec];
sorted_combined = sort(combined_mater,2);

% Get confidence region
leftCR_mater = sorted_combined(:,floor((alpha/2)*Nwalkers));
rightCR_mater = sorted_combined(:,ceil((1-alpha/2)*Nwalkers));

% Get mean parameter values
mean_mater = mean(combined_mater,2);

% Get median parameter values
median_mater = median(combined_mater,2);

% Get MLE parameters
MLE_mater = combined_mater(:,MLEind);

% Combine into a table
matertable = array2table([leftCR_mater,median_mater,mean_mater,MLE_mater,rightCR_mater],'VariableNames',{'LeftCR','Median','Mean','MLE','RightCR'},...
    'RowNames',{'Pb %','S %','Core Pb/S Ratio','Equiv diam','Lig Cov'});
display(matertable)

end

