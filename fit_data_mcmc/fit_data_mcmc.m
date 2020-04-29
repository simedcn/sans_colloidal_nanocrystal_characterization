%% Outline
% Script to fit the SANS data using the MCMC approach. This code:
%   (0) Chooses MCMC specifications
%   (1) Loads data
%   (2) Specifies materials parameters
%   (3) Gets the roots and weights needed for Gauss-Hermite quadrature
%   (4) Retrieves model and comprising functions
%   (5) Reformats data for use in the MCMC fit
%   (6) Specifies log-likelihood and log-prior
%   (7) Initialize walkers for preliminary run
%   (8) Run preliminary fit and save results
%   (9) Plot preliminary fit
%   (10) Reinitialize walkers for production run
%   (11) Run production fit and save results
%   (12) Plot production run
%   (13) Plot corner plot from production run

%% (0) User-chosen MCMC specifications

% MCMC Specifications
n_gh = 12; % Number of Gauss-Hermite quadrature nodes
Nwalkers_gen = 14^2; % Generate more walkers than necessary. Not all will satisfy constraints.
Nwalkers_target = 100;
prelim_steps_per_walker = 300;
production_steps_per_walker = 300;
stepsize = 1.25;
MCMC_burnin = 0;
manual_burnin = 0.75;
confidence_level = 0.95;
parallel_spec = false;

% Output specs (true or false)
plot_prelim = true;
plot_production = true;
plot_final_corner = true;

%% (1) Load data

datastruct = load('sans_data.mat');
data = datastruct.data;
q = data.q;                 % [218 x 1]
Idata = data.I_bkgdcorr;    % [218 x 8]
Istd = data.sigI;           % [218 x 8]
expt_sld = data.SLD;        % [1 x 8]
fd8 = data.fd8;             % [1 x 8]
sig_q = data.sig_q;         % [218 x 1]
mean_q = data.mean_q;       % [218 x 1]
shadow = data.shadow;       % [218 x 1]

%% (2) Material specifications

rho_PbS = 2.343e-6; % [Å^-2]
rho_head = 3.9e-6; % [Å^-2]
rho_OA = 0.078e-6; % [Å^-2]
rho_tol = 0.941e-6; % [Å^-2]
rho_d8 = 5.664e-6; % [Å^-2]

%% (3) Get smearing roots and weights
[x_gh_vec,w_gh_vec] = GaussHermite( n_gh ); % Function from Geert van Damme

%% (4) Get model and comprising functions

[call_model,func_struct] = define_model();

%% (5) Reformat data

Ibatch = Idata(:,:);
Istd_batch = Istd(:,:);
Nd8 = length(fd8);
expt_sldmat = repmat(expt_sld,length(q),1);
sld_data = expt_sldmat(:)';
sld_data = repmat(sld_data,n_gh,1);
q_data = repmat(q',n_gh,Nd8);
mean_q_data = repmat(mean_q',n_gh,Nd8);
sig_q_data = repmat(sig_q',n_gh,Nd8);
shadow_data = repmat(shadow',n_gh,Nd8);
x_gh = repmat(x_gh_vec,1,size(q_data,2));
w_gh = repmat(w_gh_vec,1,size(q_data,2));
Iq_mean = Ibatch(:)';
Iq_std = Istd_batch(:)';

%% (6) Make log-likelihood and log-prior functions

% Parameters
plotnames = {'r_{min} (Å)','t_{vertex} (Å)','t_{plateau} (Å)','t_{max} (Å)','f_{PbS}','f_{OA}','\phi','bkgd (cm^{-1})'};
names = {'r_min','t_vertex','t_plateau','t_max','f_PbS','f_OA','phi','bkgd'};

% Bounds
mintvertexratio = 0.2393;
maxtvertexratio = sqrt(3) - 1;
rmin_low = 10; rmin_high = 45;
tvertex_low = mintvertexratio*rmin_low; tvertex_high = maxtvertexratio*rmin_high;
tplateau_low = 5; tplateau_high = 20;
tmax_low = 0; tmax_high = 5;
fPbS_low = 0; fPbS_high = 1;
fOA_low = 0; fOA_high = 1;
phi_low = 0.005; phi_high = 0.04;
bkgd_low = 0.003; bkgd_high = 0.04;
ub = [rmin_high, tvertex_high, tplateau_high, tmax_high, fPbS_high, fOA_high, phi_high, bkgd_high];
lb = [rmin_low, tvertex_low, tplateau_low, tmax_low, fPbS_low, fOA_low, phi_low, bkgd_low];
numparams = length(ub);

% Log-prior including non-linear constraints
fhead = func_struct.fhead_eq;
nonlcon = @(x) fhead(x) >= 0 && fhead(x) <= 1 && x(5) + fhead(x) <= 1 && x(2) <= maxtvertexratio*x(1) && x(2) >= 2*mintvertexratio*x(1);
logprior = makelogprior(lb,ub,nonlcon);

% Log-likelihood
loglike = @( x ) -sum( ( Iq_mean - call_model(q_data,x,rho_PbS,rho_OA,rho_head,sld_data,mean_q_data,sig_q_data,shadow_data,x_gh,w_gh) ) .^ 2 ./ Iq_std .^ 2 );

%% (7) Initialize walkers for preliminary run

% Spread r_min and t_vertex uniformly
c_walkers_vec = linspace(rmin_low,rmin_high,14);
[c_walkers_mat,tvertex_walkers_mat] = deal(zeros(14,14));
for i = 1:14
    c_walkers_mat(:,i) = c_walkers_vec;
    tvertex_walkers_mat(i,:) = linspace(2*mintvertexratio*c_walkers_vec(i),maxtvertexratio*c_walkers_vec(i),14);
end
c_walkers = c_walkers_mat( : )';
tvertex_walkers = tvertex_walkers_mat( : )';

% Randomly sample t_plateau, t_max, f_OA, phi, and bkgd
tplateau_walkers = unifrnd(tplateau_low,tplateau_high,size(c_walkers));
tmax_walkers = unifrnd(tmax_low,tmax_high,size(c_walkers));
fOA_walkers = unifrnd(fOA_low,fOA_high,size(c_walkers));
phi_walkers = unifrnd(phi_low,phi_high,size(c_walkers));
bkgd_walkers = unifrnd(bkgd_low,bkgd_high,size(c_walkers));

% Keep track of walkers that have unphysical values of f_head
[fhead_placeholders] = zeros(size(c_walkers));
track = [];
for i = 1:length(fhead_placeholders)
    fhead_placeholders(i) = fhead([c_walkers(i),tvertex_walkers(i),tplateau_walkers(i),tmax_walkers(i),0,fOA_walkers(i)]);
    if fhead_placeholders(i) < 0 || fhead_placeholders(i) > 1
        track = [track,i];
    end
end

% Calculate f_PbS that satisfy mass balance summing to 1
fPbS_walkers = unifrnd(zeros(size(c_walkers)),1 - fhead_placeholders);

% Remove unphysical walkers
xinit = [ c_walkers; tvertex_walkers; tplateau_walkers; tmax_walkers; fPbS_walkers; fOA_walkers; phi_walkers; bkgd_walkers];
xinit(:,track) = [];
xinit = xinit(:,1:Nwalkers_target);

%% (8) Run preliminary fit and save results

% Run MCMC and time it
tic
[xfull_prelim,logP_prelim,Rej_prelim] = gwmcmc( xinit, { logprior loglike }, Nwalkers_target*prelim_steps_per_walker, 'burnin', MCMC_burnin, 'stepsize', stepsize ,'Parallel',logical(parallel_spec));
prelim_time = toc;

% Save preliminary output
[parametertable_prelim,matertable_prelim] = parameter_statistics(xfull_prelim,logP_prelim,func_struct,manual_burnin,confidence_level,names);
save('run_data_prelim.mat','parametertable_prelim','matertable_prelim','xfull_prelim','logP_prelim','Rej_prelim','prelim_time')
    
%% (9) Plot preliminary fit

if plot_prelim
    run_data = load('run_data_prelim.mat');
    xfull_prelim = run_data.xfull_prelim;
    logP_prelim = run_data.logP_prelim;
    [call_model,func_struct] = define_model();
    plot_MCMC_output(q,Ibatch,Istd_batch,expt_sld,xfull_prelim,logP_prelim,call_model,func_struct,data,parametertable_prelim,manual_burnin,false)
end

%% (10) Reinitialize walkers for production run

prelim_data = load('run_data_prelim.mat');
xwalkers = prelim_data.xfull_prelim;
logP_reload = prelim_data.logP_prelim;
xinit_restart = get_restart_walkers(xwalkers,logP_reload);

%% (11) Run production fit and save results

% Run MCMC and time it
tic
[xfull_final,logP_final,Rej_final] = gwmcmc( xinit_restart, { logprior loglike }, Nwalkers_target*production_steps_per_walker, 'burnin', MCMC_burnin, 'stepsize', stepsize ,'Parallel',logical(parallel_spec));
production_time = toc;

% Save final output
[parametertable_final,matertable_final] = parameter_statistics(xfull_final,logP_final,func_struct,manual_burnin,confidence_level,names);
save('run_data_production.mat','parametertable_final','matertable_final','xfull_final','logP_final','Rej_final','production_time')

%% (12) Plot production run

if plot_production
    run_data = load('run_data_production.mat');
    xfull_final = run_data.xfull_final;
    logP_final = run_data.logP_final;
    parametertable_final = run_data.parametertable_final;
    [call_model,func_struct] = define_model();
    plot_MCMC_output(q,Ibatch,Istd_batch,expt_sld,xfull_final,logP_final,call_model,func_struct,data,parametertable_final,manual_burnin,true)
end

%% (13) Plot corner plot from production run

if plot_final_corner
    run_data = load('run_data_production.mat');
    xfull_production = run_data.xfull_final;
    parametertable_final = run_data.parametertable_final;
    x_median = parametertable_final.Median;
    x_mle = parametertable_final.MLE;
    left_CR = parametertable_final.LeftCR;
    right_CR = parametertable_final.RightCR;
    cornerplot(xfull_production(:,:,round(manual_burnin*size(xfull_production,3))),x_median,x_mle,left_CR,right_CR,plotnames,ones(size(lb)));
end