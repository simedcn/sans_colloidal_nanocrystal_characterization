function [] = plot_MCMC_output(xdata,ydata,yerror,expt_sld,xfull,logP,call_model,func_struct,data,paramtable,burnin,production_run)

%% Useful for plotting
rho_PbS = 2.343e-6; % [Å^-2]
rho_head = 3.9e-6;
rho_OA = 0.078e-6;
rho_tol = 0.941e-6; % [Å^-2]
rho_d8 = 5.664e-6; % [Å^-2]

sig_q = data.sig_q;
mean_q = data.mean_q;
shadow = data.shadow;
n_gh = 12;
[x_gh_vec,w_gh_vec] = GaussHermite( n_gh );
Nd8 = length(expt_sld);
mean_q_data = repmat(mean_q',n_gh,1);
sig_q_data = repmat(sig_q',n_gh,1);
shadow_data = repmat(shadow',n_gh,1);
x_gh = repmat(x_gh_vec,1,size(xdata,2));
w_gh = repmat(w_gh_vec,1,size(xdata,2));

%% Extract parameters
% Time progression
final_time = size(xfull,3);
time = 1:1:final_time;
Nwalkers = size(xfull,2);

%% Burn-in plot

figure();
for i = 1:size(logP,2)
    logPwalker = squeeze(logP(2,i,:));
    plot(time,logPwalker)
    hold on
end
h = gca;
if production_run
    plot(round(burnin*final_time)*[1 1],h.YLim,'--k')
    titlespec = 'Production run';
else
    titlespec = 'Preliminary run';
end   
xlabel('t (every 10 steps per walkers)')
ylabel('Log-likelihood')
title(titlespec)

%% Fit to data

cmap = [linspace(0,1,Nd8);linspace(0.8,0,Nd8);ones(1,Nd8)]';
[~,MLEind] = max(logP(2,:,end));

% Plot data
figure();
for i = 1:Nd8
    errorbar(xdata,ydata(:,i),yerror(:,i),'o','MarkerSize',2, 'LineWidth',1.5,'Color',cmap(i,:));
    hold on
end

% Plot last position of all walkers
for i = 1:Nwalkers
    p = xfull(:,i,end);
    for j = 1:Nd8
        sld_j = expt_sld(j);
        yfit = call_model(xdata',p,rho_PbS,rho_OA,rho_head,sld_j*ones(size(xdata')),mean_q_data,sig_q_data,shadow_data,x_gh,w_gh);
        plot(xdata',yfit,'-','Color',[0 0 0 0.15],'LineWidth',2)
    end
end

% Plot best performer
xMLE = xfull(:,MLEind,end);
xmedian = paramtable.Median;
for j = 1:Nd8
    sld_j = expt_sld(j);
    ymedian = call_model(xdata',xmedian,rho_PbS,rho_OA,rho_head,sld_j*ones(size(xdata')),mean_q_data,sig_q_data,shadow_data,x_gh,w_gh);
    plot(xdata',ymedian,'-b','LineWidth',2)
    yMLE = call_model(xdata',xMLE,rho_PbS,rho_OA,rho_head,sld_j*ones(size(xdata')),mean_q_data,sig_q_data,shadow_data,x_gh,w_gh);
    plot(xdata',yMLE,'--r','LineWidth',2)
end
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'TickLength',[0.02 0.02])
set(gca,'LineWidth',1.5)
xlim([0.004 1])
ylim([0.004 110])
xlabel('q (Å^{-1})')
xticks([0.01 0.1 1])
xticklabels({'10^{-2}','10^{-1}','10^0'})
yticks([0.01 0.1 1 10])
yticklabels({'10^{-2}','10^{-1}','10^0','10^1'})
ylabel('I(q) (cm^{-1})')
title(titlespec)

%% SLD Profile

rvec_core = func_struct.rvec_core;
rvec_lig = func_struct.rvec_lig;
rvec_solvent = func_struct.rvec_solvent;
phivec_core = func_struct.phivec_core;
phivec_lig = func_struct.phivec_lig;
phivec_solvent = func_struct.phivec_solvent;
rvecfit = func_struct.rvecfit;
SLDfit = func_struct.SLDfit;
fhead = func_struct.fhead_eq;

figure()
for i = 1:Nwalkers
    xwalk = xfull(:,i,end);
    plot(rvec_core(xwalk),phivec_core(xwalk),'-','LineWidth',1.5,'color',[0 0 1 0.15])
    hold on
    plot(rvec_lig(xwalk),phivec_lig(xwalk),'-','LineWidth',1.5,'color',[1 0 0 0.15])
    plot(rvec_solvent(xwalk),phivec_solvent(xwalk),'-','LineWidth',1.5,'color',[1 0 1 0.15])
end
xlabel('Radial Distance (Å)')
ylabel('Volume Fraction (-)')
set(gcf,'position',[692 386 516 327])
set(gca,'FontWeight','normal','LineWidth',2,'FontName','Gill Sans MT','FontSize',14,'TickLength',[0.03 0.03])
xlim([0 round(max(rvecfit(xMLE)),-1)])
ylim([0 1.2])
title(titlespec)

figure()
for i = 1:Nwalkers
    xwalk = xfull(:,i,end);
    plot(rvec_core(xwalk),rho_PbS*phivec_core(xwalk)/1e-6,'-','LineWidth',1.5,'color',[0 0 1 0.15])
    hold on
    plot(rvec_lig(xwalk),[0,rho_head,rho_OA,rho_OA,0].*[0,fhead(xwalk),xwalk(6),xwalk(6),0]/1e-6,'-','LineWidth',1.5,'color',[1 0 0 0.15])
    plot(rvec_solvent(xwalk),rho_d8*phivec_solvent(xwalk)/1e-6,'-','LineWidth',1.5,'color',[1 0 1 0.15])
    plot(rvecfit(xwalk),SLDfit(xwalk,rho_PbS,rho_OA,rho_head,rho_d8)/1e-6,'-','LineWidth',1.5,'color',[0.5 0.5 0.5 0.15])
end
plot(rvecfit(xmedian),SLDfit(xmedian,rho_PbS,rho_OA,rho_head,rho_d8)/1e-6,'-k','LineWidth',1.5)
plot(rvecfit(xMLE),SLDfit(xMLE,rho_PbS,rho_OA,rho_head,rho_d8)/1e-6,'--r','LineWidth',1.5)
xlabel('Radial Distance (Å)')
ylabel('SLD (10^{-6} Å^{-2})')
set(gcf,'position',[692 386 516 327])
set(gca,'FontWeight','normal','LineWidth',2,'FontName','Gill Sans MT','FontSize',14,'TickLength',[0.03 0.03])
xlim([0 round(max(rvecfit(xMLE)),-1)])
ylim([0 6])
title(titlespec)

end