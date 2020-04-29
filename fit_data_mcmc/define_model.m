function [call_model,func_struct] = define_model()

% Base distances
rhalf = @(rmin,tvertex) rmin + tvertex/2;
rvertex = @(rmin,tvertex) rmin + tvertex;
rplateau = @(rmin,tvertex,tplateau) rmin + tvertex + tplateau;
rmax = @(rmin,tvertex,tplateau,tout) rmin + tvertex + tplateau + tout;

% Distances of each kink
r1 = @(rmin) 0;
r2 = @(rmin) rmin;
r3 = @(rmin,tvertex) rhalf(rmin,tvertex);
r4 = @(rmin,tvertex) rvertex(rmin,tvertex);
r5 = @(rmin,tvertex,tplateau) rplateau(rmin,tvertex,tplateau);
r6 = @(rmin,tvertex,tplateau,tout) rmax(rmin,tvertex,tplateau,tout);

% Needed for ligand consistency
NA = 6.022e23;
dens_OA = 0.895; % [g/cm^3]
dens_head = 1.56; % [g/cm^3]
MW_OA = 282.47; % [g/mol]
MW_head = 44.01; % [g/mol]
rveclig = @(rmin,tvertex,tplateau,tout) [r2(rmin),r3(rmin,tvertex),r4(rmin,tvertex),r5(rmin,tvertex,tplateau),r6(rmin,tvertex,tplateau,tout)];
phiveclig = @(fhead,fOA) [0,fhead,fOA,fOA,0];
m_lig = @(rmin,tvertex,tplateau,tout,fhead,fOA) (1e-8)^3*dens_OA*int_lin_profile(rveclig(rmin,tvertex,tplateau,tout),phiveclig(fhead,fOA)); % [g]
N_lig = @(rmin,tvertex,tplateau,tout,fhead,fOA) m_lig(rmin,tvertex,tplateau,tout,fhead,fOA)*NA/MW_OA; % [molecules]
m_head = @(rmin,tvertex,fhead) (1e-8)^3*dens_head*int_lin_profile([r2(rmin),r3(rmin,tvertex),r4(rmin,tvertex)],[0,fhead,0]); % [g]
N_head = @(rmin,tvertex,fhead) m_head(rmin,tvertex,fhead)*NA/MW_head; % [molecules]

% Solve for f_head
tol = 1e-4;
options = optimset('TolX',tol,'TolFun',tol);
fhead_eq = @(rmin,tvertex,tplateau,tout,fOA) fzero(@(fhead) N_lig(rmin,tvertex,tplateau,tout,fhead,fOA) - N_head(rmin,tvertex,fhead),0.5,options);

% Base volume fraction values
phi_core3 = @(fPbS) fPbS;
phi_lig3 = @(rmin,tvertex,tplateau,tout,fOA) fhead_eq(rmin,tvertex,tplateau,tout,fOA);
phi_solv3 = @(rmin,tvertex,tplateau,tout,fPbS,fOA) 1 - fPbS - fhead_eq(rmin,tvertex,tplateau,tout,fOA);
phi_lig4 = @(fOA) fOA;
phi_lig5 = @(fOA) fOA;
phi_solv4 = @(fOA) 1 - fOA;
phi_solv5 = @(fOA) 1 - fOA;

% SLD values at kinks
rho1 = @(rho_PbS) rho_PbS;
rho2 = @(rho_PbS) rho_PbS;
rho3 = @(rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_solv) phi_core3(fPbS)*rho_PbS + phi_lig3(rmin,tvertex,tplateau,tout,fOA)*rho_head + phi_solv3(rmin,tvertex,tplateau,tout,fPbS,fOA)*rho_solv;
rho4 = @(fOA,rho_OA,rho_solv) phi_lig4(fOA)*rho_OA + phi_solv4(fOA)*rho_solv;
rho5 = @(fOA,rho_OA,rho_solv) phi_lig5(fOA)*rho_OA + phi_solv5(fOA)*rho_solv;
rho6 = @(rho_solv) rho_solv;

% Form factor equations
y0 = @(q,r) -1*cos(q.*r)./(q.*r);
y1 = @(q,r) -1*cos(q.*r)./(q.*r).^2 - sin(q.*r)./(q.*r);
j1 = @(q,r) (sin(q.*r) - q.*r.*cos(q.*r))./(q.*r).^2;
alpha_in = @(rin,rout) rin./(rout - rin);
alpha_out = @(rin,rout) rout./(rout - rin);
F_flat = @(q,rin,rout,rho_in,rho_out,rho_solv) 4*pi*((rho_in - rho_solv).*(rout.^3.*j1(q,rout)./(q.*rout)));
F_lin = @(q,rin,rout,rho_in,rho_out,rho_solv) 4*pi*((rho_in - rho_solv).*(rout.^3.*j1(q,rout)./(q.*rout) - rin.^3.*j1(q,rin)./(q.*rin)) + ...
    (rho_out - rho_in).*rout.^3.*j1(q,rout)./(q.*rout) + ...
    (rho_out - rho_in).*alpha_out(rin,rout).*rout.^3.*(-1*y1(q,rout)./(q.*rout).^2 - y0(q,rout)./(q.*rout).^3) - ...
    (rho_out - rho_in).*alpha_in(rin,rout).*rin.^3.*(-1*y1(q,rin)./(q.*rin).^2 - y0(q,rin)./(q.*rin).^3));

% Form factor for each linear section
F1 = @(q,rmin,rho_PbS,rho_solv) F_flat(q,r1(rmin),r2(rmin),rho1(rho_PbS),rho2(rho_PbS),rho_solv); % 1 to 2
F2 = @(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_solv) F_lin(q,r2(rmin),r3(rmin,tvertex),rho2(rho_PbS),rho3(rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_solv),rho_solv);
F3 = @(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_OA,rho_solv) F_lin(q,r3(rmin,tvertex),r4(rmin,tvertex),rho3(rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_solv),rho4(fOA,rho_OA,rho_solv),rho_solv);
F4 = @(q,rmin,tvertex,tplateau,fOA,rho_OA,rho_solv) F_lin(q,r4(rmin,tvertex),r5(rmin,tvertex,tplateau),rho4(fOA,rho_OA,rho_solv),rho5(fOA,rho_OA,rho_solv),rho_solv);
F5 = @(q,rmin,tvertex,tplateau,tout,fOA,rho_OA,rho_solv) F_lin(q,r5(rmin,tvertex,tplateau),r6(rmin,tvertex,tplateau,tout),rho5(fOA,rho_OA,rho_solv),rho6(rho_solv),rho_solv);

% Composite Intensity
F = @(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_OA,rho_solv) ...
    F1(q,rmin,rho_PbS,rho_solv) + ...
    F2(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_solv) + ...
    F3(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_OA,rho_solv) + ...
    F4(q,rmin,tvertex,tplateau,fOA,rho_OA,rho_solv) + ...
    F5(q,rmin,tvertex,tplateau,tout,fOA,rho_OA,rho_solv);
P = @(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_OA,rho_solv) F(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_OA,rho_solv).^2;
model = @(q,rmin,tvertex,tplateau,tout,fPbS,fOA,phi,bkgd,rho_PbS,rho_head,rho_OA,rho_solv) (1e8).*(phi./((4*pi/3).*rmax(rmin,tvertex,tplateau,tout).^3)).*P(q,rmin,tvertex,tplateau,tout,fPbS,fOA,rho_PbS,rho_head,rho_OA,rho_solv) + bkgd;
model_smear = @(q,rmin,tvertex,tplateau,tout,fPbS,fOA,phi,bkgd,rho_PbS,rho_head,rho_OA,rho_solv,mean_q,sig_q,shadow,x_gh,w_gh) ...
    1/(sqrt(pi)) * sum( shadow.*( (mean_q + sqrt(2)*sig_q.*x_gh) >= 0).*model(mean_q + sqrt(2)*sig_q.*x_gh,rmin,tvertex,tplateau,tout,fPbS,fOA,phi,bkgd,rho_PbS,rho_head,rho_OA,rho_solv).*w_gh,1 );

call_model = @(q,x,rho_PbS,rho_OA,rho_head,rho_solv,mean_q,sig_q,shadow,x_gh,w_gh) model_smear(q,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),rho_PbS,rho_head,rho_OA,rho_solv,mean_q,sig_q,shadow,x_gh,w_gh);

%% Call functions of the set of parameters x. Return a data structure of function handles.
rvec_core = @(x) [r1(x(1)),r2(x(1)),r3(x(1),x(2)),r4(x(1),x(2))];
rvec_lig = @(x) [r2(x(1)),r3(x(1),x(2)),r4(x(1),x(2)),r5(x(1),x(2),x(3)),r6(x(1),x(2),x(3),x(4))];
rvec_solvent = @(x) [r2(x(1)),r3(x(1),x(2)),r4(x(1),x(2)),r5(x(1),x(2),x(3)),r6(x(1),x(2),x(3),x(4)),r6(x(1),x(2),x(3),x(4))+10];
phivec_core = @(x) [1,1,x(5),0];
phivec_lig = @(x) [0,phi_lig3(x(1),x(2),x(3),x(4),x(6)),phi_lig4(x(6)),phi_lig5(x(6)),0];
phivec_solvent = @(x) [0,phi_solv3(x(1),x(2),x(3),x(4),x(5),x(6)),phi_solv4(x(6)),phi_solv5(x(6)),1,1];

rvecfit = @(x) [r1(x(1)),r2(x(1)),r3(x(1),x(2)),r4(x(1),x(2)),r5(x(1),x(2),x(3)),r6(x(1),x(2),x(3),x(4)),r6(x(1),x(2),x(3),x(4)) + 10];
SLDfit = @(x,rho_PbS,rho_OA,rho_head,rho_solv) rho_PbS*[phivec_core(x),0,0,0] + [0,0,rho_head,rho_OA,rho_OA,0,0].*[0,0,phi_lig3(x(1),x(2),x(3),x(4),x(6)),phi_lig4(x(6)),phi_lig5(x(6)),0,0] + rho_solv*[0,phivec_solvent(x)];

call_fhead_eq = @(x) fhead_eq(x(1),x(2),x(3),x(4),x(6));

func_struct = struct();
func_struct.rvec_core = rvec_core;
func_struct.rvec_lig = rvec_lig;
func_struct.rvec_solvent = rvec_solvent;
func_struct.phivec_core = phivec_core;
func_struct.phivec_lig = phivec_lig;
func_struct.phivec_solvent = phivec_solvent;
func_struct.rvecfit = rvecfit;
func_struct.SLDfit = SLDfit;
func_struct.fhead_eq = call_fhead_eq;

end

