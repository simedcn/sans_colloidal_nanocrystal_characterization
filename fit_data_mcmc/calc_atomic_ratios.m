function [atomic_perc,lig_cov,equiv_d,core_ratio] = calc_atomic_ratios(x_mle,func_struct)

rvec_core = func_struct.rvec_core;
rvec_lig = func_struct.rvec_lig;
phivec_core = func_struct.phivec_core;
phivec_lig = func_struct.phivec_lig;

% Calculate atomic ratios from fit volume fraction profiles
dens_PbS = 7.6; % [g/cm^3]
dens_OA = 0.895; % [g/cm^3]
MW_PbS = 239.3; % [g/mol]
MW_OA = 282.47; % [g/mol]
NA = 6.022e23;

totalmass = @(rvec,phivec,dens_bulk) (1e-8)^3*dens_bulk*int_lin_profile(rvec,phivec);
mass_core = totalmass(rvec_core(x_mle),phivec_core(x_mle),dens_PbS); % [g]
mass_lig = totalmass(rvec_lig(x_mle),phivec_lig(x_mle),dens_OA); % [g]

N_PbS_core = mass_core*NA/MW_PbS;
N_OA_shell = mass_lig*NA/MW_OA;
N_Pb_core = N_PbS_core + (1/4)*N_OA_shell; % Use ligand charge balance to get correct Pb/S ratio in the core
N_S_core = N_Pb_core - (1/2)*N_OA_shell;
core_ratio = N_Pb_core/N_S_core;


Ntotal = N_Pb_core + N_S_core;
atomic_perc = 100*[(N_Pb_core)/Ntotal,(N_S_core)/Ntotal]';
core_volume = (1e7)^3*(mass_core/dens_PbS); % [nm^3]
PbS_volume = (1e7)^3*(mass_core/dens_PbS);
equiv_d = 2*((3/4/pi)*PbS_volume)^(1/3);
equiv_R_for_lig_cov = ((3/4/pi)*core_volume)^(1/3);
surf_area_sphere = 4*pi*(equiv_R_for_lig_cov)^2;
max_sphericity = 0.9212;
surf_area = surf_area_sphere/max_sphericity;
lig_cov = N_OA_shell/surf_area; % [ligands/nm^2]

end

