%% Inversi_TETM_DirectJ_Broyden_clean_plots_v2.m
% Clean English release script for the representative 2-D MT joint TE-TM inversion demo.
% This script loads fixed synthetic data and cached direct-sensitivity Jacobians,
% runs the controlled Broyden inversion, prints a compact Command Window summary,
% and displays model/RMS/data-fit figures.
%
% Paper-specific text reports, CSV summaries, helper-code copying, and report exports
% are intentionally omitted from this release version.
%
clc;
clearvars;
close all;
t_total_start = tic;

% Release-path setup: make the script independent of MATLAB Current Folder.
this_file = mfilename('fullpath');
if isempty(this_file)
    repo_root = pwd;
else
    repo_root = fileparts(fileparts(this_file));  % this file is expected in repo_root/src
end


%% ===================== USER SETTINGS =====================
use_topo            = true;
time_sign           = +1;

% Shared earth model.
% TE_NDZ includes air and earth cells. With the 1:3 air-to-earth split,
% the TE earth part is TE_NDZb = 38, equal to TM_NDZ.
TE_NDZ              = 51;
TM_NDZ              = 38;
TE_ZMIN             = -50000;
TM_ZMIN             = -50000;

use_phase_in_radian = true;
param_mode          = 'cell';
include_air_tm      = false;

phase_mode_te       = 'raw';
phase_mode_tm       = 'raw';
phase_plot_mode     = 'inversion';   % 'inversion' | 'raw_display'

% Select one model type:
% model_type = 'homogen';
% model_type = 'layered_10_100';
% model_type = 'layered_100_10';
% model_type = 'vertical_10_100';
% model_type = 'vertical_100_10';
model_type = 'layered_100_10';

mode_name = 'Joint_TE_TM';
run_tag   = 'DirectJ_Broyden';

% Controlled Broyden update settings.
% The initial Jacobian is loaded from the cached direct-sensitivity files.
use_broyden         = true;
broyden_start_outer = 1;
broyden_start_rms   = inf;
broyden_sec_err_max = 0.15;

lambda_max          = 1e12;
tol_rms_joint       = 1.0;
tol_outer_improve   = 7e-4;
max_ls              = 16;
print_every         = 10;

% Clean release options
show_figures        = true;    % display model, RMS, and data-fit figures
save_results        = true;   % save only a compact .mat result if true
save_figures        = false;   % save figures only if true

% Direct-sensitivity Jacobian settings
% For the first run of a new model, set force_build_initial_jacobian = true.
% After a valid direct-J file has been generated, cached Jacobians can be loaded.
% Set force_build_initial_jacobian = false to load the cached Jacobian directly.

force_rebuild_each_outer     = false;
force_build_initial_jacobian = false;
load_only_jacobian           = true;
validate_direct_jacobian = false;   % true only for checking selected columns against finite differences
n_validate_cols              = 6;
rhs_chunk_size               = 100;

%% ===================== JOINT TE-TM PRESETS BY MODEL =====================
switch lower(strtrim(model_type))

    case 'homogen'
        case_name   = 'synthetic_JointTETM_homogeneous_DirectJ_Broyden';
        save_prefix = 'TETM_homogeneous_DirectJ_Broyden';
        model_label = 'Homogeneous';

        n_outer           = 8;
        inner_iter        = 10;

        lambda0           = 3.0e4;
        lambda_min        = 25;
        step_max          = 0.025;

        alpha_smallness   = 2.0e-2;
        wx_reg            = 3.5;
        wz_reg            = 3.5;
        reg_eps           = 2.0e-3;

        mode_weight_te    = 1.0;
        mode_weight_tm    = 1.0;

    case 'layered_10_100'
        case_name   = 'synthetic_JointTETM_layer_10_100_DirectJ_Broyden';
        save_prefix = 'TETM_layer_10_100_DirectJ_Broyden';
        model_label = 'Layered 10-100 \Omega m';

        n_outer           = 10;
        inner_iter        = 10;

        lambda0           = 2.5e4;
        lambda_min        = 8.0;
        step_max          = 0.030;

        alpha_smallness   = 1.2e-2;
        wx_reg            = 3.4;
        wz_reg            = 1.25;
        reg_eps           = 1.2e-3;

        mode_weight_te    = 0.85;
        mode_weight_tm    = 1.00;

    case 'layered_100_10'
        case_name   = 'synthetic_JointTETM_layer_100_10_DirectJ_Broyden';
        save_prefix = 'TETM_layer_100_10_DirectJ_Broyden';
        model_label = 'Layered 100-10 \Omega m';

        n_outer           = 10;
        inner_iter        = 10;

        lambda0           = 2.8e4;
        lambda_min        = 10.0;
        step_max          = 0.028;

        alpha_smallness   = 1.6e-2;
        wx_reg            = 3.2;
        wz_reg            = 1.35;
        reg_eps           = 1.5e-3;

        mode_weight_te    = 0.90;
        mode_weight_tm    = 1.00;

    case 'vertical_10_100'
        case_name   = 'synthetic_JointTETM_vertical_10_100_DirectJ_Broyden';
        save_prefix = 'TETM_vertical_10_100_DirectJ_Broyden';
        model_label = 'Vertical contact 10-100 \Omega m';

        n_outer           = 10;
        inner_iter        = 10;

        lambda0           = 1.2e4;
        lambda_min        = 3.5;
        step_max          = 0.040;

        alpha_smallness   = 6.0e-3;
        wx_reg            = 0.70;
        wz_reg            = 2.60;
        reg_eps           = 8.0e-4;

        mode_weight_te    = 0.90;
        mode_weight_tm    = 1.00;

    case 'vertical_100_10'
        case_name   = 'synthetic_JointTETM_vertical_100_10_DirectJ_Broyden';
        save_prefix = 'TETM_vertical_100_10_DirectJ_Broyden';
        model_label = 'Vertical contact 100-10 \Omega m';

        n_outer           = 10;
        inner_iter        = 10;

        lambda0           = 1.5e4;
        lambda_min        = 4.0;
        step_max          = 0.038;

        alpha_smallness   = 7.0e-3;
        wx_reg            = 0.80;
        wz_reg            = 3.00;
        reg_eps           = 1.0e-3;

        mode_weight_te    = 0.90;
        mode_weight_tm    = 1.00;

    otherwise
        error('Unknown model_type: %s', model_type);
end

%% ===================== OUTPUT FOLDER =====================
output_dir = fullfile(repo_root, ['outputs_' case_name]);
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end
fprintf('Output folder: %s\n', output_dir);
fprintf('Run tag: %s | use_broyden = %d\n', run_tag, use_broyden);

%% ===================== FILE SETTINGS =====================
switch lower(strtrim(model_type))
    case 'homogen'
        input_case_dir = 'ascii_homogeneous';

    case 'layered_100_10'
        input_case_dir = 'ascii_layered_100_10';

    case 'layered_10_100'
        input_case_dir = 'ascii_layered_10_100';

    case 'vertical_100_10'
        input_case_dir = 'ascii_vertical_100_10';

    case 'vertical_10_100'
        input_case_dir = 'ascii_vertical_10_100';

    otherwise
        error('Unknown model_type: %s', model_type);
end

input_dir = fullfile(repo_root, 'data', input_case_dir);
jac_dir   = fullfile(repo_root, 'data', 'jacobian');

if exist(input_dir, 'dir') ~= 7
    error('Input folder was not found: %s', input_dir);
end

if exist(jac_dir, 'dir') ~= 7
    error('Jacobian folder was not found: %s', jac_dir);
end

if use_topo
    file_topo = pick_first_existing({ ...
        fullfile(input_dir,'topomesh_Sint_topo.txt'), ...
        fullfile(input_dir,'topomesh_sint_topo.txt'), ...
        fullfile(input_dir,'topomesh_topo.txt'), ...
        fullfile(repo_root,'topomesh_Sint_topo.txt'), ...
        fullfile(repo_root,'topomesh_topo.txt')});
    file_xst  = pick_first_existing({ ...
        fullfile(input_dir,'x_st_Sint_topo.txt'), ...
        fullfile(input_dir,'x_st_sint_topo.txt'), ...
        fullfile(input_dir,'x_st_topo.txt'), ...
        fullfile(repo_root,'x_st_Sint_topo.txt'), ...
        fullfile(repo_root,'x_st_topo.txt')});
    file_zst  = pick_first_existing({ ...
        fullfile(input_dir,'z_st_Sint_topo.txt'), ...
        fullfile(input_dir,'z_st_sint_topo.txt'), ...
        fullfile(input_dir,'z_st_topo.txt'), ...
        fullfile(repo_root,'z_st_Sint_topo.txt'), ...
        fullfile(repo_root,'z_st_topo.txt')});
else
    file_topo = pick_first_existing({ ...
        fullfile(input_dir,'topomesh_Sint_datar.txt'), ...
        fullfile(input_dir,'topomesh_sint_datar.txt'), ...
        fullfile(input_dir,'topomesh_datar.txt'), ...
        fullfile(repo_root,'topomesh_Sint_datar.txt'), ...
        fullfile(repo_root,'topomesh_datar.txt')});
    file_xst  = pick_first_existing({ ...
        fullfile(input_dir,'x_st_Sint_datar.txt'), ...
        fullfile(input_dir,'x_st_sint_datar.txt'), ...
        fullfile(input_dir,'x_st_datar.txt'), ...
        fullfile(repo_root,'x_st_Sint_datar.txt'), ...
        fullfile(repo_root,'x_st_datar.txt')});
    file_zst  = pick_first_existing({ ...
        fullfile(input_dir,'z_st_Sint_datar.txt'), ...
        fullfile(input_dir,'z_st_sint_datar.txt'), ...
        fullfile(input_dir,'z_st_datar.txt'), ...
        fullfile(repo_root,'z_st_Sint_datar.txt'), ...
        fullfile(repo_root,'z_st_datar.txt')});
end

file_freq = pick_first_existing({ ...
    fullfile(input_dir,'frek_buatan2.txt'), ...
    fullfile(input_dir,'frek.txt'), ...
    fullfile(input_dir,'frek_sintetik.txt'), ...
    fullfile(repo_root,'frek_buatan2.txt'), ...
    fullfile(repo_root,'frek.txt')});

% TE files
file_te_rhoa_obs  = pick_first_existing({fullfile(input_dir,'rhoastTE_obs.txt'), fullfile(repo_root,'rhoastTE_obs.txt')});
file_te_phi_obs   = pick_first_existing({fullfile(input_dir,'fasastTE_obs.txt'), fullfile(input_dir,'fasastTE_obs_deg.txt'), fullfile(repo_root,'fasastTE_obs.txt')});
file_te_sigma_rho = pick_first_existing({fullfile(input_dir,'sigma_logrhoastTE.txt'), fullfile(repo_root,'sigma_logrhoastTE.txt')});
file_te_sigma_phi = pick_first_existing({fullfile(input_dir,'sigma_phistTE_deg.txt'), fullfile(repo_root,'sigma_phistTE_deg.txt')});
file_te_J         = pick_first_existing({ ...
    fullfile(jac_dir, sprintf('J_TE_direct_shared_%s.mat', model_type)), ...
    fullfile(repo_root, sprintf('J_TE_direct_shared_%s.mat', model_type)), ...
    sprintf('J_TE_direct_shared_%s.mat', model_type)});
file_te_mcurrent  = fullfile(repo_root, 'm_current_TE_sint_cell.txt');

% TM files
file_tm_rhoa_obs  = pick_first_existing({fullfile(input_dir,'rhoastTM_obs.txt'), fullfile(repo_root,'rhoastTM_obs.txt')});
file_tm_phi_obs   = pick_first_existing({fullfile(input_dir,'fasastTM_obs.txt'), fullfile(input_dir,'fasastTM_obs_deg.txt'), fullfile(repo_root,'fasastTM_obs.txt')});
file_tm_sigma_rho = pick_first_existing({fullfile(input_dir,'sigma_logrhoastTM.txt'), fullfile(repo_root,'sigma_logrhoastTM.txt')});
file_tm_sigma_phi = pick_first_existing({fullfile(input_dir,'sigma_phistTM_deg.txt'), fullfile(repo_root,'sigma_phistTM_deg.txt')});
file_tm_J         = pick_first_existing({ ...
    fullfile(jac_dir, sprintf('J_TM_direct_shared_%s.mat', model_type)), ...
    fullfile(repo_root, sprintf('J_TM_direct_shared_%s.mat', model_type)), ...
    sprintf('J_TM_direct_shared_%s.mat', model_type)});
file_tm_mcurrent  = fullfile(repo_root, 'm_current_TM_sint_cell.txt');

fprintf('Input directory : %s\n', input_dir);
fprintf('Jacobian dir    : %s\n', jac_dir);
fprintf('Topography file : %s\n', file_topo);
fprintf('TE Jacobian     : %s\n', file_te_J);
fprintf('TM Jacobian     : %s\n', file_tm_J);

%% ===================== LOAD COMMON =====================
assert(exist(file_topo,'file')==2, 'Topography file was not found.');
assert(exist(file_xst,'file')==2,  'Station x-coordinate file was not found.');
assert(exist(file_freq,'file')==2, 'Frequency file was not found.');

topo  = load(file_topo);
assert(size(topo,2) >= 2, 'topomesh must contain two columns: [x z].');
xtopo = topo(:,1);
ztopo = topo(:,2);

x_st = load(file_xst);
x_st = x_st(:);

if exist(file_zst,'file')==2
    z_st = load(file_zst);
    z_st = z_st(:);
else
    z_st = interp1(xtopo, ztopo, x_st, 'linear', 'extrap');
end

datafreq = load(file_freq);
F = datafreq(:,1);
F = F(:);
NF = numel(F);

NDX  = numel(xtopo)-1;
n_st = numel(x_st);

%% ===================== STATION NAMES =====================
% Station names must follow the column order of the data and x_st.
nm = {'MT-01','MT-02','MT-03','MT-04','MT-05','MT-06', ...
      'MT-07','MT-08','MT-09','MT-10','MT-11'}';

% Example for non-sequential station order:
% nm = {'MT-01','MT-03','MT-07','MT-02','MT-05','MT-11', ...
%       'MT-04','MT-08','MT-10','MT-06','MT-09'}';

if numel(nm) ~= n_st
    error('Number of station names (%d) does not match number of stations (%d).', numel(nm), n_st);
end
nm = nm(:);

%% ===================== CHECK SHARED MODEL SIZE =====================
a = 1;
b = 3;
TE_NDZa = round(a*TE_NDZ/(a+b));
TE_NDZb = TE_NDZ - TE_NDZa;

if TE_NDZb ~= TM_NDZ
    error('TE_NDZb=%d and TM_NDZ=%d must be equal for the shared earth model.', TE_NDZb, TM_NDZ);
end

Nmodel = NDX * TM_NDZ;
fprintf('Shared joint model size = %d (NDX=%d, NDZ=%d)\n', Nmodel, NDX, TM_NDZ);

%% ===================== LOAD TE DATA =====================
assert(exist(file_te_rhoa_obs,'file')==2,  'TE apparent-resistivity observations were not found.');
assert(exist(file_te_phi_obs,'file')==2,   'TE phase observations were not found.');

TE = struct();
TE.rhoa_obs_raw    = load(file_te_rhoa_obs);
TE.phi_obs_deg_raw = load(file_te_phi_obs);

if exist(file_te_sigma_rho,'file')==2
    TE.sigma_log_rhoa = load(file_te_sigma_rho);
else
    TE.sigma_log_rhoa = (0.10/log(10)) * ones(size(TE.rhoa_obs_raw));
end

if exist(file_te_sigma_phi,'file')==2
    TE.sigma_phi_deg = load(file_te_sigma_phi);
else
    TE.sigma_phi_deg = 2.0 * ones(size(TE.phi_obs_deg_raw));
end

assert(size(TE.rhoa_obs_raw,1)==NF, 'TE observations: frequency count mismatch.');
assert(isequal(size(TE.phi_obs_deg_raw), size(TE.rhoa_obs_raw)), 'TE phase observation size mismatch.');
assert(size(TE.rhoa_obs_raw,2)==n_st, 'TE station count mismatch.');

TE.rhoa_obs_raw(TE.rhoa_obs_raw<=0) = realmin;
TE.sigma_log_rhoa(~isfinite(TE.sigma_log_rhoa) | TE.sigma_log_rhoa<=0) = 0.02;
TE.sigma_phi_deg(~isfinite(TE.sigma_phi_deg) | TE.sigma_phi_deg<=0)    = 2.0;

TE.sigma_log_rhoa = max(TE.sigma_log_rhoa, 0.02);
TE.sigma_phi_deg  = max(TE.sigma_phi_deg,  2.0);

TE.rhoa_obs    = TE.rhoa_obs_raw;
TE.phi_obs_deg = apply_phase_mode_deg(TE.phi_obs_deg_raw, phase_mode_te);

if use_phase_in_radian
    TE.phi_obs   = TE.phi_obs_deg*pi/180;
    TE.sigma_phi = max(TE.sigma_phi_deg*pi/180, 1e-6);
else
    TE.phi_obs   = TE.phi_obs_deg;
    TE.sigma_phi = max(TE.sigma_phi_deg, 1e-6);
end

TE.dobs_rho = reshape(log10(TE.rhoa_obs)', [], 1);
TE.dobs_phi = reshape(TE.phi_obs',        [], 1);
TE.sig_rho  = reshape(TE.sigma_log_rhoa', [], 1);
TE.sig_phi  = reshape(TE.sigma_phi',      [], 1);
TE.Ndata    = numel(TE.dobs_rho) + numel(TE.dobs_phi);

%% ===================== LOAD TM DATA =====================
assert(exist(file_tm_rhoa_obs,'file')==2,  'TM apparent-resistivity observations were not found.');
assert(exist(file_tm_phi_obs,'file')==2,   'TM phase observations were not found.');

TM = struct();
TM.rhoa_obs_raw    = load(file_tm_rhoa_obs);
TM.phi_obs_deg_raw = load(file_tm_phi_obs);

if exist(file_tm_sigma_rho,'file')==2
    TM.sigma_log_rhoa = load(file_tm_sigma_rho);
else
    TM.sigma_log_rhoa = (0.10/log(10)) * ones(size(TM.rhoa_obs_raw));
end

if exist(file_tm_sigma_phi,'file')==2
    TM.sigma_phi_deg = load(file_tm_sigma_phi);
else
    TM.sigma_phi_deg = 2.0 * ones(size(TM.phi_obs_deg_raw));
end

assert(size(TM.rhoa_obs_raw,1)==NF, 'TM observations: frequency count mismatch.');
assert(isequal(size(TM.phi_obs_deg_raw), size(TM.rhoa_obs_raw)), 'TM phase observation size mismatch.');
assert(size(TM.rhoa_obs_raw,2)==n_st, 'TM station count mismatch.');

TM.rhoa_obs_raw(TM.rhoa_obs_raw<=0) = realmin;
TM.sigma_log_rhoa(~isfinite(TM.sigma_log_rhoa) | TM.sigma_log_rhoa<=0) = 0.02;
TM.sigma_phi_deg(~isfinite(TM.sigma_phi_deg) | TM.sigma_phi_deg<=0)    = 2.0;

TM.sigma_log_rhoa = max(TM.sigma_log_rhoa, 0.02);
TM.sigma_phi_deg  = max(TM.sigma_phi_deg,  2.0);

TM.rhoa_obs    = TM.rhoa_obs_raw;
TM.phi_obs_deg = apply_phase_mode_deg(TM.phi_obs_deg_raw, phase_mode_tm);

if use_phase_in_radian
    TM.phi_obs   = TM.phi_obs_deg*pi/180;
    TM.sigma_phi = max(TM.sigma_phi_deg*pi/180, 1e-6);
else
    TM.phi_obs   = TM.phi_obs_deg;
    TM.sigma_phi = max(TM.sigma_phi_deg, 1e-6);
end

TM.dobs_rho = reshape(log10(TM.rhoa_obs)', [], 1);
TM.dobs_phi = reshape(TM.phi_obs',        [], 1);
TM.sig_rho  = reshape(TM.sigma_log_rhoa', [], 1);
TM.sig_phi  = reshape(TM.sigma_phi',      [], 1);
TM.Ndata    = numel(TM.dobs_rho) + numel(TM.dobs_phi);

Ndata_joint = TE.Ndata + TM.Ndata;

%% ===================== INITIAL MODEL =====================
m     = ones(Nmodel,1)*log10(31.6);
m_ref = m;

% Bounds used for the 1-1000 Ohm m synthetic benchmark.
% They prevent unrealistically large values produced only to reduce RMS.
m_min = 0.0;   % log10(1 Ohm m)
m_max = 3.0;   % log10(1000 Ohm m)

%% ===================== LOAD / BUILD DIRECT-SENSITIVITY JACOBIANS =====================
[J_te, J_te_base, built_te0, metaJ_te, elapsed_jacobian_te_s] = load_or_build_direct_jacobian( ...
    file_te_J, m, TE.Ndata, Nmodel, load_only_jacobian, force_build_initial_jacobian, 'TE', ...
    xtopo, ztopo, datafreq, x_st, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
    time_sign, use_topo, include_air_tm, use_phase_in_radian, param_mode, ...
    model_type, rhs_chunk_size, validate_direct_jacobian, n_validate_cols);

[J_tm, J_tm_base, built_tm0, metaJ_tm, elapsed_jacobian_tm_s] = load_or_build_direct_jacobian( ...
    file_tm_J, m, TM.Ndata, Nmodel, load_only_jacobian, force_build_initial_jacobian, 'TM', ...
    xtopo, ztopo, datafreq, x_st, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
    time_sign, use_topo, include_air_tm, use_phase_in_radian, param_mode, ...
    model_type, rhs_chunk_size, validate_direct_jacobian, n_validate_cols);

elapsed_jacobian_total_s = elapsed_jacobian_te_s + elapsed_jacobian_tm_s;
W_te = spdiags(1./[TE.sig_rho; TE.sig_phi], 0, TE.Ndata, TE.Ndata);
W_tm = spdiags(1./[TM.sig_rho; TM.sig_phi], 0, TM.Ndata, TM.Ndata);

Jw_te_base = sqrt(mode_weight_te) * (W_te * J_te);
Jw_tm_base = sqrt(mode_weight_tm) * (W_tm * J_tm);

Jw_te = Jw_te_base;
Jw_tm = Jw_tm_base;

%% ===================== REGULARIZATION =====================
R = inv_occam_TM2(NDX, Nmodel, 'param_mode', param_mode, 'wx', wx_reg, 'wz', wz_reg);
R = sparse(double(R));
Ireg = speye(Nmodel);
RtR  = (R'*R) + (alpha_smallness + reg_eps)*Ireg;

%% ===================== HISTORIES =====================
best_obj = inf;
best_rms_joint = inf;
best_m   = m;

outer_rms_joint_hist = nan(n_outer,1);
outer_rms_te_hist    = nan(n_outer,1);
outer_rms_tm_hist    = nan(n_outer,1);

inner_rms_joint_hist = nan(n_outer, inner_iter);
inner_rms_te_hist    = nan(n_outer, inner_iter);
inner_rms_tm_hist    = nan(n_outer, inner_iter);
inner_lambda_hist    = nan(n_outer, inner_iter);

rms_joint_hist = [];
rms_te_hist    = [];
rms_tm_hist    = [];
lambda_hist    = [];

lambda = lambda0;
stop_all = false;

%% ===================== RUN COUNTERS =====================
n_step_accepted         = 0;
n_step_rejected         = 0;
n_broyden_update_te     = 0;
n_broyden_update_tm     = 0;
n_broyden_rejected_te   = 0;
n_broyden_rejected_tm   = 0;
n_reset_to_baseJ_te     = 0;
n_reset_to_baseJ_tm     = 0;
n_direct_jacobian_build_te = built_te0;
n_direct_jacobian_build_tm = built_tm0;

% Legacy counters are retained for backward compatibility;
% their values now refer to full direct-sensitivity Jacobian builds.
n_full_jacobian_build_te = built_te0;
n_full_jacobian_build_tm = built_tm0;

sec_err_te_hist = [];
sec_err_tm_hist = [];
alpha_hist      = [];
obj_hist        = [];

rms_te_unweighted_hist = [];
rms_tm_unweighted_hist = [];

t_inversion_start = tic;

%% ===================== OUTER LOOP =====================
for io = 1:n_outer
    fprintf('\n================ OUTER %d / %d ================\n', io, n_outer);

    if io > 1 && force_rebuild_each_outer
        if load_only_jacobian
            warning('force_rebuild_each_outer is ignored because load_only_jacobian=true.');
        else
            [J_te, J_te_base, built_te, metaJ_te, elapsed_te_tmp] = load_or_build_direct_jacobian( ...
                file_te_J, m, TE.Ndata, Nmodel, false, true, 'TE', ...
                xtopo, ztopo, datafreq, x_st, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
                time_sign, use_topo, include_air_tm, use_phase_in_radian, param_mode, ...
                model_type, rhs_chunk_size, false, n_validate_cols);
            [J_tm, J_tm_base, built_tm, metaJ_tm, elapsed_tm_tmp] = load_or_build_direct_jacobian( ...
                file_tm_J, m, TM.Ndata, Nmodel, false, true, 'TM', ...
                xtopo, ztopo, datafreq, x_st, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
                time_sign, use_topo, include_air_tm, use_phase_in_radian, param_mode, ...
                model_type, rhs_chunk_size, false, n_validate_cols);

            elapsed_jacobian_te_s = elapsed_jacobian_te_s + elapsed_te_tmp;
            elapsed_jacobian_tm_s = elapsed_jacobian_tm_s + elapsed_tm_tmp;
            elapsed_jacobian_total_s = elapsed_jacobian_te_s + elapsed_jacobian_tm_s;

            n_direct_jacobian_build_te = n_direct_jacobian_build_te + built_te;
            n_direct_jacobian_build_tm = n_direct_jacobian_build_tm + built_tm;
            n_full_jacobian_build_te = n_full_jacobian_build_te + built_te;
            n_full_jacobian_build_tm = n_full_jacobian_build_tm + built_tm;

            Jw_te_base = sqrt(mode_weight_te) * (W_te * J_te);
            Jw_tm_base = sqrt(mode_weight_tm) * (W_tm * J_tm);
            Jw_te = Jw_te_base;
            Jw_tm = Jw_tm_base;
        end
    end

    for ii = 1:inner_iter
        [r_joint, r_te, r_tm] = residual_joint_weighted( ...
            m, xtopo, ztopo, datafreq, x_st, ...
            TE, TM, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
            time_sign, use_topo, include_air_tm, use_phase_in_radian, ...
            phase_mode_te, phase_mode_tm, param_mode, ...
            mode_weight_te, mode_weight_tm);

        obj_d = 0.5*(r_joint' * r_joint);
        dmref = (m - m_ref);
        rough = R*dmref;
        small = sqrt(alpha_smallness + reg_eps)*dmref;
        obj_m = 0.5*lambda*(rough'*rough + small'*small);
        obj   = obj_d + obj_m;

        rms_joint = sqrt(mean(r_joint.^2));
        rms_te    = sqrt(mean(r_te.^2));
        rms_tm    = sqrt(mean(r_tm.^2));

        inner_rms_joint_hist(io,ii) = rms_joint;
        inner_rms_te_hist(io,ii)    = rms_te;
        inner_rms_tm_hist(io,ii)    = rms_tm;
        inner_lambda_hist(io,ii)    = lambda;

        rms_joint_hist(end+1,1) = rms_joint;
        rms_te_hist(end+1,1)    = rms_te;
        rms_tm_hist(end+1,1)    = rms_tm;
        lambda_hist(end+1,1)    = lambda;

        if obj < best_obj
            best_obj = obj;
            best_rms_joint = rms_joint;
            best_m = m;
        end

        if mod(ii,print_every)==1 || ii==1 || ii==inner_iter
            fprintf('Outer %d | Inner %2d | RMSjoint=%.6f | RMSTE=%.6f | RMSTM=%.6f | Obj=%.6e | lambda=%.3e\n', ...
                io, ii, rms_joint, rms_te, rms_tm, obj, lambda);
        end

        Jw_joint = [Jw_te; Jw_tm];
        g = Jw_joint.' * r_joint + lambda*(RtR*(m-m_ref));
        A = Jw_joint.' * Jw_joint + lambda*RtR;
        dm = -A\g;

        if any(~isfinite(dm))
            warning('dm contains NaN/Inf. Increasing lambda and resetting the Jacobian.');
            lambda = min(lambda*10, lambda_max);
            Jw_te = Jw_te_base;
            Jw_tm = Jw_tm_base;
            n_reset_to_baseJ_te = n_reset_to_baseJ_te + 1;
            n_reset_to_baseJ_tm = n_reset_to_baseJ_tm + 1;
            continue
        end

        maxabs = max(abs(dm));
        if maxabs > step_max
            dm = dm * (step_max/maxabs);
        end

        accepted = false;
        alpha = 1.0;

        for ils = 1:max_ls
            m_try = min(max(m + alpha*dm, m_min), m_max);

            [r_joint_try, ~, ~] = residual_joint_weighted( ...
                m_try, xtopo, ztopo, datafreq, x_st, ...
                TE, TM, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
                time_sign, use_topo, include_air_tm, use_phase_in_radian, ...
                phase_mode_te, phase_mode_tm, param_mode, ...
                mode_weight_te, mode_weight_tm);

            obj_d_try = 0.5*(r_joint_try' * r_joint_try);
            dmref_try = (m_try - m_ref);
            rough_try = R*dmref_try;
            small_try = sqrt(alpha_smallness + reg_eps)*dmref_try;
            obj_m_try = 0.5*lambda*(rough_try'*rough_try + small_try'*small_try);
            obj_try   = obj_d_try + obj_m_try;

            if obj_try < obj
                accepted = true;
                break
            else
                alpha = alpha/2;
            end
        end

        if accepted
            n_step_accepted = n_step_accepted + 1;
            alpha_hist(end+1,1) = alpha;
            obj_hist(end+1,1)   = obj_try;

            dm_acc = m_try - m;
            dr_acc = r_joint_try - r_joint;
            m = m_try;

            if obj_try < best_obj
                best_obj = obj_try;
                best_rms_joint = sqrt(mean(r_joint_try.^2));
                best_m = m_try;
            end

            if use_broyden && io >= broyden_start_outer && rms_joint < broyden_start_rms
                dr_te = dr_acc(1:numel(r_te));
                dr_tm = dr_acc(numel(r_te)+1:end);

                denom = dm_acc'*dm_acc + eps;

                sec_te = norm(dr_te - Jw_te*dm_acc) / max(norm(dr_te), eps);
                sec_tm = norm(dr_tm - Jw_tm*dm_acc) / max(norm(dr_tm), eps);

                sec_err_te_hist(end+1,1) = sec_te;
                sec_err_tm_hist(end+1,1) = sec_tm;

                if sec_te <= broyden_sec_err_max
                    Jw_te = Jw_te + ((dr_te - Jw_te*dm_acc) * dm_acc') / denom;
                    n_broyden_update_te = n_broyden_update_te + 1;
                else
                    Jw_te = Jw_te_base;
                    n_broyden_rejected_te = n_broyden_rejected_te + 1;
                    n_reset_to_baseJ_te = n_reset_to_baseJ_te + 1;
                end

                if sec_tm <= broyden_sec_err_max
                    Jw_tm = Jw_tm + ((dr_tm - Jw_tm*dm_acc) * dm_acc') / denom;
                    n_broyden_update_tm = n_broyden_update_tm + 1;
                else
                    Jw_tm = Jw_tm_base;
                    n_broyden_rejected_tm = n_broyden_rejected_tm + 1;
                    n_reset_to_baseJ_tm = n_reset_to_baseJ_tm + 1;
                end
            end

            lambda = max(lambda/1.5, lambda_min);

            if sqrt(mean(r_joint_try.^2)) < tol_rms_joint
                fprintf('Stop: joint RMS is below tolerance.\n');
                stop_all = true;
                break
            end
        else
            n_step_rejected = n_step_rejected + 1;
            lambda = min(lambda*3, lambda_max);
            Jw_te = Jw_te_base;
            Jw_tm = Jw_tm_base;
            n_reset_to_baseJ_te = n_reset_to_baseJ_te + 1;
            n_reset_to_baseJ_tm = n_reset_to_baseJ_tm + 1;
            fprintf('  Step rejected -> lambda=%.3e\n', lambda);
        end
    end

    good_idx = isfinite(inner_rms_joint_hist(io,:));
    if any(good_idx)
        outer_rms_joint_hist(io) = min(inner_rms_joint_hist(io,good_idx));
        outer_rms_te_hist(io)    = min(inner_rms_te_hist(io,good_idx));
        outer_rms_tm_hist(io)    = min(inner_rms_tm_hist(io,good_idx));
    end

    if io > 1 && isfinite(outer_rms_joint_hist(io-1)) && isfinite(outer_rms_joint_hist(io))
        improve = outer_rms_joint_hist(io-1) - outer_rms_joint_hist(io);
        if improve < tol_outer_improve
            fprintf('Stop outer loop: improvement is below tolerance.\n');
            break
        end
    end

    if stop_all
        break
    end
end

fprintf('\nBest joint RMS = %.6f | Best Obj = %.6e\n', best_rms_joint, best_obj);
m = best_m;

% Recompute RMS using the final/best model that is actually used.
[r_final_joint_vec, r_final_te_vec, r_final_tm_vec] = residual_joint_weighted( ...
    m, xtopo, ztopo, datafreq, x_st, ...
    TE, TM, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
    time_sign, use_topo, include_air_tm, use_phase_in_radian, ...
    phase_mode_te, phase_mode_tm, param_mode, ...
    mode_weight_te, mode_weight_tm);
rms_final_joint_model = sqrt(mean(r_final_joint_vec.^2));
rms_final_te_model    = sqrt(mean(r_final_te_vec.^2));
rms_final_tm_model    = sqrt(mean(r_final_tm_vec.^2));
fprintf('Final saved-model RMS: Joint=%.6f | TE=%.6f | TM=%.6f\n', ...
    rms_final_joint_model, rms_final_te_model, rms_final_tm_model);

%% ===================== RUN SUMMARY =====================
elapsed_inversion_s = toc(t_inversion_start);
elapsed_total_s = toc(t_total_start);

lambda_hist    = lambda_hist(isfinite(lambda_hist));
rms_joint_hist = rms_joint_hist(isfinite(rms_joint_hist));
rms_te_hist    = rms_te_hist(isfinite(rms_te_hist));
rms_tm_hist    = rms_tm_hist(isfinite(rms_tm_hist));

if isempty(rms_joint_hist)
    rms_initial_joint = NaN;
else
    rms_initial_joint = rms_joint_hist(1);
end
if isempty(rms_te_hist)
    rms_initial_te = NaN;
else
    rms_initial_te = rms_te_hist(1);
end
if isempty(rms_tm_hist)
    rms_initial_tm = NaN;
else
    rms_initial_tm = rms_tm_hist(1);
end
if isempty(lambda_hist)
    lambda_initial = NaN;
    lambda_final   = NaN;
else
    lambda_initial = lambda_hist(1);
    lambda_final   = lambda_hist(end);
end

rms_final_joint = rms_final_joint_model;
rms_final_te    = rms_final_te_model;
rms_final_tm    = rms_final_tm_model;
model_rho = 10.^m;

fprintf('\n==================== RELEASE RUN SUMMARY ====================\n');
fprintf('Case                         : %s\n', case_name);
fprintf('Model type                   : %s\n', model_type);
fprintf('Mode                         : %s\n', mode_name);
fprintf('Cached TE Jacobian           : %s\n', file_te_J);
fprintf('Cached TM Jacobian           : %s\n', file_tm_J);
fprintf('Data size                    : TE=%d, TM=%d, joint=%d\n', TE.Ndata, TM.Ndata, Ndata_joint);
fprintf('Model size                   : Nmodel=%d (NDX=%d, NDZ=%d)\n', Nmodel, NDX, TM_NDZ);
fprintf('Joint RMS                    : %.8f -> %.8f\n', rms_initial_joint, rms_final_joint);
fprintf('TE RMS                       : %.8f -> %.8f\n', rms_initial_te, rms_final_te);
fprintf('TM RMS                       : %.8f -> %.8f\n', rms_initial_tm, rms_final_tm);
fprintf('Best joint RMS               : %.8f\n', best_rms_joint);
fprintf('Lambda                       : %.6e -> %.6e\n', lambda_initial, lambda_final);
fprintf('Accepted/rejected steps       : %d / %d\n', n_step_accepted, n_step_rejected);
fprintf('Broyden updates TE/TM         : %d / %d\n', n_broyden_update_te, n_broyden_update_tm);
fprintf('Broyden rejections TE/TM      : %d / %d\n', n_broyden_rejected_te, n_broyden_rejected_tm);
fprintf('Final rho range              : %.6g to %.6g ohm.m\n', min(model_rho), max(model_rho));
fprintf('Final rho mean +/- std       : %.6g +/- %.6g ohm.m\n', mean(model_rho), std(model_rho));
fprintf('Jacobian load/build time      : %.3f s\n', elapsed_jacobian_total_s);
fprintf('Inversion time                : %.3f s\n', elapsed_inversion_s);
fprintf('Total script time             : %.3f s (%.3f min)\n', elapsed_total_s, elapsed_total_s/60);
fprintf('=============================================================\n');

if save_results
    if exist(output_dir, 'dir') ~= 7
        mkdir(output_dir);
    end
    result = struct();
    result.case_name = case_name;
    result.model_type = model_type;
    result.m_final = m;
    result.rho_final = model_rho;
    result.RMSjoint_hist = rms_joint_hist;
    result.RMSTE_hist = rms_te_hist;
    result.RMSTM_hist = rms_tm_hist;
    result.lambda_hist = lambda_hist;
    result.summary = struct( ...
        'rms_initial_joint', rms_initial_joint, ...
        'rms_final_joint', rms_final_joint, ...
        'rms_final_te', rms_final_te, ...
        'rms_final_tm', rms_final_tm, ...
        'best_rms_joint', best_rms_joint, ...
        'elapsed_total_s', elapsed_total_s);
    save(fullfile(output_dir, sprintf('result_%s.mat', model_type)), 'result', '-v7.3');
end

%% ===================== FINAL FORWARD + FIGURES =====================
if show_figures
    try
        [te_rhoa_cal, te_phi_cal_raw, out_te] = call_forward_te_shared( ...
            m, xtopo, ztopo, datafreq, TE_NDZ, ...
            'time_sign', time_sign, 'use_topo', use_topo, ...
            'x_st', x_st, 'ZMIN', TE_ZMIN);

        [tm_rhoa_cal, tm_phi_cal_raw, out_tm] = call_forward_tm_shared( ...
            m, xtopo, ztopo, datafreq, TM_NDZ, ...
            'time_sign', time_sign, 'use_topo', use_topo, ...
            'x_st', x_st, 'ZMIN', TM_ZMIN, ...
            'include_air', include_air_tm, 'param_mode', param_mode);

        te_rhoa_cal(te_rhoa_cal<=0) = realmin;
        tm_rhoa_cal(tm_rhoa_cal<=0) = realmin;
        te_rhocal = log10(te_rhoa_cal);
        tm_rhocal = log10(tm_rhoa_cal);

        switch lower(phase_plot_mode)
            case 'inversion'
                te_phi_cal_plot = apply_phase_mode_deg(te_phi_cal_raw, phase_mode_te);
                te_phi_obs_plot = TE.phi_obs_deg;
                tm_phi_cal_plot = apply_phase_mode_deg(tm_phi_cal_raw, phase_mode_tm);
                tm_phi_obs_plot = TM.phi_obs_deg;
            case 'raw_display'
                te_phi_cal_plot = te_phi_cal_raw;
                te_phi_obs_plot = TE.phi_obs_deg_raw;
                tm_phi_cal_plot = tm_phi_cal_raw;
                tm_phi_obs_plot = TM.phi_obs_deg_raw;
            otherwise
                error('phase_plot_mode must be either inversion or raw_display.');
        end

        plot_model_section_smooth(m, out_te, out_tm, xtopo, ztopo, x_st, nm, save_prefix, model_label, output_dir, save_figures);
        plot_histories_joint(outer_rms_joint_hist, outer_rms_te_hist, outer_rms_tm_hist, ...
            lambda_hist, rms_joint_hist, save_prefix, output_dir, save_figures);
        plot_datafit_joint(F, te_rhocal, log10(TE.rhoa_obs), te_phi_cal_plot, te_phi_obs_plot, ...
            tm_rhocal, log10(TM.rhoa_obs), tm_phi_cal_plot, tm_phi_obs_plot, nm, save_prefix, output_dir, save_figures);

    catch ME
        warning('Final forward calculation or plotting failed: %s', ME.message);
    end
end

%% ===================== LOCAL FUNCTIONS =====================

function [J, J_base, built_flag, metaJ, elapsed_jacobian_s] = load_or_build_direct_jacobian( ...
    file_J, m, Ndata, Nmodel, load_only, force_build, label, ...
    xtopo, ztopo, datafreq, x_st, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
    time_sign, use_topo, include_air_tm, use_phase_in_radian, param_mode, ...
    model_type, rhs_chunk_size, validate_direct_jacobian, n_validate_cols)

has_valid_J = false;
J = [];
metaJ = struct();
elapsed_jacobian_s = 0;
built_flag = 0;

if exist(file_J,'file')==2 && ~force_build
    S = load(file_J);
    if isfield(S,'J')
        Jtmp = double(S.J);
        ok_size = (size(Jtmp,1)==Ndata) && (size(Jtmp,2)==Nmodel);
        if ok_size
            J = Jtmp;
            has_valid_J = true;
            if isfield(S,'meta')
                metaJ = S.meta;
                if isfield(metaJ,'elapsed_jacobian_s')
                    elapsed_jacobian_s = metaJ.elapsed_jacobian_s;
                end
            end
            fprintf('Using cached direct-sensitivity Jacobian %s from %s\n', label, file_J);
        else
            fprintf('Direct-J %s was rejected because of a size mismatch. size(J)=[%d %d], expected=[%d %d].\n', ...
                label, size(Jtmp,1), size(Jtmp,2), Ndata, Nmodel);
        end
    end
end

if ~has_valid_J
    if load_only
        error('Direct-J %s is not valid and load_only_jacobian=true. Set load_only_jacobian=false or build the Jacobian first.', label);
    end

    if exist(file_J,'file')==2 && force_build
        delete(file_J);
    end

    fprintf('\nBuilding direct-sensitivity Jacobian %s ...\n', label);
    tjac = tic;

    if strcmpi(label,'TE')
        if exist('build_J_TE_direct_sensitivity_shared_cell','file') ~= 2
            error('build_J_TE_direct_sensitivity_shared_cell.m was not found on the MATLAB path.');
        end
        [J, metaJ, out0] = build_J_TE_direct_sensitivity_shared_cell( ...
            m, xtopo, ztopo, datafreq, TE_NDZ, ...
            'time_sign', time_sign, ...
            'use_topo', use_topo, ...
            'x_st', x_st, ...
            'ZMIN', TE_ZMIN, ...
            'use_phase_in_radian', use_phase_in_radian, ...
            'rhs_chunk_size', rhs_chunk_size, ...
            'model_type', model_type, ...
            'validate_against_fd', validate_direct_jacobian, ...
            'n_validate_cols', n_validate_cols, ...
            'verbose', true);
    else
        if exist('build_J_TM_direct_sensitivity_shared_cell','file') == 2
            fbuild_tm = 'build_J_TM_direct_sensitivity_shared_cell';
        elseif exist('build_J_TM_direct_sensitivity_cell','file') == 2
            fbuild_tm = 'build_J_TM_direct_sensitivity_cell';
        else
            error('build_J_TM_direct_sensitivity_shared_cell.m / build_J_TM_direct_sensitivity_cell.m was not found on the MATLAB path.');
        end

        [J, out0, timingJ, metaJ] = feval(fbuild_tm, ...
            m, xtopo, ztopo, datafreq, TM_NDZ, ...
            'time_sign', time_sign, ...
            'use_topo', use_topo, ...
            'x_st', x_st, ...
            'ZMIN', TM_ZMIN, ...
            'include_air', include_air_tm, ...
            'param_mode', param_mode, ...
            'use_phase_in_radian', use_phase_in_radian, ...
            'rhs_chunk_size', rhs_chunk_size);
        metaJ.model_type = model_type;
        metaJ.timing = timingJ;
    end

    elapsed_jacobian_s = toc(tjac);

    if size(J,1) ~= Ndata || size(J,2) ~= Nmodel
        error('Built Direct-J %s has a size mismatch. size(J)=[%d %d], expected=[%d %d].', ...
            label, size(J,1), size(J,2), Ndata, Nmodel);
    end

    metaJ.elapsed_jacobian_s = elapsed_jacobian_s;
    metaJ.file_J = file_J;
    metaJ.created_at = datestr(now);
    meta = metaJ; %#ok<NASGU>
    save(file_J, 'J', 'meta', 'out0', '-v7.3');

    built_flag = 1;
    fprintf('Direct-J %s build completed: %s | %.3f s (%.3f min)\n', ...
        label, file_J, elapsed_jacobian_s, elapsed_jacobian_s/60);
end

J_base = J;
end

function [r_joint, r_te_w, r_tm_w] = residual_joint_weighted( ...
    m, xtopo, ztopo, datafreq, x_st, ...
    TE, TM, TE_NDZ, TM_NDZ, TE_ZMIN, TM_ZMIN, ...
    time_sign, use_topo, include_air_tm, use_phase_in_radian, ...
    phase_mode_te, phase_mode_tm, param_mode, ...
    mode_weight_te, mode_weight_tm)

r_te = residual_te_single( ...
    m, xtopo, ztopo, datafreq, TE_NDZ, TE_ZMIN, ...
    time_sign, use_topo, x_st, use_phase_in_radian, phase_mode_te, TE);

r_tm = residual_tm_single( ...
    m, xtopo, ztopo, datafreq, TM_NDZ, TM_ZMIN, ...
    time_sign, use_topo, include_air_tm, x_st, use_phase_in_radian, phase_mode_tm, param_mode, TM);

r_te_w = sqrt(mode_weight_te) * r_te;
r_tm_w = sqrt(mode_weight_tm) * r_tm;
r_joint = [r_te_w; r_tm_w];
end

function r = residual_te_single( ...
    m, xtopo, ztopo, datafreq, NDZ, ZMIN, ...
    time_sign, use_topo, x_st, use_phase_in_radian, phase_mode, TE)

[rhoa_cal, phi_cal_deg_raw] = call_forward_te_shared( ...
    m, xtopo, ztopo, datafreq, NDZ, ...
    'time_sign', time_sign, 'use_topo', use_topo, 'x_st', x_st, 'ZMIN', ZMIN);

rhoa_cal(rhoa_cal<=0) = realmin;
dcal_rho = reshape(log10(rhoa_cal)', [], 1);

phi_cal_deg = apply_phase_mode_deg(phi_cal_deg_raw, phase_mode);
if use_phase_in_radian
    dcal_phi = reshape((phi_cal_deg*pi/180)', [], 1);
else
    dcal_phi = reshape(phi_cal_deg', [], 1);
end

dphi = wrap_phase_diff(dcal_phi - TE.dobs_phi, use_phase_in_radian);
r_rho = (dcal_rho - TE.dobs_rho) ./ TE.sig_rho;
r_phi = dphi ./ TE.sig_phi;
r = [r_rho; r_phi];
end

function r = residual_tm_single( ...
    m, xtopo, ztopo, datafreq, NDZ, ZMIN, ...
    time_sign, use_topo, include_air, x_st, use_phase_in_radian, phase_mode, param_mode, TM)

[rhoa_cal, phi_cal_deg_raw] = call_forward_tm_shared( ...
    m, xtopo, ztopo, datafreq, NDZ, ...
    'time_sign', time_sign, 'use_topo', use_topo, 'include_air', include_air, ...
    'x_st', x_st, 'ZMIN', ZMIN, 'param_mode', param_mode);

rhoa_cal(rhoa_cal<=0) = realmin;
dcal_rho = reshape(log10(rhoa_cal)', [], 1);

phi_cal_deg = apply_phase_mode_deg(phi_cal_deg_raw, phase_mode);
if use_phase_in_radian
    dcal_phi = reshape((phi_cal_deg*pi/180)', [], 1);
else
    dcal_phi = reshape(phi_cal_deg', [], 1);
end

dphi = wrap_phase_diff(dcal_phi - TM.dobs_phi, use_phase_in_radian);
r_rho = (dcal_rho - TM.dobs_rho) ./ TM.sig_rho;
r_phi = dphi ./ TM.sig_phi;
r = [r_rho; r_phi];
end

function [rhoa, phi_deg, out] = call_forward_te_shared(m, xtopo, ztopo, datafreq, NDZ, varargin)
if exist('fungsi_TE_mod_cell_shared','file') == 2
    fwd_name = 'fungsi_TE_mod_cell_shared';
elseif exist('fungsi_TE_mod_cell','file') == 2
    fwd_name = 'fungsi_TE_mod_cell';
elseif exist('fungsi_TE_mod','file') == 2
    fwd_name = 'fungsi_TE_mod';
else
    error('TE forward function was not found.');
end

try
    [rhoa, phi_deg, out] = feval(fwd_name, m, xtopo, ztopo, datafreq, NDZ, varargin{:});
catch ME
    if contains(ME.message, 'Too many output') || contains(ME.message, 'Too many output arguments')
        [rhoa, phi_deg] = feval(fwd_name, m, xtopo, ztopo, datafreq, NDZ, varargin{:});
        out = struct();
    else
        rethrow(ME);
    end
end
end

function [rhoa, phi_deg, out] = call_forward_tm_shared(m, xtopo, ztopo, datafreq, NDZ, varargin)
if exist('fungsi_TM_mod_shared','file') == 2
    fwd_name = 'fungsi_TM_mod_shared';
elseif exist('fungsi_TM_mod','file') == 2
    fwd_name = 'fungsi_TM_mod';
else
    error('TM forward function was not found.');
end

try
    [rhoa, phi_deg, out] = feval(fwd_name, m, xtopo, ztopo, datafreq, NDZ, varargin{:});
catch ME
    if contains(ME.message, 'Too many output') || contains(ME.message, 'Too many output arguments')
        [rhoa, phi_deg] = feval(fwd_name, m, xtopo, ztopo, datafreq, NDZ, varargin{:});
        out = struct();
    else
        rethrow(ME);
    end
end
end

function d = wrap_phase_diff(d, use_rad)
if use_rad
    d = angle(exp(1i*d));
else
    d = mod(d + 180, 360) - 180;
end
end

function phi_out_deg = apply_phase_mode_deg(phi_in_deg, mode)
wrap_deg = @(x) mod(x + 180, 360) - 180;
switch lower(strtrim(mode))
    case 'raw'
        phi_out_deg = phi_in_deg;
    case 'shift180'
        phi_out_deg = wrap_deg(phi_in_deg + 180);
    case 'neg'
        phi_out_deg = wrap_deg(-phi_in_deg);
    case 'neg_shift180'
        phi_out_deg = wrap_deg(-phi_in_deg + 180);
    otherwise
        error('Unknown phase_mode: %s', mode);
end
end

function f = pick_first_existing(cands)
for k = 1:numel(cands)
    if exist(cands{k}, 'file') == 2
        f = cands{k};
        return
    end
end
error('None of the candidate files were found: %s', strjoin(cands, ', '));
end

function plot_model_section_smooth(m, out_te, out_tm, xtopo, ztopo, x_st, nm, save_prefix, model_label, output_dir, save_figures)
fig = figure('Color',[0.94 0.94 0.94], 'Units','pixels','Position',[120 60 760 650]);
ax = axes('Parent',fig);
hold(ax,'on');

[X, Z, NE, model_id, m_tri_log10] = select_plot_mesh(m, out_te, out_tm);

xc_tri = zeros(numel(model_id),1);
zc_tri = zeros(numel(model_id),1);
for k = 1:numel(model_id)
    e = model_id(k);
    xc_tri(k) = mean(X(NE(:,e)));
    zc_tri(k) = mean(Z(NE(:,e)));
end

xq = linspace(-500, 8000, 420);
zq = linspace(-5000, 5000, 420);
[Xq,Zq] = meshgrid(xq,zq);

Fint = scatteredInterpolant(xc_tri, zc_tri, m_tri_log10(:), 'natural', 'nearest');
Mq = Fint(Xq,Zq);

zsurf_q = interp1(xtopo(:), ztopo(:), xq, 'linear', 'extrap');
for ix = 1:numel(xq)
    Mq(Zq(:,ix) < zsurf_q(ix), ix) = NaN;
end

contourf(ax, Xq, Zq, Mq, 45, 'LineStyle','none');

if exist('cmap2.mat','file') == 2
    load cmap2
    colormap(ax, cmap2);
else
    colormap(ax, flipud(jet));
end
clim(ax, [0 3]);

% Mask the air region so it is not colored.
air_color = get(ax, 'Color');
z_air_top = min(zq);
x_air = [xq, fliplr(xq)];
z_air = [zsurf_q, z_air_top*ones(size(xq))];
patch('Parent', ax, 'XData', x_air, 'YData', z_air, ...
      'FaceColor', air_color, 'EdgeColor', 'none');

% Topography line. Comment this line if no surface line is desired.
plot(ax, xtopo, ztopo, 'k-', 'LineWidth', 2.0);

x_st_plot = x_st(:);
n_st_plot = numel(x_st_plot);
if numel(nm) ~= n_st_plot
    error('Number of station names (%d) does not match number of plotted stations (%d).', numel(nm), n_st_plot);
end

zsurf_st = interp1(xtopo(:), ztopo(:), x_st_plot, 'linear', 'extrap');
marker_offset = 120;
label_gap     = 650;
x_shift       = 150;

z_marker = zsurf_st - marker_offset;
z_label  = z_marker - label_gap;
x_label  = x_st_plot + x_shift;

plot(ax, x_st_plot, z_marker, 'vk', 'MarkerFaceColor','k', 'MarkerSize',7, 'LineWidth',0.8);
for i = 1:n_st_plot
    text(ax, x_label(i), z_label(i), nm{i}, ...
        'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',11, 'FontWeight','normal', 'FontName','cambria', ...
        'Interpreter','none', 'Clipping','off', 'Color','k');
end

set(ax, 'YDir','reverse', 'FontSize',14, 'FontName','cambria', ...
    'LineWidth',1.2, 'Box','on', 'TickDir','in', 'Layer','top');

xlim(ax, [-500 8000]);
ylim(ax, [-5000 5000]);
xticks(ax, 0:2000:8000);
yticks(ax, [-5000 -3000 -1000 0 1000 3000 5000]);

xlabel(ax, 'X (m)', 'FontSize',16, 'FontWeight','bold', 'FontName','cambria');
ylabel(ax, 'Depth (m)', 'FontSize',16, 'FontWeight','bold', 'FontName','cambria');
title(ax, model_label, 'FontSize',18, 'FontWeight','bold', 'FontName','cambria', 'Interpreter','tex');

ax.Position = [0.12 0.10 0.63 0.82];

cb = colorbar(ax);
cb.Position = [0.80 0.10 0.035 0.82];
ticks_ohm = [1 3.1 10 31.6 100 316 1000];
cb.Ticks = log10(ticks_ohm);
cb.TickLabels = arrayfun(@(v)sprintf('%g',v), ticks_ohm, 'UniformOutput', false);
cb.FontSize = 12;
cb.FontName = 'cambria';
title(cb, '\rho (\Omega m)', 'FontSize',14, 'FontWeight','bold', 'FontName','cambria');

if save_figures
    saveas(fig, fullfile(output_dir, sprintf('Model_%s.png',save_prefix)));
end
if save_figures
    saveas(fig, fullfile(output_dir, sprintf('Model_%s.fig',save_prefix)));
end
end

function [X, Z, NE, model_id, m_tri_log10] = select_plot_mesh(m, out_te, out_tm)
% Prefer the TE mesh because it usually follows topography and the shared earth model.
if isstruct(out_te) && isfield(out_te,'X') && isfield(out_te,'Z') && isfield(out_te,'NE') && isfield(out_te,'model_id')
    X = out_te.X;
    Z = out_te.Z;
    NE = out_te.NE;
    model_id = out_te.model_id(:);
    if numel(model_id) == 4*numel(m)
        m_tri_log10 = repelem(m(:),4);
        return
    elseif isfield(out_te,'sigma')
        rho_tri = 1 ./ out_te.sigma(model_id);
        m_tri_log10 = log10(rho_tri(:));
        return
    end
end

if isstruct(out_tm) && isfield(out_tm,'X') && isfield(out_tm,'Z') && isfield(out_tm,'NE') && isfield(out_tm,'model_id')
    X = out_tm.X;
    Z = out_tm.Z;
    NE = out_tm.NE;
    model_id = out_tm.model_id(:);
    if numel(model_id) == 4*numel(m)
        m_tri_log10 = repelem(m(:),4);
        return
    elseif isfield(out_tm,'sigma')
        rho_tri = 1 ./ out_tm.sigma(model_id);
        m_tri_log10 = log10(rho_tri(:));
        return
    end
end

error('Cannot build model plot: out_te/out_tm does not contain valid X, Z, NE, and model_id fields.');
end

function plot_histories_joint(outer_joint, outer_te, outer_tm, lambda_hist, rms_joint_hist, save_prefix, output_dir, save_figures)
figure('Color','w');
plot(outer_joint,'-o','LineWidth',2); hold on;
plot(outer_te,'--','LineWidth',2);
plot(outer_tm,'-.','LineWidth',2);
grid on; box on;
xlabel('Outer iteration');
ylabel('Best RMS');
legend('Joint','TE','TM','location','best');
title('Convergence RMS ', 'FontWeight','bold');
set(gca,'fontsize',14,'fontname','cambria');
if save_figures
    saveas(gcf, fullfile(output_dir, sprintf('OuterRMS_%s.png',save_prefix)));
end
if save_figures
    saveas(gcf, fullfile(output_dir, sprintf('OuterRMS_%s.fig',save_prefix)));
end

figure('Color','w');
semilogy(lambda_hist,'LineWidth',2);
grid on; box on;
xlabel('Iteration');
ylabel('\lambda');
title('Lambda history', 'FontWeight','bold');
set(gca,'fontsize',14,'fontname','cambria');
if save_figures
    saveas(gcf, fullfile(output_dir, sprintf('Lambda_%s.png',save_prefix)));
end
if save_figures
    saveas(gcf, fullfile(output_dir, sprintf('Lambda_%s.fig',save_prefix)));
end

figure('Color','w');
plot(rms_joint_hist,'LineWidth',3);
grid on; box on;
xlabel('Iteration');
ylabel('RMS joint');
title('Weighted RMS Joint', 'FontWeight','bold');
set(gca,'fontsize',15,'fontname','cambria');
if save_figures
    saveas(gcf, fullfile(output_dir, sprintf('Misfit_%s.png',save_prefix)));
end
if save_figures
    saveas(gcf, fullfile(output_dir, sprintf('Misfit_%s.fig',save_prefix)));
end
end

function plot_datafit_joint(F, te_cal_logrho, te_obs_logrho, te_phi_cal, te_phi_obs, ...
                               tm_cal_logrho, tm_obs_logrho, tm_phi_cal, tm_phi_obs, nm, save_prefix, output_dir, save_figures)
n_st = size(te_cal_logrho,2);
if numel(nm) ~= n_st
    nm = arrayfun(@(k) sprintf('MT-%02d',k), 1:n_st, 'UniformOutput', false).';
end

col_te = [0 0.4470 0.7410];
col_tm = [1 0 0];

for i = 1:n_st
    fig = figure('Color','w','Units','pixels','Position',[180 60 760 900]);

    ax1 = subplot(2,1,1);
    semilogx(F, te_cal_logrho(:,i), '-', 'Color', col_te, 'LineWidth', 2.5); hold on;
    semilogx(F, te_obs_logrho(:,i), 'o', 'Color', col_te, 'MarkerFaceColor', col_te, 'LineWidth', 1.2, 'MarkerSize', 4);
    semilogx(F, tm_cal_logrho(:,i), '-', 'Color', col_tm, 'LineWidth', 2.5);
    semilogx(F, tm_obs_logrho(:,i), 'o', 'Color', col_tm, 'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 4);

    xlim([1e-2 1e3]);
    set(gca,'xdir','reverse');
    xticks([1e-2 1e-1 1 10 100 1000]);
    xticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'});
    ylim([-1 4]);
    yticks([-1 0 1 2 3 4]);

    title(sprintf('Resistivity - %s', nm{i}), 'fontweight','bold');
    ylabel('Log App.Res (\Omega m)');
    xlabel('Freq. (Hz)');
    legend('TE cal','TE obs','TM cal','TM obs','Location','northwest','Orientation','horizontal');
    box on; grid on;
    set(gca,'fontsize',14,'fontname','cambria');

    ax2 = subplot(2,1,2);
    semilogx(F, te_phi_cal(:,i), '-', 'Color', col_te, 'LineWidth', 2.5); hold on;
    semilogx(F, te_phi_obs(:,i), 'o', 'Color', col_te, 'MarkerFaceColor', col_te, 'LineWidth', 1.2, 'MarkerSize', 4);
    semilogx(F, tm_phi_cal(:,i), '-', 'Color', col_tm, 'LineWidth', 2.5);
    semilogx(F, tm_phi_obs(:,i), 'o', 'Color', col_tm, 'MarkerFaceColor', 'none', 'LineWidth', 1.2, 'MarkerSize', 4);

    xlim([1e-2 1e3]);
    set(gca,'xdir','reverse');
    xticks([1e-2 1e-1 1 10 100 1000]);
    xticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'});
    ylim([-250 250]);
    yticks([-200 -100 0 100 200]);

    title(sprintf('Phase - %s', nm{i}), 'fontweight','bold');
    ylabel('Phase (\circ)');
    xlabel('Freq. (Hz)');
    legend('TE cal','TE obs','TM cal','TM obs','Location','northwest','Orientation','horizontal');
    box on; grid on;
    set(gca,'fontsize',14,'fontname','cambria');

    set(ax1, 'Position', [0.14 0.46 0.72 0.50]);
    set(ax2, 'Position', [0.14 0.10 0.72 0.24]);

    if save_figures
        saveas(fig, fullfile(output_dir, sprintf('FitGabung_%s_%02d.png', save_prefix, i)));
    end
    if save_figures
        saveas(fig, fullfile(output_dir, sprintf('FitGabung_%s_%02d.fig', save_prefix, i)));
    end
end
end
