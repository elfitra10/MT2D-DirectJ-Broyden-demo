% Clean release copy. This function builds the TE direct-sensitivity Jacobian.

function [J, meta, out0, cache] = build_J_TE_direct_sensitivity_shared_cell(m_log10rho_cell, xtopo, ztopo, datafreq, NDZ, varargin)
%BUILD_J_TE_DIRECT_SENSITIVITY_SHARED_CELL
% TE direct-sensitivity Jacobian 
%
% Main syntax:
%   [J, meta, out0] = build_J_TE_direct_sensitivity_cell(m0, xtopo, ztopo, datafreq, NDZ, ...)
%
% The Jacobian follows this data ordering:
%   d_rho = reshape(log10(rhoa)', [], 1);
%   d_phi = reshape(phi', [], 1);
%   d     = [d_rho; d_phi];
%
% Sensitivity equation:
%   K(m) u = rhs(m)
%   K du/dm_j = - dM/dm_j * E
%
% This file requires funcDZearth_shared.m and konektivitas.m on the MATLAB path.

%% ===================== OPTIONS =====================
opt = struct();
opt.ZMIN                = -50000;
opt.use_topo            = any(abs(ztopo(:)) > 1e-9);
opt.time_sign           = +1;
opt.use_phase_in_radian = true;
opt.x_st                = [];
opt.rho_air             = 1e8;
opt.rho_bg              = 100;
opt.rhs_chunk_size      = 100;
opt.verbose             = true;
opt.model_type          = 'unknown';
opt.validate_against_fd = false;
opt.n_validate_cols     = 6;
opt.validation_seed     = 7;
opt.delta_rel_fd        = 1e-2;
opt.delta_abs_fd        = 5e-3;

if mod(numel(varargin),2) ~= 0
    error('Optional arguments must be provided as paired ''key'', value entries.');
end
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k+1};
    switch key
        case "zmin"
            opt.ZMIN = val;
        case "use_topo"
            opt.use_topo = logical(val);
        case "time_sign"
            opt.time_sign = val;
        case "use_phase_in_radian"
            opt.use_phase_in_radian = logical(val);
        case "x_st"
            opt.x_st = val;
        case "rho_air"
            opt.rho_air = val;
        case {"rho_bg","rho_hs"}
            opt.rho_bg = val;
        case "rhs_chunk_size"
            opt.rhs_chunk_size = val;
        case "verbose"
            opt.verbose = logical(val);
        case "model_type"
            opt.model_type = char(val);
        case {"validate","validate_against_fd"}
            opt.validate_against_fd = logical(val);
        case "n_validate_cols"
            opt.n_validate_cols = val;
        case "validation_seed"
            opt.validation_seed = val;
        case "delta_rel_fd"
            opt.delta_rel_fd = val;
        case "delta_abs_fd"
            opt.delta_abs_fd = val;
        otherwise
            error('Unknown optional key: %s', key);
    end
end

assert(exist('funcDZearth_shared','file')==2, 'funcDZearth_shared.m was not found on the MATLAB path.');
assert(exist('konektivitas','file')==2, 'konektivitas.m was not found on the MATLAB path.');

m_log10rho_cell = m_log10rho_cell(:);
xtopo = xtopo(:);
ztopo = ztopo(:);
F = datafreq(:,1);
F = F(:);
NF = numel(F);

if isempty(opt.x_st)
    error('x_st must be provided to build the station Jacobian.');
end
x_st = opt.x_st(:);
n_st = numel(x_st);

%% ===================== BUILD CACHE / MESH =====================
t_all = tic;
cache = build_te_cache_direct_local(xtopo, ztopo, NDZ, opt.ZMIN, opt.use_topo);
NDX = cache.NDX;
NDZb = cache.NDZb;
Nmodel = NDX * NDZb;
Ndata  = 2 * NF * n_st;

if numel(m_log10rho_cell) ~= Nmodel
    error('Model size mismatch: numel(m)=%d, expected Nmodel=%d.', numel(m_log10rho_cell), Nmodel);
end

if opt.verbose
    fprintf('\n=== BUILD DIRECT-SENSITIVITY JACOBIAN TE ===\n');
    fprintf('model_type = %s\n', opt.model_type);
    fprintf('NF         = %d\n', NF);
    fprintf('nStation   = %d\n', n_st);
    fprintf('NDX        = %d\n', NDX);
    fprintf('NDZ        = %d\n', NDZ);
    fprintf('NDZb       = %d\n', NDZb);
    fprintf('Nmodel     = %d\n', Nmodel);
    fprintf('Ndata      = %d\n', Ndata);
    fprintf('rhs chunk  = %d\n', opt.rhs_chunk_size);
end

%% ===================== BUILD J =====================
t_jac = tic;
[J, out0, timing] = build_J_TE_direct_core_local( ...
    m_log10rho_cell, cache, datafreq, x_st, opt.time_sign, opt.use_phase_in_radian, ...
    opt.rho_air, opt.rho_bg, opt.rhs_chunk_size, opt.verbose);
elapsed_jacobian_s = toc(t_jac);

if opt.verbose
    fprintf('\nDirect-sensitivity Jacobian completed.\n');
    fprintf('elapsed_jacobian_s = %.6f s (%.3f min)\n', elapsed_jacobian_s, elapsed_jacobian_s/60);
end

%% ===================== VALIDATION OPTIONAL =====================
validation = struct();
validation.enabled = opt.validate_against_fd;

if opt.validate_against_fd
    if opt.verbose
        fprintf('\n=== VALIDATION VS FINITE DIFFERENCE ===\n');
    end

    rng(opt.validation_seed);
    cols = randperm(Nmodel, min(opt.n_validate_cols, Nmodel));
    relerr  = nan(numel(cols),1);
    corrval = nan(numel(cols),1);
    norm_fd = nan(numel(cols),1);
    norm_ds = nan(numel(cols),1);

    [rhoa0_fd, phi0_deg_fd] = call_forward_te_direct_local(m_log10rho_cell, xtopo, ztopo, datafreq, NDZ, ...
        'time_sign', opt.time_sign, 'use_topo', opt.use_topo, 'x_st', x_st, 'ZMIN', opt.ZMIN, ...
        'rho_air', opt.rho_air, 'rho_bg', opt.rho_bg);

    rhoa0_fd(rhoa0_fd <= 0) = realmin;
    d0_rho = reshape(log10(rhoa0_fd)', [], 1);
    if opt.use_phase_in_radian
        d0_phi = reshape((phi0_deg_fd*pi/180)', [], 1);
    else
        d0_phi = reshape(phi0_deg_fd', [], 1);
    end
    d0 = [d0_rho; d0_phi];

    for kk = 1:numel(cols)
        ip = cols(kk);
        del = max(opt.delta_abs_fd, opt.delta_rel_fd*abs(m_log10rho_cell(ip)));

        m1 = m_log10rho_cell;
        m1(ip) = m1(ip) + del;

        [rhoa1, phi1_deg] = call_forward_te_direct_local(m1, xtopo, ztopo, datafreq, NDZ, ...
            'time_sign', opt.time_sign, 'use_topo', opt.use_topo, 'x_st', x_st, 'ZMIN', opt.ZMIN, ...
            'rho_air', opt.rho_air, 'rho_bg', opt.rho_bg);

        rhoa1(rhoa1 <= 0) = realmin;
        d1_rho = reshape(log10(rhoa1)', [], 1);
        if opt.use_phase_in_radian
            d1_phi = reshape((phi1_deg*pi/180)', [], 1);
        else
            d1_phi = reshape(phi1_deg', [], 1);
        end
        d1 = [d1_rho; d1_phi];

        Jfd = (d1 - d0) / del;
        Jds = J(:,ip);

        relerr(kk) = norm(Jds - Jfd) / max(norm(Jfd), eps);
        norm_fd(kk) = norm(Jfd);
        norm_ds(kk) = norm(Jds);
        C = corrcoef(double(Jfd), double(Jds));
        if numel(C) == 4
            corrval(kk) = C(1,2);
        end

        if opt.verbose
            fprintf('col %5d | relerr = %.4e | corr = %.4f | normFD = %.4e | normDS = %.4e\n', ...
                ip, relerr(kk), corrval(kk), norm_fd(kk), norm_ds(kk));
        end
    end

    validation.cols = cols(:);
    validation.relerr = relerr;
    validation.corrval = corrval;
    validation.norm_fd = norm_fd;
    validation.norm_ds = norm_ds;
    validation.mean_relerr = mean(relerr,'omitnan');
    validation.median_relerr = median(relerr,'omitnan');

    if opt.verbose
        fprintf('Mean relerr   = %.4e\n', validation.mean_relerr);
        fprintf('Median relerr = %.4e\n', validation.median_relerr);
    end
end

%% ===================== META =====================
meta = struct();
meta.mode                  = 'TE';
meta.model_type            = opt.model_type;
meta.jacobian_type         = 'direct_sensitivity';
meta.param_mode            = 'cell';
meta.use_topo              = opt.use_topo;
meta.NDZ                   = NDZ;
meta.ZMIN                  = opt.ZMIN;
meta.NDX                   = NDX;
meta.NDZb                  = NDZb;
meta.Nmodel                = Nmodel;
meta.Ndata                 = Ndata;
meta.NF                    = NF;
meta.n_station             = n_st;
meta.time_sign             = opt.time_sign;
meta.use_phase_in_radian   = opt.use_phase_in_radian;
meta.rho_air               = opt.rho_air;
meta.rho_bg                = opt.rho_bg;
meta.xtopo                 = xtopo(:);
meta.ztopo                 = ztopo(:);
meta.x_st                  = x_st(:);
meta.F                     = F(:);
meta.rhs_chunk_size        = opt.rhs_chunk_size;
meta.elapsed_jacobian_s    = elapsed_jacobian_s;
meta.elapsed_total_s       = toc(t_all);
meta.timing                = timing;
meta.validation            = validation;
end

%% ========================================================================
function cache = build_te_cache_direct_local(xtopo, ztopo, NDZ, ZMIN, use_topo)
xtopo = xtopo(:);
ztopo = ztopo(:);

NDX = length(xtopo)-1;
a = 1; b = 3;
NDZa = round(a*NDZ/(a+b));
NDZb = NDZ - NDZa;
NGX = NDX + 1;
NGZ = NDZ + 1;
NTE = NDX * NDZ * 4;
NN  = (NGX*NGZ) + (NDX*NDZ);
DX = xtopo(:);

DZa = logspace(log10(25), log10(abs(ZMIN)), NDZa).';
DZa(end) = abs(ZMIN);

DZb = funcDZearth_shared(NDZb);
DZb = DZb(:);
while ~isempty(DZb) && abs(DZb(1)) < 1e-12
    DZb(1) = [];
end
DZ = [-flipud(DZa); 0; DZb];

k = 0;
X = zeros(NN,1);
Z = zeros(NN,1);
idxpl = zeros(NGX,NGZ);

for j = 1:NGZ
    for i = 1:NGX
        k = k+1;
        X(k) = DX(i);
        Z(k) = DZ(j);
        idxpl(i,j) = k;
    end
    k = k + NDX;
end

k = NGX;
for j = 1:(NGZ-1)
    for i = 1:(NGX-1)
        k = k+1;
        X(k) = (DX(i)+DX(i+1))/2;
        Z(k) = (DZ(j)+DZ(j+1))/2;
    end
    k = k + NGX;
end

Zref = Z;
if use_topo
    zt = interp1(xtopo(:), ztopo(:), X, 'makima', 'extrap');
else
    zt = zeros(size(X));
end

ZTOP = min(Zref);
ZBOT = max(Zref);
w = zeros(size(Zref));
idxE = (Zref >= 0);
w(idxE) = 1 - (Zref(idxE) / ZBOT);
idxA = (Zref < 0);
w(idxA) = (Zref(idxA) - ZTOP) / (0 - ZTOP);
w = max(0, min(1, w));
Z = Zref + w .* zt;

NE = konektivitas(NTE, NGX, NDX, NDZ);
zc = (Zref(NE(1,:)) + Zref(NE(2,:)) + Zref(NE(3,:))) / 3;
model_id = find(zc > 0);

n1 = NE(1,:).';
n2 = NE(2,:).';
n3 = NE(3,:).';
be1 = Z(n2) - Z(n3);
be2 = Z(n3) - Z(n1);
be3 = Z(n1) - Z(n2);
ce1 = X(n3) - X(n2);
ce2 = X(n1) - X(n3);
ce3 = X(n2) - X(n1);
deltae = abs(0.5*(be1.*ce2 - be2.*ce1));
if any(deltae < 1e-20)
    error('An element with near-zero area was found. Check the mesh/topography.');
end
Bmat = [be1 be2 be3];
Cmat = [ce1 ce2 ce3];
nodes = [n1 n2 n3];

I = zeros(9*NTE,1);
J = zeros(9*NTE,1);
Vkg = zeros(9*NTE,1);
mass_base = zeros(9*NTE,1);
idx = 0;
for aij = 1:3
    for bij = 1:3
        idx = idx + 1;
        s = (idx-1)*NTE + (1:NTE);
        I(s) = nodes(:,aij);
        J(s) = nodes(:,bij);
        Vkg(s) = (Bmat(:,aij).*Bmat(:,bij) + Cmat(:,aij).*Cmat(:,bij)) ./ (4*deltae);
        mass_base(s) = (1 + (aij==bij)) .* (deltae/12);
    end
end

Kg = sparse(I, J, Vkg, NN, NN);

bot_nodes = idxpl(:,NGZ);
Irb = zeros(4*(NGX-1),1);
Jrb = zeros(4*(NGX-1),1);
Vrb = zeros(4*(NGX-1),1);
kk = 0;
for ee = 1:(NGX-1)
    nb1 = bot_nodes(ee);
    nb2 = bot_nodes(ee+1);
    Ledge = hypot(X(nb2)-X(nb1), Z(nb2)-Z(nb1));
    base = Ledge/6;
    kk=kk+1; Irb(kk)=nb1; Jrb(kk)=nb1; Vrb(kk)=2*base;
    kk=kk+1; Irb(kk)=nb1; Jrb(kk)=nb2; Vrb(kk)=1*base;
    kk=kk+1; Irb(kk)=nb2; Jrb(kk)=nb1; Vrb(kk)=1*base;
    kk=kk+1; Irb(kk)=nb2; Jrb(kk)=nb2; Vrb(kk)=2*base;
end
Rb = sparse(Irb, Jrb, Vrb, NN, NN);

D = (1:NGX).';
U = (NGX+1:NN).';
Pdir = ones(NGX,1);

Kg_uu = Kg(U,U);
Kg_ud = Kg(U,D);
Rb_uu = Rb(U,U);
Rb_ud = Rb(U,D);

j0 = find(DZ==0, 1, 'first');
if isempty(j0) || j0<=1 || j0>=numel(DZ)
    error('DZ==0 was not found or is invalid.');
end
kS = idxpl(:, j0);
kU = idxpl(:, j0-1);
kD = idxpl(:, j0+1);

cell2tri = cell(NDX*NDZb,1);
for izb = 1:NDZb
    j_box_global = NDZa + izb;
    for ix = 1:NDX
        cell_id = (izb-1)*NDX + ix;
        box_id = (j_box_global-1)*NDX + ix;
        tri_ids = (box_id-1)*4 + (1:4);
        cell2tri{cell_id} = tri_ids(:);
    end
end

g2u = zeros(NN,1);
g2u(U) = 1:numel(U);

cache = struct();
cache.NDX = NDX;
cache.NDZ = NDZ;
cache.NDZa = NDZa;
cache.NDZb = NDZb;
cache.NGX = NGX;
cache.NGZ = NGZ;
cache.NTE = NTE;
cache.NN = NN;
cache.DX = DX;
cache.DZ = DZ;
cache.X = X;
cache.Z = Z;
cache.Zref = Zref;
cache.idxpl = idxpl;
cache.NE = NE;
cache.model_id = model_id;
cache.cell2tri = cell2tri;
cache.I = I;
cache.J = J;
cache.Vkg = Vkg;
cache.mass_base = mass_base;
cache.Kg_uu = Kg_uu;
cache.Kg_ud = Kg_ud;
cache.Rb_uu = Rb_uu;
cache.Rb_ud = Rb_ud;
cache.D = D;
cache.U = U;
cache.Pdir = Pdir;
cache.kS = kS;
cache.kU = kU;
cache.kD = kD;
cache.g2u = g2u;
end

%% ========================================================================
function [J, out0, timing] = build_J_TE_direct_core_local( ...
    m_log10rho_cell, cache, datafreq, x_st, time_sign, use_phase_in_radian, ...
    rho_air, rho_bg, rhs_chunk_size, verbose)

mu0 = 4*pi*1e-7;
sigma_air = 1/rho_air;
sigma_hs  = 1/rho_bg;
F = datafreq(:,1);
F = F(:);
NF = numel(F);
x_st = x_st(:);
n_st = numel(x_st);

NDX = cache.NDX;
NDZb = cache.NDZb;
Nmodel = NDX*NDZb;
Ndata = 2*NF*n_st;
NTE = cache.NTE;
NN = cache.NN;
U = cache.U;
D = cache.D;
Pdir = cache.Pdir;
DX = cache.DX;
Z = cache.Z;
model_id = cache.model_id;
cell2tri = cache.cell2tri;

ix_st = zeros(n_st,1);
for ss = 1:n_st
    [~, ix_st(ss)] = min(abs(DX - x_st(ss)));
end

tri_global_to_local = zeros(NTE,1);
tri_global_to_local(model_id) = 1:numel(model_id);
m_tri_log10 = zeros(numel(model_id),1);
for ic = 1:Nmodel
    tri_ids = cell2tri{ic};
    tri_loc = tri_global_to_local(tri_ids);
    m_tri_log10(tri_loc) = m_log10rho_cell(ic);
end

sigma = sigma_air * ones(NTE,1);
rho_earth_tri = 10.^m_tri_log10(:);
sigma(model_id) = 1 ./ rho_earth_tri;

sigma_rep = repmat(sigma, 9, 1);
M = sparse(cache.I, cache.J, cache.mass_base .* sigma_rep, NN, NN);
M_uu = M(U,U);
M_ud = M(U,D);

Kg_uu = cache.Kg_uu;
Kg_ud = cache.Kg_ud;
Rb_uu = cache.Rb_uu;
Rb_ud = cache.Rb_ud;
kS_all = cache.kS;
kU_all = cache.kU;
kD_all = cache.kD;

kS = kS_all(ix_st);
kU = kU_all(ix_st);
kD = kD_all(ix_st);
rowS = cache.g2u(kS);
rowU = cache.g2u(kU);
rowD = cache.g2u(kD);
if any(rowS==0) || any(rowU==0) || any(rowD==0)
    error('The station surface/upper/lower node is on the Dirichlet/top boundary. Check the mesh.');
end

dz_st = Z(kD) - Z(kU);

J = zeros(Ndata, Nmodel, 'double');
Rhoapp  = zeros(NF, n_st);
FasaImp = zeros(NF, n_st);

freq_time = zeros(NF,1);
factor_solve_time = zeros(NF,1);
rhs_assemble_time = zeros(NF,1);
fill_time = zeros(NF,1);

if verbose
    fprintf('\nBuilding direct-sensitivity J per frequency...\n');
end

for ifr = 1:NF
    tfreq = tic;
    freq = F(ifr);
    omega = 2*pi*freq;
    omega_miu = omega * mu0;
    beta  = time_sign * 1i * omega_miu;
    alpha = 1 / beta;

    kbot = sqrt(beta * sigma_hs);
    if real(kbot) < 0 || (abs(real(kbot))<1e-14 && imag(kbot)<0)
        kbot = -kbot;
    end
    gamma = alpha * kbot;

    Kuu = alpha*Kg_uu + M_uu + gamma*Rb_uu;
    rhs = -(alpha*Kg_ud + M_ud + gamma*Rb_ud) * Pdir;

    tfac = tic;
    Kfac = decomposition(Kuu, 'lu');
    Eu = Kfac \ rhs;
    factor_solve_time(ifr) = toc(tfac);

    E = complex(zeros(NN,1));
    E(U) = Eu;
    E(D) = Pdir;

    dE_dz = (E(kD_all) - E(kU_all)) ./ (Z(kD_all) - Z(kU_all));
    Hx_all = (1i/omega_miu) .* dE_dz;
    Zte_all = E(kS_all) ./ Hx_all;

    Zte_st = Zte_all(ix_st);
    Hx_st = Hx_all(ix_st);
    E_st = E(kS);

    Rhoapp(ifr,:)  = ((abs(Zte_st).^2) ./ omega_miu).';
    FasaImp(ifr,:) = (angle(Zte_st)*180/pi).';

    trhs = tic;
    B = assemble_rhs_direct_all_cells_local(E, sigma, cache);
    rhs_assemble_time(ifr) = toc(trhs);

    % Ordering follows reshape(A',[],1): [all stations at freq1; all stations at freq2; ...]
    rho_rows = (ifr-1)*n_st + (1:n_st).';
    phi_rows = NF*n_st + rho_rows;

    tfill = tic;
    for c1 = 1:rhs_chunk_size:Nmodel
        c2 = min(c1 + rhs_chunk_size - 1, Nmodel);
        cols = c1:c2;

        DU = Kfac \ B(:,cols);

        dE_S = DU(rowS,:);
        dE_U = DU(rowU,:);
        dE_D = DU(rowD,:);

        dHx = (1i/omega_miu) .* ((dE_D - dE_U) ./ dz_st);
        dZ  = (dE_S .* Hx_st - E_st .* dHx) ./ (Hx_st.^2);

        dlogrho = (2/log(10)) * real(dZ ./ Zte_st);
        dphi_rad = imag(dZ ./ Zte_st);
        if use_phase_in_radian
            dphi = dphi_rad;
        else
            dphi = dphi_rad * 180/pi;
        end

        J(rho_rows, cols) = dlogrho;
        J(phi_rows, cols) = dphi;
    end
    fill_time(ifr) = toc(tfill);

    freq_time(ifr) = toc(tfreq);
    if verbose
        fprintf('  freq %3d/%3d | f=%10.4g Hz | %.2f s\n', ifr, NF, freq, freq_time(ifr));
    end
end

out0 = struct();
out0.Rhoapp = Rhoapp;
out0.PhiDeg = FasaImp;
out0.DX = DX;
out0.x_st = x_st;
out0.model_id = model_id;
out0.NGX = cache.NGX;
out0.NTE = cache.NTE;
out0.NN = cache.NN;
out0.X = cache.X;
out0.Z = cache.Z;
out0.DZ = cache.DZ;
out0.NE = cache.NE;
out0.idxpl = cache.idxpl;
out0.Zref = cache.Zref;
out0.sigma = sigma;
out0.cell2tri = cell2tri;

timing = struct();
timing.freq_time = freq_time;
timing.factor_solve_time = factor_solve_time;
timing.rhs_assemble_time = rhs_assemble_time;
timing.fill_time = fill_time;
timing.total_freq_time = sum(freq_time);
timing.total_factor_solve_time = sum(factor_solve_time);
timing.total_rhs_assemble_time = sum(rhs_assemble_time);
timing.total_fill_time = sum(fill_time);
end

%% ========================================================================
function B = assemble_rhs_direct_all_cells_local(E, sigma, cache)
Nmodel = cache.NDX * cache.NDZb;
NTE = cache.NTE;
U = cache.U;
g2u = cache.g2u;
maxnnz = 36 * Nmodel;
Irow = zeros(maxnnz,1);
Jcol = zeros(maxnnz,1);
Vval = complex(zeros(maxnnz,1));
ptr = 0;

for ic = 1:Nmodel
    tri_ids = cache.cell2tri{ic};
    for tt = 1:numel(tri_ids)
        e = tri_ids(tt);
        dsig_dm = -log(10) * sigma(e);
        idxs = e + (0:8)'*NTE;
        rowsG = cache.I(idxs);
        colsG = cache.J(idxs);
        vals = -cache.mass_base(idxs) .* dsig_dm .* E(colsG);

        for q = 1:numel(idxs)
            ru = g2u(rowsG(q));
            if ru > 0
                ptr = ptr + 1;
                Irow(ptr) = ru;
                Jcol(ptr) = ic;
                Vval(ptr) = vals(q);
            end
        end
    end
end

B = sparse(Irow(1:ptr), Jcol(1:ptr), Vval(1:ptr), numel(U), Nmodel);
end

%% ========================================================================
function [rhoa, phi_deg, out] = call_forward_te_direct_local(m, xtopo, ztopo, datafreq, NDZ, varargin)
if exist('fungsi_TE_mod_cell_shared','file') == 2
    fwd_name = 'fungsi_TE_mod_cell_shared';
elseif exist('fungsi_TE_mod_cell','file') == 2
    fwd_name = 'fungsi_TE_mod_cell';
elseif exist('fungsi_TE_mod_cell_lap','file') == 2
    fwd_name = 'fungsi_TE_mod_cell_lap';
elseif exist('fungsi_TE_mod','file') == 2
    fwd_name = 'fungsi_TE_mod';
else
    error('TE forward solver was not found: fungsi_TE_mod_cell_shared / fungsi_TE_mod_cell / fungsi_TE_mod_cell_lap / fungsi_TE_mod.');
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
