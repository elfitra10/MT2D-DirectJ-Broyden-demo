% Clean release copy. This function builds the TM direct-sensitivity Jacobian.

function [J, out0, timing, meta] = build_J_TM_direct_sensitivity_shared_cell(m_log10rho_cell, xtopo, ztopo, datafreq, NDZ, varargin)
%BUILD_J_TM_DIRECT_SENSITIVITY_SHARED_CELL
% Direct-sensitivity Jacobian for 2-D TM MT (Hy formulation) with cell parameterization.
%
% This function follows the matrix structure used in fungsi_TM_mod.m:
%   K(m) H = b(m)
% and computes sensitivity from
%   K dH/dm = - (dK/dm) H
% for free nodes after Dirichlet elimination.
%
% Output J order is consistent with reshape(K',[],1):
%   [log10(rhoa) station/frequency-major; phase station/frequency-major]
%
% Required local/user functions in MATLAB path:
%   funcDZ.m (optional; fallback available)
%
% Notes:
%   - Designed primarily for param_mode='cell' and include_air=false.
%   - The post-processing derivative includes sigma_below sensitivity in Ex.

p = inputParser;
p.addParameter('time_sign', 1, @isnumeric);
p.addParameter('ZMIN', -50000, @isnumeric);
p.addParameter('use_topo', any(abs(ztopo(:)) > 1e-9), @(x)islogical(x)||isnumeric(x));
p.addParameter('x_st', [], @isnumeric);
p.addParameter('rho_air', 1e8, @isnumeric);
p.addParameter('rho_bg', 100, @isnumeric);
p.addParameter('include_air', false, @(x)islogical(x)||isnumeric(x));
p.addParameter('add_side_robin', false, @(x)islogical(x)||isnumeric(x));
p.addParameter('param_mode', 'cell', @(s)ischar(s)||isstring(s));
p.addParameter('use_phase_in_radian', true, @(x)islogical(x)||isnumeric(x));
p.addParameter('rhs_chunk_size', 100, @isnumeric);
p.parse(varargin{:});
opt = p.Results;
opt.use_topo = logical(opt.use_topo);
opt.include_air = logical(opt.include_air);
opt.add_side_robin = logical(opt.add_side_robin);
opt.use_phase_in_radian = logical(opt.use_phase_in_radian);
opt.param_mode = lower(string(opt.param_mode));

if opt.param_mode ~= "cell"
    error('build_J_TM_direct_sensitivity_cell currently supports only param_mode="cell".');
end
if opt.include_air
    error('This TM direct-sensitivity version is designed for include_air=false, as used in the synthetic TM setup.');
end

mu0 = 4*pi*1e-7;
time_sign = opt.time_sign;
rho_air = opt.rho_air;
rho_bg = opt.rho_bg;
sigma_hs = 1/rho_bg;
F = datafreq;
if ~isvector(F), F = F(:,1); end
F = F(:);
NF = numel(F);
x_st = opt.x_st(:);

cache = build_tm_cache_direct(xtopo, ztopo, NDZ, opt.ZMIN, opt.use_topo, opt.include_air);
NDX = cache.NDX;
NDZ_eff = cache.NDZ_eff;
Nmodel = NDX * NDZ_eff;

m = m_log10rho_cell(:);
if numel(m) ~= Nmodel
    error('Model size mismatch: numel(m)=%d, expected TM cell Nmodel=%d.', numel(m), Nmodel);
end

if isempty(x_st)
    x_st = cache.DX(:);
end
n_st = numel(x_st);
Ndata = 2 * NF * n_st;

% Station nearest grid index
ix_st = zeros(n_st,1);
for ss = 1:n_st
    [~, ix_st(ss)] = min(abs(cache.DX - x_st(ss)));
end

% Model to element rho/sigma
rho_cell = 10.^m;
rho_elem = repelem(rho_cell(:), 4);          % NTE x 1 for include_air=false
sigma_elem = 1 ./ rho_elem;

% sigma below surface and its derivative wrt model cells
[sigma_below_all, dSigBelow_all] = sigma_below_and_derivative_TM(cache, sigma_elem, m);
sig_st = sigma_below_all(ix_st).';           % n_st x 1
Dsig_st_dm = dSigBelow_all(ix_st, :);        % n_st x Nmodel sparse

% Dirichlet/free partition
D = cache.ND;
U = cache.U;
Pdir = ones(numel(D),1);
g2u = cache.g2u;
NN = cache.NN;
NTE = cache.NTE;

% Surface derivative nodes/weights
[der_nodes, der_w] = surface_derivative_weights(cache, ix_st);
k0 = der_nodes(:,1);
row0 = g2u(k0);
row1 = g2u(der_nodes(:,2));
row2 = g2u(der_nodes(:,3));

% Allocate outputs
J = zeros(Ndata, Nmodel, 'double');
Rhoapp  = zeros(NF, n_st);
FasaImp = zeros(NF, n_st);
Impedansi = complex(zeros(NF, n_st));
Ex_store = complex(zeros(NF, n_st));
MatH = complex(zeros(NN, NF));

freq_time = zeros(NF,1);
factor_solve_time = zeros(NF,1);
rhs_assemble_time = zeros(NF,1);
fill_time = zeros(NF,1);

fprintf('\nBuilding TM direct-sensitivity Jacobian per frequency...\n');

for ifr = 1:NF
    tfreq = tic;
    freq = F(ifr);
    omega = 2*pi*freq;
    omega_miu = omega*mu0;
    beta = time_sign * 1i * omega_miu;

    % Matrix before Dirichlet elimination: stiffness(rho) + beta*mass + Robin
    rho_rep = repmat(rho_elem(:), 9, 1);
    Kfull = sparse(cache.I, cache.J, cache.Vstiff .* rho_rep + beta*cache.Vmass, NN, NN);

    kbot = sqrt(beta * sigma_hs);
    if real(kbot) < 0 || (abs(real(kbot)) < 1e-14 && imag(kbot) < 0)
        kbot = -kbot;
    end
    gamma = rho_bg * kbot;
    Kfull = add_robin_line(Kfull, cache.X, cache.Z, cache.idxpl(:,cache.NGZ), gamma);
    if opt.add_side_robin
        Kfull = add_robin_line(Kfull, cache.X, cache.Z, cache.idxpl(1,:).', gamma);
        Kfull = add_robin_line(Kfull, cache.X, cache.Z, cache.idxpl(end,:).', gamma);
    end

    Kuu = Kfull(U,U);
    rhs = -Kfull(U,D) * Pdir;

    tfac = tic;
    Kfac = decomposition(Kuu, 'lu');
    Hu = Kfac \ rhs;
    factor_solve_time(ifr) = toc(tfac);

    H = complex(zeros(NN,1));
    H(U) = Hu;
    H(D) = Pdir;
    MatH(:,ifr) = H;

    % Base response at stations
    H0 = H(k0);
    dHdz0 = der_w(:,1).*H(der_nodes(:,1)) + der_w(:,2).*H(der_nodes(:,2)) + der_w(:,3).*H(der_nodes(:,3));
    Ex0 = -(1 ./ sig_st) .* dHdz0;
    Ztm = Ex0 ./ H0;

    Rhoapp(ifr,:) = ((abs(Ztm).^2) ./ omega_miu).';
    FasaImp(ifr,:) = (mod(angle(Ztm)*180/pi + 180, 360) - 180).';
    Impedansi(ifr,:) = Ztm.';
    Ex_store(ifr,:) = Ex0.';

    % Sensitivity RHS: B = -dK/dm * H, restricted to free rows
    trhs = tic;
    Bsen = assemble_rhs_TM_direct_all_cells(H, rho_elem, cache);
    rhs_assemble_time(ifr) = toc(trhs);

    rho_rows = (ifr-1)*n_st + (1:n_st).';
    phi_rows = NF*n_st + rho_rows;

    tfill = tic;
    for c1 = 1:opt.rhs_chunk_size:Nmodel
        c2 = min(c1 + opt.rhs_chunk_size - 1, Nmodel);
        cols = c1:c2;

        DU = Kfac \ Bsen(:,cols);

        dH0 = get_dH_at_rows(DU, row0);
        dH1 = get_dH_at_rows(DU, row1);
        dH2 = get_dH_at_rows(DU, row2);

        dDdz = der_w(:,1).*dH0 + der_w(:,2).*dH1 + der_w(:,3).*dH2;

        dSigma = full(Dsig_st_dm(:,cols));
        dInvSigma = - dSigma ./ (sig_st.^2);
        dEx = -(dInvSigma .* dHdz0 + (1 ./ sig_st) .* dDdz);

        dZ = (dEx .* H0 - Ex0 .* dH0) ./ (H0.^2);

        dlogrho = (2/log(10)) * real(dZ ./ Ztm);
        dphi_rad = imag(dZ ./ Ztm);
        if opt.use_phase_in_radian
            dphi = dphi_rad;
        else
            dphi = dphi_rad * 180/pi;
        end

        J(rho_rows, cols) = dlogrho;
        J(phi_rows, cols) = dphi;
    end
    fill_time(ifr) = toc(tfill);
    freq_time(ifr) = toc(tfreq);

    fprintf('  freq %3d/%3d | f=%10.4g Hz | %.2f s\n', ifr, NF, freq, freq_time(ifr));
end

out0 = struct();
out0.Rhoapp = Rhoapp;
out0.PhiDeg = FasaImp;
out0.Impedansi = Impedansi;
out0.Ex_store = Ex_store;
out0.MatH = MatH;
out0.DX = cache.DX;
out0.DZ = cache.DZ;
out0.X = cache.X;
out0.Z = cache.Z;
out0.Zref = cache.Zref;
out0.NE = cache.NE;
out0.idxpl = cache.idxpl;
out0.x_st = x_st;
out0.ix_st = ix_st;
out0.model_id = cache.model_id;
out0.NGX = cache.NGX;
out0.NGZ = cache.NGZ;
out0.NTE = cache.NTE;
out0.NN = cache.NN;
out0.NDX = cache.NDX;
out0.j0 = cache.j0;
out0.include_air = opt.include_air;
out0.use_topo = opt.use_topo;
out0.time_sign = time_sign;
out0.ZMIN = opt.ZMIN;
out0.rho_bg = rho_bg;
out0.rho_air = rho_air;
out0.sigma = sigma_elem(:);
out0.sigma_below = sigma_below_all(:);

timing = struct();
timing.freq_time = freq_time;
timing.factor_solve_time = factor_solve_time;
timing.rhs_assemble_time = rhs_assemble_time;
timing.fill_time = fill_time;
timing.total_freq_time = sum(freq_time);
timing.total_factor_solve_time = sum(factor_solve_time);
timing.total_rhs_assemble_time = sum(rhs_assemble_time);
timing.total_fill_time = sum(fill_time);

meta = struct();
meta.mode = 'TM';
meta.jacobian_type = 'direct_sensitivity';
meta.param_mode = 'cell';
meta.use_topo = opt.use_topo;
meta.include_air = opt.include_air;
meta.add_side_robin = opt.add_side_robin;
meta.NDZ = NDZ;
meta.NDZ_eff = cache.NDZ_eff;
meta.ZMIN = opt.ZMIN;
meta.NDX = NDX;
meta.Nmodel = Nmodel;
meta.Ndata = Ndata;
meta.NF = NF;
meta.n_station = n_st;
meta.time_sign = time_sign;
meta.use_phase_in_radian = opt.use_phase_in_radian;
meta.rho_air = rho_air;
meta.rho_bg = rho_bg;
meta.xtopo = xtopo(:);
meta.ztopo = ztopo(:);
meta.x_st = x_st(:);
meta.F = F(:);
meta.rhs_chunk_size = opt.rhs_chunk_size;
meta.timing = timing;
end

%% ========================================================================
function cache = build_tm_cache_direct(xtopo, ztopo, NDZ, ZMIN, use_topo, include_air)
xtopo = xtopo(:); ztopo = ztopo(:);
DX = xtopo(:);
NDX = numel(DX)-1;

if include_air
    error('build_tm_cache_direct: include_air=true is not enabled for direct-J TM.');
else
    DZ = safe_funcDZ(NDZ);
    DZ = DZ(:);
    if abs(DZ(1)) > 1e-12
        error('DZ(1) must be 0 for TM without air.');
    end
    if any(diff(DZ) <= 0)
        error('DZ must increase monotonically.');
    end
    NDZ_eff = numel(DZ)-1;
    NDZa = 0;
    NDZb_eff = NDZ_eff;
end

NGX = NDX + 1;
NGZ = numel(DZ);
NTE = NDX * NDZ_eff * 4;
NN  = (NGX * NGZ) + (NDX * NDZ_eff);

% Nodes
k = 0;
X = zeros(NN,1); Z = zeros(NN,1); idxpl = zeros(NGX,NGZ);
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

% Connectivity
NE = zeros(3,NTE);
k = 0; l = 0;
for j = 1:NDZ_eff
    for i = 1:NDX
        k = k+1; l = l+1;
        N1 = k; N2 = k+1; N3 = k+NGX; N4 = k+NGX+NDX; N5 = N4+1;
        mm = (l-1)*4;
        NE(:,mm+1) = [N1; N2; N3];
        NE(:,mm+2) = [N2; N5; N3];
        NE(:,mm+3) = [N5; N4; N3];
        NE(:,mm+4) = [N4; N1; N3];
    end
    k = k + NGX;
end

Zref = Z;
if use_topo
    zt = interp1(xtopo(:), ztopo(:), X, 'makima', 'extrap');
else
    zt = zeros(size(X));
end
ZBOT = max(Zref);
w = 1 - (Zref / ZBOT);
w = max(0,min(1,w));
Z = Zref + w.*zt;

% Geometry and sparse triplets
n1 = NE(1,:).'; n2 = NE(2,:).'; n3 = NE(3,:).';
be1 = Z(n2)-Z(n3); be2 = Z(n3)-Z(n1); be3 = Z(n1)-Z(n2);
ce1 = X(n3)-X(n2); ce2 = X(n1)-X(n3); ce3 = X(n2)-X(n1);
A_all = abs(0.5*(be1.*ce2 - be2.*ce1));
if any(A_all < 1e-20), error('A degenerate element was found.'); end
Bmat = [be1 be2 be3];
Cmat = [ce1 ce2 ce3];
nodes = [n1 n2 n3];

I = zeros(9*NTE,1); J = zeros(9*NTE,1); Vstiff = zeros(9*NTE,1); Vmass = zeros(9*NTE,1);
idx = 0;
for a = 1:3
    for b = 1:3
        idx = idx + 1;
        s = (idx-1)*NTE + (1:NTE);
        I(s) = nodes(:,a);
        J(s) = nodes(:,b);
        Vstiff(s) = (Bmat(:,a).*Bmat(:,b) + Cmat(:,a).*Cmat(:,b)) ./ (4*A_all);
        Vmass(s) = (1 + (a==b)) .* (A_all/12);
    end
end

j0 = find(abs(DZ) < 1e-12,1,'first');
if isempty(j0), error('DZ=0 was not found.'); end
if j0 >= NGZ, error('The surface node is located at the lower grid boundary.'); end

ND = idxpl(:,1);     % top boundary Dirichlet
allnodes = (1:NN).';
U = setdiff(allnodes, ND);
g2u = zeros(NN,1);
g2u(U) = 1:numel(U);

model_id = (1:NTE).';
cell2tri = cell(NDX*NDZ_eff,1);
for iz = 1:NDZ_eff
    for ix = 1:NDX
        cid = (iz-1)*NDX + ix;
        box_id = (iz-1)*NDX + ix;
        cell2tri{cid} = (box_id-1)*4 + (1:4);
    end
end

cache = struct();
cache.DX = DX; cache.DZ = DZ; cache.NDX = NDX; cache.NGX = NGX; cache.NGZ = NGZ;
cache.NDZ_eff = NDZ_eff; cache.NDZa = NDZa; cache.NDZb_eff = NDZb_eff;
cache.NTE = NTE; cache.NN = NN; cache.X = X; cache.Z = Z; cache.Zref = Zref;
cache.idxpl = idxpl; cache.NE = NE; cache.A_all = A_all(:); cache.model_id = model_id;
cache.I = I; cache.J = J; cache.Vstiff = Vstiff; cache.Vmass = Vmass;
cache.j0 = j0; cache.ND = ND; cache.U = U; cache.g2u = g2u; cache.cell2tri = cell2tri;
end

%% ========================================================================
function [sigma_below, dSigBelow] = sigma_below_and_derivative_TM(cache, sigma_elem, m_cell)
NDX = cache.NDX; NGX = cache.NGX; j0 = cache.j0;
Nmodel = NDX*cache.NDZ_eff;
if j0 ~= 1
    error('The current direct TM sigma_below derivative assumes include_air=false, so j0=1.');
end

sigma_cell = zeros(1,NDX);
for i = 1:NDX
    l = (j0-1)*NDX + i;
    tri = (l-1)*4 + (1:4);
    sigma_cell(i) = mean(sigma_elem(tri));
end
sigma_below = zeros(1,NGX);
for i = 1:NGX
    if i == 1
        sigma_below(i) = sigma_cell(1);
    elseif i == NGX
        sigma_below(i) = sigma_cell(end);
    else
        sigma_below(i) = 0.5*(sigma_cell(i-1)+sigma_cell(i));
    end
end

% derivative wrt m_cell for surface row only
rows = []; cols = []; vals = [];
for i = 1:NGX
    if i == 1
        c = 1; w = 1;
        rows(end+1) = i; cols(end+1) = c; vals(end+1) = w*(-log(10))*sigma_cell(c);
    elseif i == NGX
        c = NDX; w = 1;
        rows(end+1) = i; cols(end+1) = c; vals(end+1) = w*(-log(10))*sigma_cell(c);
    else
        c1 = i-1; c2 = i;
        rows(end+1:end+2) = [i i];
        cols(end+1:end+2) = [c1 c2];
        vals(end+1:end+2) = 0.5*(-log(10))*[sigma_cell(c1) sigma_cell(c2)];
    end
end
dSigBelow = sparse(rows, cols, vals, NGX, Nmodel);
end

%% ========================================================================
function B = assemble_rhs_TM_direct_all_cells(H, rho_elem, cache)
% B(:,ic) = -dK/dm_ic * H, restricted to free rows.
Nmodel = cache.NDX * cache.NDZ_eff;
NTE = cache.NTE;
g2u = cache.g2u;
maxnnz = 36 * Nmodel;
Irow = zeros(maxnnz,1); Jcol = zeros(maxnnz,1); Vval = complex(zeros(maxnnz,1));
ptr = 0;

for ic = 1:Nmodel
    tri_ids = cache.cell2tri{ic};
    for tt = 1:numel(tri_ids)
        e = tri_ids(tt);
        drho_dm = log(10) * rho_elem(e);
        idxs = e + (0:8)'*NTE;
        rowsG = cache.I(idxs);
        colsG = cache.J(idxs);
        vals = -cache.Vstiff(idxs) .* drho_dm .* H(colsG);
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
B = sparse(Irow(1:ptr), Jcol(1:ptr), Vval(1:ptr), numel(cache.U), Nmodel);
end

%% ========================================================================
function [nodes, weights] = surface_derivative_weights(cache, ix_st)
n = numel(ix_st);
nodes = zeros(n,3);
weights = zeros(n,3);
for s = 1:n
    i = ix_st(s);
    k0 = cache.idxpl(i, cache.j0);
    if (cache.j0 + 2) <= cache.NGZ
        k1 = cache.idxpl(i, cache.j0+1);
        k2 = cache.idxpl(i, cache.j0+2);
        z0 = cache.Z(k0); z1 = cache.Z(k1); z2 = cache.Z(k2);
        w0 = deriv2_weights_at_surface(z0,z1,z2,1);
        w1 = deriv2_weights_at_surface(z0,z1,z2,2);
        w2 = deriv2_weights_at_surface(z0,z1,z2,3);
        nodes(s,:) = [k0 k1 k2];
        weights(s,:) = [w0 w1 w2];
    else
        k1 = cache.idxpl(i, cache.j0+1);
        dz = cache.Z(k1) - cache.Z(k0);
        nodes(s,:) = [k0 k1 k1];
        weights(s,:) = [-1/dz 1/dz 0];
    end
end
end

function w = deriv2_weights_at_surface(z0,z1,z2,whichH)
dz1 = z1-z0; dz2 = z2-z0;
K = [dz1, dz1^2/2; dz2, dz2^2/2];
if whichH == 1
    Y = [-1; -1];
elseif whichH == 2
    Y = [1; 0];
else
    Y = [0; 1];
end
coef = K \ Y;
w = coef(1);
end

function dH = get_dH_at_rows(DU, rows)
n = numel(rows); nc = size(DU,2);
dH = zeros(n,nc);
mask = rows > 0;
if any(mask)
    dH(mask,:) = DU(rows(mask),:);
end
end

function XK = add_robin_line(XK, X, Z, line_nodes, gamma)
line_nodes = line_nodes(:);
for ee = 1:(numel(line_nodes)-1)
    n1 = line_nodes(ee); n2 = line_nodes(ee+1);
    Ledge = hypot(X(n2)-X(n1), Z(n2)-Z(n1));
    base = gamma * Ledge / 6;
    XK(n1,n1) = XK(n1,n1) + 2*base;
    XK(n1,n2) = XK(n1,n2) + 1*base;
    XK(n2,n1) = XK(n2,n1) + 1*base;
    XK(n2,n2) = XK(n2,n2) + 2*base;
end
end

function DZ = safe_funcDZ(NDZ)
% For joint TE-TM inversion, the TM earth grid must match the shared TE earth grid.
% Therefore, funcDZearth_shared is used with priority.
if exist('funcDZearth_shared','file') == 2
    DZ = funcDZearth_shared(NDZ);
elseif exist('funcDZ','file') == 2
    DZ = funcDZ(NDZ);
else
    ZMAX = 550000;
    zpos = logspace(log10(25), log10(ZMAX), NDZ).';
    zpos(end) = ZMAX;
    DZ = [0; zpos];
end
DZ = DZ(:);
end
