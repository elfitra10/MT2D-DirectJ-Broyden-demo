function [rhoa_st, phi_st_deg, out] = fungsi_TM_mod_shared(m_log10rho, xtopo, ztopo, datafreq, NDZ, varargin)
% fungsi_TM_mod_shared
% Forward 2-D MT TM-mode modelling (Hy formulation, triangular FEM).
% MAIN INPUTS
%   m_log10rho : log10(rho), size depends on param_mode
%   xtopo      : surface x-node coordinates (NGX x 1)
%   ztopo      : topographic elevation at xtopo (NGX x 1)
%   datafreq   : frequency vector or matrix; column 1 is used
%   NDZ        : number of vertical earth cells, excluding air cells
%
% OPTIONAL ARGUMENTS ('key',value)
%   'time_sign'      : +1 (default) or -1
%   'ZMIN'           : upper air boundary if include_air=true (default -50000)
%   'use_topo'       : true/false (default is inferred from ztopo)
%   'x_st'           : station positions (default: passed station coordinates, if available)
%   'rho_air'        : air resistivity (default 1e8 ohm.m)
%   'rho_bg'         : bottom half-space resistivity for the Robin boundary condition (default 100)
%   'include_air'    : true/false (default false)
%   'add_side_robin' : true/false (default false)


% ===================== PARSE OPTIONS =====================
opt = struct();
opt.time_sign      = 1;
opt.ZMIN           = -50000;
opt.use_topo       = any(abs(ztopo(:)) > 1e-9);
opt.x_st           = [];
opt.rho_air        = 1e8;
opt.rho_bg         = 100;
opt.include_air    = false;
opt.add_side_robin = false;
opt.param_mode     = 'auto';

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('Optional arguments must be provided as paired ''key'', value entries.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k+1};
        switch key
            case "time_sign",      opt.time_sign = val;
            case "zmin",           opt.ZMIN = val;
            case "use_topo",       opt.use_topo = logical(val);
            case "x_st",           opt.x_st = val;
            case "rho_air",        opt.rho_air = val;
            case {"rho_bg","rho_hs"}, opt.rho_bg = val;
            case "include_air",    opt.include_air = logical(val);
            case "add_side_robin", opt.add_side_robin = logical(val);
            case "param_mode",     opt.param_mode = lower(string(val));
            otherwise
                % kept for backward compatibility with legacy callers
        end
    end
end

% ===================== CONSTANTS =====================
mu0       = 4*pi*1e-7;
time_sign = opt.time_sign;
rho_air   = opt.rho_air;
rho_bg    = opt.rho_bg;
sigma_air = 1/rho_air;
sigma_hs  = 1/rho_bg;
use_topo    = opt.use_topo;
include_air = opt.include_air;

F = datafreq;
if ~isvector(F), F = F(:,1); end
F = F(:);
NF = numel(F);

% ===================== MESH CACHE =====================
persistent cache
need_rebuild = true;

if ~isempty(cache)
    try
        need_rebuild = (cache.NDZ_input ~= NDZ) || ...
                       ~isequal(cache.xtopo, xtopo(:)) || ...
                       ~isequal(cache.ztopo, ztopo(:)) || ...
                       (cache.use_topo ~= use_topo) || ...
                       (cache.include_air ~= include_air) || ...
                       (cache.ZMIN ~= opt.ZMIN);
    catch
        need_rebuild = true;
    end
end

if need_rebuild
    cache = struct();
    cache.xtopo       = xtopo(:);
    cache.ztopo       = ztopo(:);
    cache.NDZ_input   = NDZ;
    cache.use_topo    = use_topo;
    cache.include_air = include_air;
    cache.ZMIN        = opt.ZMIN;

    DX  = xtopo(:);
    NDX = numel(DX)-1;

    if include_air
        DZb = funcDZearth_shared(NDZ);
        DZb = DZb(:);
        if abs(DZb(1)) > 1e-12
            error('funcDZ(NDZ) must return DZ starting from 0.');
        end

        ZMIN = opt.ZMIN;
        if ZMIN >= 0
            error('ZMIN must be negative when include_air=true.');
        end

        a = 1; 
        b = 3;
        NDZa_guess = max(1, round(a*NDZ/(a+b)));
        dz1 = max(1, DZb(min(2,end)) - DZb(1));
        DZa = logspace(log10(dz1), log10(abs(ZMIN)), NDZa_guess).';
        DZa(end) = abs(ZMIN);

        DZ = [-flipud(cumsum(DZa)); DZb];
        DZ(abs(DZ) < 1e-12) = 0;
        iz0 = find(abs(DZ) < 1e-12);
        if numel(iz0) ~= 1
            DZ = unique(DZ,'stable');
        end

        NDZ_eff = numel(DZ)-1;
        NDZa = find(DZ < 0, 1, 'last');
        if isempty(NDZa), NDZa = 0; end
        NDZb_eff = NDZ_eff - NDZa;
    else
        DZ = funcDZearth_shared(NDZ);
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

    % ===================== NODE COORDINATES =====================
    k = 0;
    X = zeros(NN,1);
    Z = zeros(NN,1);
    idxpl = zeros(NGX, NGZ);

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

    % ===================== CONNECTIVITY =====================
    NE = zeros(3, NTE);
    k = 0; l = 0;
    for j = 1:NDZ_eff
        for i = 1:NDX
            k = k+1;
            l = l+1;
            N1 = k;
            N2 = k+1;
            N3 = k+NGX;
            N4 = k+NGX+NDX;
            N5 = N4+1;

            m = (l-1)*4;
            NE(:,m+1) = [N1; N2; N3];
            NE(:,m+2) = [N2; N5; N3];
            NE(:,m+3) = [N5; N4; N3];
            NE(:,m+4) = [N4; N1; N3];
        end
        k = k + NGX;
    end

    Zref = Z;

    if use_topo
        zt = interp1(xtopo(:), ztopo(:), X, 'makima', 'extrap');
    else
        zt = zeros(size(X));
    end

    if include_air
        ZTOP = min(Zref);
        ZBOT = max(Zref);
        w = zeros(size(Zref));
        idxE = (Zref >= 0);
        w(idxE) = 1 - (Zref(idxE) / ZBOT);
        idxA = (Zref < 0);
        w(idxA) = (Zref(idxA) - ZTOP) / (0 - ZTOP);
    else
        ZBOT = max(Zref);
        w = 1 - (Zref / ZBOT);
    end
    w = max(0, min(1, w));
    Z = Zref + w .* zt;

    % ===================== GEOMETRY =====================
    n1 = NE(1,:).'; n2 = NE(2,:).'; n3 = NE(3,:).';
    be1 = Z(n2) - Z(n3);
    be2 = Z(n3) - Z(n1);
    be3 = Z(n1) - Z(n2);
    ce1 = X(n3) - X(n2);
    ce2 = X(n1) - X(n3);
    ce3 = X(n2) - X(n1);

    A_all = abs(0.5*(be1.*ce2 - be2.*ce1));
    if any(A_all < 1e-20)
        error('A degenerate element was found (area ~ 0).');
    end

    be_all = [be1.'; be2.'; be3.'];
    ce_all = [ce1.'; ce2.'; ce3.'];

    % ===================== EARTH-MODEL ELEMENT INDICES =====================
    zc_ref = (Zref(n1) + Zref(n2) + Zref(n3))/3;
    if include_air
        model_id = find(zc_ref > 0);
    else
        model_id = (1:NTE).';
    end

    % ===================== SURFACE INDEX =====================
    j0 = find(abs(DZ) < 1e-12, 1, 'first');
    if isempty(j0)
        error('DZ=0 (surface) was not found.');
    end
    if j0 >= NGZ
        error('The surface node is located at the lower grid boundary.');
    end

    % ===================== DEFAULT STATIONS =====================
    basePath = fileparts(mfilename('fullpath'));
    if use_topo
        cand = {fullfile(basePath,'x_st_topo.txt'), fullfile(basePath,'x_st_topo_Patuha.txt')};
    else
        cand = {fullfile(basePath,'x_st_datar.txt'), fullfile(basePath,'x_st_sint_datar.txt')};
    end

    xst_def = DX;
    for ic = 1:numel(cand)
        if isfile(cand{ic})
            xst_def = load(cand{ic});
            break;
        end
    end
    xst_def = xst_def(:);

    ix_def = zeros(numel(xst_def),1);
    for ss = 1:numel(xst_def)
        [~, ix_def(ss)] = min(abs(DX - xst_def(ss)));
    end

    % ===================== SAVE CACHE =====================
    cache.DX = DX;
    cache.DZ = DZ;
    cache.NDX = NDX;
    cache.NGX = NGX;
    cache.NGZ = NGZ;
    cache.NDZ_eff = NDZ_eff;
    cache.NDZa = NDZa;
    cache.NDZb_eff = NDZb_eff;
    cache.NTE = NTE;
    cache.NN  = NN;
    cache.X = X;
    cache.Z = Z;
    cache.Zref = Zref;
    cache.idxpl = idxpl;
    cache.NE = NE;
    cache.be_all = be_all;
    cache.ce_all = ce_all;
    cache.A_all = A_all;
    cache.model_id = model_id;
    cache.j0 = j0;
    cache.x_st_def = xst_def;
    cache.ix_st_def = ix_def;
end

% ===================== UNPACK CACHE =====================
DX       = cache.DX;
DZ       = cache.DZ;
X        = cache.X;
Z        = cache.Z;
Zref     = cache.Zref;
NE       = cache.NE;
idxpl    = cache.idxpl;
be_all   = cache.be_all;
ce_all   = cache.ce_all;
A_all    = cache.A_all;
model_id = cache.model_id;
j0       = cache.j0;
NGX      = cache.NGX;
NGZ      = cache.NGZ;
NTE      = cache.NTE;
NN       = cache.NN;
NDX      = cache.NDX;
NDZ_eff  = cache.NDZ_eff;

% ===================== MODEL -> SIGMA =====================
if isempty(m_log10rho)
    error('m_log10rho is empty.');
end
m = m_log10rho(:);

if include_air
    NDZb_eff = cache.NDZb_eff;
    nCellEarth = NDX * NDZb_eff;
else
    nCellEarth = NDX * NDZ_eff;
end
nTriEarth = numel(model_id);

mode = opt.param_mode;
if mode == "auto"
    if numel(m) == nCellEarth
        mode = "cell";
    elseif numel(m) == nTriEarth
        mode = "tri";
    else
        error('Model length mismatch: numel(m)=%d, expected %d (cell) or %d (tri).', numel(m), nCellEarth, nTriEarth);
    end
end

if mode == "cell"
    rho_cell = 10.^m;
    rho_earth_tri = expand_cell_to_tri(rho_cell, NDX, NDZ_eff, DZ, include_air);
    if numel(rho_earth_tri) ~= nTriEarth
        error('Cell-to-triangle expansion failed: tri=%d, expected=%d.', numel(rho_earth_tri), nTriEarth);
    end
else
    rho_earth_tri = 10.^m;
end

sigma = sigma_air * ones(NTE,1);
if ~include_air
    sigma(:) = 1;
end
sigma(model_id) = 1 ./ rho_earth_tri;

% average conductivity immediately below the surface for Ex
sigma_cell = zeros(1, NDX);
for i = 1:NDX
    l = (j0-1)*NDX + i;
    tri = (l-1)*4 + (1:4);
    sigma_cell(i) = mean(sigma(tri));
end

sigma_below = zeros(1, NGX);
for i = 1:NGX
    if i == 1
        sigma_below(i) = sigma_cell(1);
    elseif i == NGX
        sigma_below(i) = sigma_cell(end);
    else
        sigma_below(i) = 0.5 * (sigma_cell(i-1) + sigma_cell(i));
    end
end

% ===================== STATION INDICES =====================
if ~isempty(opt.x_st)
    xst = opt.x_st(:);
    ix = zeros(numel(xst),1);
    for ss = 1:numel(xst)
        [~, ix(ss)] = min(abs(DX - xst(ss)));
    end
else
    xst = cache.x_st_def;
    ix  = cache.ix_st_def;
end

% ===================== SOLVE PER FREQUENCY =====================
MatH      = complex(zeros(NN, NF));
Impedansi = complex(zeros(NF, NGX));
Rhoapp    = zeros(NF, NGX);
FasaImp   = zeros(NF, NGX);
Ex_store  = complex(zeros(NF, NGX));

ND = idxpl(:,1);        % top-most boundary nodes
P  = ones(numel(ND),1); % Hy = 1

for ifr = 1:NF
    freq = F(ifr);
    omega = 2*pi*freq;
    omega_miu = omega*mu0;
    beta = time_sign * 1i * omega_miu;

    rho_elem = zeros(NTE,1);
    rho_elem(model_id) = 1 ./ sigma(model_id);
    if include_air
        rho_elem(zref_centroid(NE,Zref) < 0) = rho_air;
    end

    est_nnz = 9*NTE + 8*(NGX+NGZ);
    ii = zeros(est_nnz,1);
    jj = zeros(est_nnz,1);
    vv = complex(zeros(est_nnz,1));
    nnz_cnt = 0;

    for e = 1:NTE
        loc = NE(:,e);
        be = be_all(:,e);
        ce = ce_all(:,e);
        A  = A_all(e);
        rhoe = rho_elem(e);

        for a = 1:3
            na = loc(a);
            for b = a:3
                nb = loc(b);
                Ke = rhoe * (be(a)*be(b) + ce(a)*ce(b)) / (4*A);
                Me = beta * (1 + (a==b)) * (A/12);
                val = Ke + Me;

                nnz_cnt = nnz_cnt + 1;
                ii(nnz_cnt) = na; jj(nnz_cnt) = nb; vv(nnz_cnt) = val;
                if a ~= b
                    nnz_cnt = nnz_cnt + 1;
                    ii(nnz_cnt) = nb; jj(nnz_cnt) = na; vv(nnz_cnt) = val;
                end
            end
        end
    end

    XK = sparse(ii(1:nnz_cnt), jj(1:nnz_cnt), vv(1:nnz_cnt), NN, NN);
    B  = complex(zeros(NN,1));

    % bottom Robin boundary condition
    kbot = sqrt(beta * sigma_hs);
    if real(kbot) < 0 || (abs(real(kbot))<1e-14 && imag(kbot)<0)
        kbot = -kbot;
    end
    gamma = rho_bg * kbot;

    bot_nodes = idxpl(:,NGZ);
    XK = add_robin_line(XK, X, Z, bot_nodes, gamma);

    if opt.add_side_robin
        left_nodes  = idxpl(1,:).';
        right_nodes = idxpl(end,:).';
        XK = add_robin_line(XK, X, Z, left_nodes,  gamma);
        XK = add_robin_line(XK, X, Z, right_nodes, gamma);
    end

    % Dirichlet condition at the top boundary
    B = B - XK(:,ND)*P;
    B(ND) = P;
    XK(:,ND) = 0;
    XK(ND,:) = 0;
    XK = XK + sparse(ND, ND, ones(numel(ND),1), NN, NN);

    H = XK \ B;
    MatH(:,ifr) = H;

    Z_tm    = complex(zeros(1,NGX));
    rhoa_tm = zeros(1,NGX);
    phi_tm  = zeros(1,NGX);
    Ex_line = complex(zeros(1,NGX));

    for i = 1:NGX
        k0 = idxpl(i,j0);
        if (j0+2) <= NGZ
            k1 = idxpl(i,j0+1);
            k2 = idxpl(i,j0+2);
            dH_dz = deriv2_at_surface(H(k0), H(k1), H(k2), Z(k0), Z(k1), Z(k2));
        else
            k1 = idxpl(i,j0+1);
            dH_dz = (H(k1)-H(k0)) / (Z(k1)-Z(k0));
        end

        Ex = -(1/sigma_below(i)) * dH_dz;     % FIX: sign Ex
        Ex_line(i) = Ex;
        Z_tm(i)   = Ex / H(k0);
        rhoa_tm(i)= (abs(Z_tm(i))^2) / omega_miu;
        
        phi_tm(i) = angle(Z_tm(i)) * 180/pi;

        % Optional: keep phase stable in the [-180, 180] interval
        phi_tm(i) = mod(phi_tm(i)+180, 360) - 180;
    end

    Impedansi(ifr,:) = Z_tm;
    Rhoapp(ifr,:)    = rhoa_tm;
    FasaImp(ifr,:)   = phi_tm;
    Ex_store(ifr,:)  = Ex_line;
end

rhoa_st    = Rhoapp(:, ix);
phi_st_deg = FasaImp(:, ix);

out = struct();
out.Rhoapp      = Rhoapp;
out.PhiDeg      = FasaImp;
out.Impedansi   = Impedansi;
out.Ex_store    = Ex_store;
out.MatH        = MatH;
out.DX          = DX;
out.DZ          = DZ;
out.X           = X;
out.Z           = Z;
out.Zref        = Zref;
out.NE          = NE;
out.idxpl       = idxpl;
out.x_st        = xst;
out.ix_st       = ix;
out.model_id    = model_id;
out.NGX         = NGX;
out.NGZ         = NGZ;
out.NTE         = NTE;
out.NN          = NN;
out.NDX         = NDX;
out.j0          = j0;
out.include_air = include_air;
out.use_topo    = use_topo;
out.time_sign   = time_sign;
out.ZMIN        = opt.ZMIN;
out.rho_bg      = rho_bg;
out.rho_air     = rho_air;
out.sigma       = sigma(:);
out.sigma_below = sigma_below(:);
end

%% ================= LOCAL FUNCTIONS =================
function DZ = safe_funcDZ(NDZ)
if exist('funcDZ','file') == 2
    DZ = funcDZ(NDZ);
else
    ZMAX = 550000;
    zpos = logspace(log10(25), log10(ZMAX), NDZ).';
    zpos(end) = ZMAX;
    DZ = [0; zpos];
end
DZ = DZ(:);
end

function XK = add_robin_line(XK, X, Z, line_nodes, gamma)
for ee = 1:(numel(line_nodes)-1)
    n1 = line_nodes(ee);
    n2 = line_nodes(ee+1);
    Ledge = hypot(X(n2)-X(n1), Z(n2)-Z(n1));
    base  = gamma * Ledge / 6;
    XK(n1,n1) = XK(n1,n1) + 2*base;
    XK(n1,n2) = XK(n1,n2) + 1*base;
    XK(n2,n1) = XK(n2,n1) + 1*base;
    XK(n2,n2) = XK(n2,n2) + 2*base;
end
end

function dHdz = deriv2_at_surface(H0,H1,H2,z0,z1,z2)
dz1 = z1 - z0;
dz2 = z2 - z0;
Y = [H1-H0; H2-H0];
K = [dz1, dz1^2/2; dz2, dz2^2/2];
coef = K\Y;
dHdz = coef(1);
end

function zc = zref_centroid(NE,Zref)
zc = (Zref(NE(1,:)) + Zref(NE(2,:)) + Zref(NE(3,:)))/3;
zc = zc(:);
end

function rho_tri = expand_cell_to_tri(rho_cell, NDX, NDZ_eff, DZ, include_air)
if ~include_air
    rho_tri = repelem(rho_cell(:), 4);
    return;
end

j0 = find(abs(DZ)<1e-12,1,'first');
NDZa = j0 - 1;
NDZb = (numel(DZ)-1) - NDZa;
if numel(rho_cell) ~= NDX*NDZb
    error('For include_air=true, rho_cell must have size NDX*NDZb.');
end
rho_tri = zeros(4*NDX*NDZb,1);
idx = 0;
for jb = 1:NDZb
    for i = 1:NDX
        idx = idx + 1;
        rhoe = rho_cell((jb-1)*NDX + i);
        rho_tri((idx-1)*4 + (1:4)) = rhoe;
    end
end
end
