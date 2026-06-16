% Clean release copy. TE forward solver; station coordinates are passed through x_st.

function [rhoa_st, phi_st_deg, out] = fungsi_TE_mod_cell_shared(m_log10rho_cell, xtopo, ztopo, datafreq, NDZ, varargin)
% fungsi_TE_mod_cell_shared
% Release version: station coordinates are passed through 'x_st', x_st;
% this function no longer reads x_st_Sint_topo.txt or x_st_Sint_datar.txt internally.
% Model input: one log10(rho) value per earth cell.
%
% Output remains compatible with the inversion and plotting routines:


% ===================== PARSE OPTIONAL ARGS =====================
opt = struct();
opt.time_sign = 1;
opt.ZMIN      = -50000;
opt.use_topo  = any(abs(ztopo(:)) > 1e-9);
opt.x_st      = [];
opt.rho_air   = 1e8;
opt.rho_bg    = 100;

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('Optional arguments must be provided as paired ''key'', value entries.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k+1};
        switch key
            case "time_sign"
                opt.time_sign = val;
            case "zmin"
                opt.ZMIN = val;
            case "use_topo"
                opt.use_topo = logical(val);
            case "x_st"
                opt.x_st = val;
            case "rho_air"
                opt.rho_air = val;
            case "rho_bg"
                opt.rho_bg = val;
            case "rho_hs"
                opt.rho_bg = val;
        end
    end
end

% ===================== CONSTANTS AND INPUTS =====================
time_sign = opt.time_sign;
mu0 = 4*pi*1e-7;

rho_air   = opt.rho_air;
rho_bg    = opt.rho_bg;
sigma_air = 1/rho_air;
sigma_hs  = 1/rho_bg;

F = datafreq(:,1);
F = F(:);
NF = numel(F);

use_topo = opt.use_topo;

% ===================== CACHE =====================
persistent cache
need_rebuild = true;

required_fields = {'NDX','NDZ','NDZa','NDZb','NGX','NGZ','NTE','NN', ...
                   'DX','DZ','X','Z','Zref','idxpl','NE','model_id', ...
                   'cell2tri','I','J','Vkg','mass_base', ...
                   'Kg_uu','Kg_ud','Rb_uu','Rb_ud', ...
                   'D','U','Pdir','kS','kU','kD', ...
                   'x_st_def','ix_st_def', ...
                   'xtopo','ztopo','use_topo','ZMIN'};

    if ~isempty(cache)
        try
            has_all_fields = all(isfield(cache, required_fields));
    
            if ~has_all_fields
                need_rebuild = true;
            else
                need_rebuild = (cache.NDZ ~= NDZ) || ...
                               ~isequal(cache.xtopo, xtopo(:)) || ...
                               ~isequal(cache.ztopo, ztopo(:)) || ...
                               (cache.use_topo ~= use_topo) || ...
                               (cache.ZMIN ~= opt.ZMIN);
            end
        catch
            need_rebuild = true;
        end
    end

if need_rebuild
    cache = struct();
    cache.xtopo    = xtopo(:);
    cache.ztopo    = ztopo(:);
    cache.NDZ      = NDZ;
    cache.use_topo = use_topo;
    cache.ZMIN     = opt.ZMIN;

    ZMIN = opt.ZMIN;

    NDX = length(xtopo)-1;
    a = 1; b = 3;
    NDZa = round(a*NDZ/(a+b));
    NDZb = NDZ-NDZa;

    NGX = NDX+1;
    NGZ = NDZ+1;

    NTE = NDX*NDZ*4;
    NN  = (NGX*NGZ)+(NDX*NDZ);

    DX = xtopo(:);

    % Air layer
    DZa = logspace(log10(25), log10(abs(ZMIN)), NDZa).';
    DZa(end) = abs(ZMIN);

   % earth layer using the shared grid
    DZb = funcDZearth_shared(NDZb);
    DZb = DZb(:);
    while ~isempty(DZb) && abs(DZb(1)) < 1e-12
        DZb(1) = [];
    end
    DZ = [-flipud(DZa); 0; DZb];

    % nodes
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

    % midpoints
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

    % Topographic shift
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

    % Connectivity
    NE = konektivitas(NTE, NGX, NDX, NDZ);

    % earth triangles by centroid in reference grid
    zc = (Zref(NE(1,:)) + Zref(NE(2,:)) + Zref(NE(3,:))) / 3;
    model_id = find(zc > 0);

    % Geometry
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
    
    % store geometry coefficients in three-column matrices
    Bmat = [be1 be2 be3];
    Cmat = [ce1 ce2 ce3];
    
    % sparse assembly triplets
    I = zeros(9*NTE,1);
    J = zeros(9*NTE,1);
    Vkg = zeros(9*NTE,1);
    mass_base = zeros(9*NTE,1);
    
    nodes = [n1 n2 n3];
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

    % Robin bottom
    bot_nodes = idxpl(:,NGZ);
    Irb = zeros(4*(NGX-1),1);
    Jrb = zeros(4*(NGX-1),1);
    Vrb = zeros(4*(NGX-1),1);
    kk = 0;
    for ee = 1:(NGX-1)
        nb1 = bot_nodes(ee);
        nb2 = bot_nodes(ee+1);

        Ledge = hypot(X(nb2)-X(nb1), Z(nb2)-Z(nb1));
        base  = Ledge/6;

        kk=kk+1; Irb(kk)=nb1; Jrb(kk)=nb1; Vrb(kk)=2*base;
        kk=kk+1; Irb(kk)=nb1; Jrb(kk)=nb2; Vrb(kk)=1*base;
        kk=kk+1; Irb(kk)=nb2; Jrb(kk)=nb1; Vrb(kk)=1*base;
        kk=kk+1; Irb(kk)=nb2; Jrb(kk)=nb2; Vrb(kk)=2*base;
    end
    Rb = sparse(Irb, Jrb, Vrb, NN, NN);

    % Top Dirichlet boundary
    D = (1:NGX).';
    U = (NGX+1:NN).';
    Pdir = ones(NGX,1);

    Kg_uu = Kg(U,U);
    Kg_ud = Kg(U,D);
    Rb_uu = Rb(U,U);
    Rb_ud = Rb(U,D);

    % surface extraction
    j0 = find(DZ==0, 1, 'first');
    if isempty(j0) || j0<=1 || j0>=numel(DZ)
        error('DZ==0 was not found or is invalid.');
    end
    kS = idxpl(:, j0);
    kU = idxpl(:, j0-1);
    kD = idxpl(:, j0+1);

    % default station coordinates
    % Release-clean version: do not load x_st_Sint_*.txt inside this function.
    % Station coordinates must be passed from the main inversion script using
    % the optional argument: 'x_st', x_st. If not provided, all surface nodes
    % are used as a fallback.
    if ~isempty(opt.x_st)
        xst_def = opt.x_st(:);
    else
        xst_def = DX(:);
    end

    ix_def = zeros(numel(xst_def),1);
    for ss = 1:numel(xst_def)
        [~, ix_def(ss)] = min(abs(DX - xst_def(ss)));
    end

    % ==== mapping CELL -> 4 TRIANGLES ====
    % earth layer dimulai pada box-row pertama setelah air
    % total box row earth layer = NDZb
    cell2tri = cell(NDX*NDZb,1);

    for izb = 1:NDZb
        j_box_global = NDZa + izb;   % global box row in the full mesh
        for ix = 1:NDX
            cell_id = (izb-1)*NDX + ix;
            box_id = (j_box_global-1)*NDX + ix;
            tri_ids = (box_id-1)*4 + (1:4);
            cell2tri{cell_id} = tri_ids(:);
        end
    end

    cache.NDX = NDX;
    cache.NDZ = NDZ;
    cache.NDZa = NDZa;
    cache.NDZb = NDZb;
    cache.NGX = NGX;
    cache.NGZ = NGZ;
    cache.NTE = NTE;
    cache.NN  = NN;
    cache.DX  = DX;
    cache.DZ  = DZ;
    cache.X   = X;
    cache.Z   = Z;
    cache.Zref = Zref;
    cache.idxpl = idxpl;
    cache.NE   = NE;
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
    cache.x_st_def = xst_def;
    cache.ix_st_def = ix_def;
end

% ===================== UNPACK CACHE =====================
NDX = cache.NDX;
NDZb = cache.NDZb;
NTE = cache.NTE;
NN  = cache.NN;
DX  = cache.DX;
Z   = cache.Z;
model_id = cache.model_id;
cell2tri = cache.cell2tri;

% ===================== MODEL CHECK =====================
Nmodel_cell = NDX*NDZb;
if numel(m_log10rho_cell) ~= Nmodel_cell
    error('m_cell length (%d) must equal NDX*NDZb (%d).', numel(m_log10rho_cell), Nmodel_cell);
end

% ===================== MAP CELL -> TRI =====================
m_tri_log10 = zeros(numel(model_id),1);

% peta global triangle id -> local earth index
tri_global_to_local = zeros(NTE,1);
tri_global_to_local(model_id) = 1:numel(model_id);

for ic = 1:Nmodel_cell
    tri_ids = cell2tri{ic};
    tri_loc = tri_global_to_local(tri_ids);
    m_tri_log10(tri_loc) = m_log10rho_cell(ic);
end

% ===================== SIGMA PER ELEM =====================
sigma = sigma_air * ones(NTE,1);
rho_earth_tri = 10.^m_tri_log10(:);
sigma(model_id) = 1 ./ rho_earth_tri;

% ===================== MASS MATRIX =====================
sigma_rep = repmat(sigma, 9, 1);
M = sparse(cache.I, cache.J, cache.mass_base .* sigma_rep, NN, NN);

U = cache.U; D = cache.D; Pdir = cache.Pdir;
M_uu = M(U,U);
M_ud = M(U,D);

Kg_uu = cache.Kg_uu; Kg_ud = cache.Kg_ud;
Rb_uu = cache.Rb_uu; Rb_ud = cache.Rb_ud;

kS = cache.kS; kU = cache.kU; kD = cache.kD;

% ===================== SOLVE PER FREQ =====================
Rhoapp  = zeros(NF, numel(DX));
FasaImp = zeros(NF, numel(DX));

for ifr = 1:NF
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

    Eu = Kuu \ rhs;

    E = complex(zeros(NN,1));
    E(U) = Eu;
    E(D) = Pdir;

    dz = Z(kD) - Z(kU);
    dE_dz = (E(kD) - E(kU)) ./ dz;

    Hx = (1i/omega_miu) .* dE_dz;
    Zte = E(kS) ./ Hx;

    Rhoapp(ifr,:)  = (abs(Zte).^2) ./ omega_miu;
    FasaImp(ifr,:) = angle(Zte)*180/pi;
end

% ===================== SAMPLE TO STATIONS =====================
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

rhoa_st    = Rhoapp(:, ix);
phi_st_deg = FasaImp(:, ix);

% ===================== OPTIONAL OUTPUT =====================
out = struct();
out.Rhoapp   = Rhoapp;
out.PhiDeg   = FasaImp;
out.DX       = DX;
out.x_st     = xst;
out.model_id = model_id;
out.NGX      = cache.NGX;
out.NTE      = cache.NTE;
out.NN       = cache.NN;
out.X        = cache.X;
out.Z        = cache.Z;
out.DZ       = cache.DZ;
out.NE       = cache.NE;
out.idxpl    = cache.idxpl;
out.Zref     = cache.Zref;
out.sigma    = sigma;
out.cell2tri = cell2tri;
end