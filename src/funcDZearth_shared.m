function DZb = funcDZearth_shared(NDZb)
%FUNCDZEARTH_SHARED Build the shared earth-depth node grid for TE and TM.
% The output starts at zero and increases monotonically with depth.

NDrp = 20;
NDrg = NDZb - NDrp;
tebal_rp = 25;

BRp  = 5000;
ZMAX = 550000;
BRg  = ZMAX - BRp;

Zrp = get_z_new(tebal_rp, BRp, NDrp);
tebal_rg = Zrp(end) - Zrp(end-1);

Zrg = get_z_new(tebal_rg, BRg, NDrg+2) + BRp;

DZb = [Zrp; Zrg(2:end)];
DZb = DZb(:);

if abs(DZb(1)) > 1e-12
    error('funcDZearth_shared must return zero as the first node.');
end
if any(diff(DZb) <= 0)
    error('funcDZearth_shared must return monotonically increasing nodes.');
end
end
