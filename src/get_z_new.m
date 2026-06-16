function Z_all = get_z_new(z_1,z_end,nlayer)

t = z_1;
Zl = z_end;
l = nlayer - 1;
b0 = 5;

for it = 1:1e2
    Zcal = t*(1-b0^l)/(1-b0);
    J = (Zcal-l*t*b0^(l-1))/(1-b0);

    deld = Zl - Zcal;

    if deld^2 < 1e-7
        break
    end

    delb = deld/J;
    b0 = b0 + delb;
end

Z_all = [0 cumsum(t*b0.^((1:l)-1))]';

