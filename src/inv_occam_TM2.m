function R = inv_occam_TM2(NDX, nmodel, varargin)
%INV_OCCAM_TM2 First-order anisotropic Occam regularization matrix.
% R*m contains differences between neighboring model parameters.
% Horizontal differences are weighted by wx and vertical differences by wz.

p = inputParser;
p.addParameter('param_mode','auto',@(s)ischar(s)||isstring(s));
p.addParameter('wx',2.0,@isnumeric);
p.addParameter('wz',1.0,@isnumeric);
p.parse(varargin{:});

mode = lower(string(p.Results.param_mode));
wx   = p.Results.wx;
wz   = p.Results.wz;

if mode == "auto"
    if mod(nmodel,4*NDX)==0
        mode = "tri";
    elseif mod(nmodel,NDX)==0
        mode = "cell";
    else
        error('Cannot infer param_mode from NDX and nmodel.');
    end
end

switch mode
    case "cell"
        NDZ = round(nmodel / NDX);
        if NDX*NDZ ~= nmodel
            error('Inconsistent cell-model size.');
        end

        rows_i = [];
        cols_j = [];
        vals_v = [];
        row = 0;

        for iz = 1:NDZ
            for ix = 1:(NDX-1)
                id1 = (iz-1)*NDX + ix;
                id2 = id1 + 1;
                row = row + 1;
                rows_i(end+1:end+2) = [row row];
                cols_j(end+1:end+2) = [id1 id2];
                vals_v(end+1:end+2) = [-wx wx];
            end
        end

        for iz = 1:(NDZ-1)
            for ix = 1:NDX
                id1 = (iz-1)*NDX + ix;
                id2 = id1 + NDX;
                row = row + 1;
                rows_i(end+1:end+2) = [row row];
                cols_j(end+1:end+2) = [id1 id2];
                vals_v(end+1:end+2) = [-wz wz];
            end
        end

        R = sparse(rows_i, cols_j, vals_v, row, nmodel);

    case "tri"
        rows_i = [];
        cols_j = [];
        vals_v = [];
        row = 0;

        for k = 1:(nmodel-1)
            row = row + 1;
            rows_i(end+1:end+2) = [row row];
            cols_j(end+1:end+2) = [k k+1];
            vals_v(end+1:end+2) = [-1 1];
        end

        R = sparse(rows_i, cols_j, vals_v, row, nmodel);

    otherwise
        error('param_mode must be cell, tri, or auto.');
end
end
