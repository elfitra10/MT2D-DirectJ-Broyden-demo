function NE = konektivitas(NTE, NGX, NDX, NDZ)
%KONEKTIVITAS Build triangular-element connectivity for a structured grid.
% Each rectangular cell is split into four triangular elements.

k = 0;
l = 0;
NE = zeros(3,NTE);
NEXZ = zeros(5,1);

for j = 1:NDZ
    for i = 1:NDX
        k = k+1;
        l = l+1;

        NEXZ(1) = k;
        NEXZ(2) = k+1;
        NEXZ(3) = k+NGX;
        NEXZ(4) = k+NGX+NDX;
        NEXZ(5) = NEXZ(4)+1;

        m = (l-1)*4;
        NE(:,m+1) = [NEXZ(1); NEXZ(2); NEXZ(3)];
        NE(:,m+2) = [NEXZ(2); NEXZ(5); NEXZ(3)];
        NE(:,m+3) = [NEXZ(5); NEXZ(4); NEXZ(3)];
        NE(:,m+4) = [NEXZ(4); NEXZ(1); NEXZ(3)];
    end
    k = k+NGX;
end
end
