
%---------------------------------------------------------------%
%             CSVWRITER - .DAT Factory for PGF                  %
%---------------------------------------------------------------%

Z = hn;   % escolho aqui a matriz referente a Z

nc = size(Z,1);

Zc = Z(:,1);

for i = 2:nc
    Zc = [Zc ; Z(:,i)];
end

X1 = ones(nc,1);

for i = 2:nc
    X1 = [X1 ; X1(1:nc)*i];
end

X2 = linspace(1,nc,nc)';

for i = 2:nc;
    X2 = [X2 ; X2(1:nc)];
end

ZV = [X1 X2 Zc];

save ztemp.dat ZV -ASCII