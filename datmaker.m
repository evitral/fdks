
%---------------------------------------------------------------%
%             CSVWRITER - .DAT Factory for PGF                  %
%---------------------------------------------------------------%

% fid = fopen('L2np129.dat','w');
% for t = 1:400
% fprintf(fid, '%6.0f %12.4e \n',t,L2np129(t));
% end
% fclose(fid);

Z = hy;   % escolho aqui a matriz referente a Z

dX = 1;

nc = size(Z,1);

Zc = Z(:,1);

for i = 2:nc
    Zc = [Zc ; Z(:,i)];
end

X1 = zeros(nc,1);
Xones = ones(nc,1);

for i = 1:(nc-1)
    X1 = [X1 ; Xones(1:nc)*dX*i];
end

X2 = linspace(0,(dX*nc-dX),nc)';

for i = 2:nc;
    X2 = [X2 ; X2(1:nc)];
end

%X1 = X1 - 255;
%X2 = X2 - 255;

X1 = X1 - 128;
X2 = X2 - 128;

ZV = [X1 X2 Zc];

save h41960fft.dat ZV -ASCII