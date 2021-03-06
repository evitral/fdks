
%---------------------------------------------------------------%
% Finite Differences Numerical Scheme for Solving a Generalized %
% Kuramoto-Sivashinsky Equation Modelling Nano-Patterning of    %
% Surfaces by Ion Sputtering                                    %
%---------------------------------------------------------------%

clear all
clear variables

fid = fopen('log10L1.txt','w');
fprintf(fid, ['log10L1' '\n']);

condName = 'PRdt1d0np128dx1d0rand.mat';
condNameA = strcat('a',condName);
condNameB = strcat('b',condName);

disp('********************************')
disp(['   FDKS Simulation'])
disp(['   ',condName])
disp(['   Start: ',datestr(now)])

ld = 2;   % ld = 1: load file

if ld == 1
    
    load(condName)
    t1 = t;

else    

    t1 = 1;
    
%---------------------------------------------------------------%
%                Parameters of the dynamics                     %
%---------------------------------------------------------------%

alpha = 0.15;      % unclear physical meaning (B&H, alpha = 0)
teta  = 0.5236;    % angle of incidence of the beam (normal, teta=0)

s = sin(teta);     % sin teta
c = cos(teta);     % cos teta
a = 2;             % penetration depth (nm)
mu = 0.5;          % width of the energy distribution (nm)
amu = a/mu;
J = 10;            % flux of bombarding ions (nm^-2 s^-1)
ei = 500;          % energy carried by the ions (eV)
p = 0.3;           % depends of surface binding energy and scattering
                   % cross section (nm^4/eV)
                   
F = (J*ei*p/sqrt(2*pi()*mu))*exp(-amu^2*c^2/2); % VERIFICAR UNIDADES

T = 500;           % temperature (K)
dE = 1.25;         % surface diffusion activate energy

%dt = 0.01;        % time step (s) ???

%---------------------------------------------------------------%
%                           Mesh                                %
%---------------------------------------------------------------%

L = 128;

np = 128;

dX = L/np;
dY = dX;

%---------------------------------------------------------------%
%                    Transformed parameters                     %
%---------------------------------------------------------------%

%dtau = 2*F*amu^2*dt/a;
dtau = 1.0;
ntau = 7000;
tmax = ntau*dtau;

Kt = (amu^2*0.54*10^16/T)*exp(-dE*1.16*10^4/T)*exp(-amu^2*s^2/2);
Kt = 5.0;
mut = 2*s^2-c^2-(amu*s*c)^2;

nuxt = c*(3*s^2-c^2-(amu*s*c)^2); 

%---------------------------------------------------------------%
%                      Useful constants                         %
%---------------------------------------------------------------%

Dxx = (c^2-4*s^2+2*amu^2*s^2*(c^2-(2/3)*s^2)+(amu^4/3)*s^4*c^2);

Dxy = 2*(c^2-2*s^2+(amu*s*c)^2); % 2 times Dxy actually

alphat = (dtau/4)*alpha;

beta = (Dxx+Kt)/(dX^4);
beta2 = (dtau/2)*beta;  
gamma = Kt/(dY^4);
gamma2 = (dtau/2)*gamma;

delta = mut/dX^2;   % d^2h/dX^2
zeta = c^2/dY^2;   % d^2h/dY^2

eta = nuxt/(2*dX^2);   % (dh/dX)^2
psi = c^3/(2*dY^2);   % (dh/dY)^2

omega = (Dxy-2*Kt)/(dX^2*dY^2);   % d^4h/dX^2dY^2
phi = c^2/dY^4;   % d^4h/dY^4

deck = [delta,zeta,eta,psi,omega,phi];

sg2at = 6*gamma2+alphat;

npm = np-1;
npmm = np-2;
npmmm = np-3;

%---------------------------------------------------------------%
%                           Operators                           %
%---------------------------------------------------------------%

E = zeros(np,np);     % identity matrix
for i = 1:np
    E(i,i) = 1;
end

Lx = zeros(np,np);    % operator delta tau * lambda X
Lx(1:2,1:4) = beta2*[-6 4 -1 0;4 -6 4 -1]+alphat*[-1 0 0 0;0 -1 0 0];
Lx(1:2,(np-1):np) = beta2*[-1 4;0 -1];
Lx((np-1):np,(np-3):np) = beta2*[-1 4 -6 4;0 -1 4 -6]+...
    alphat*[0 0 -1 0;0 0 0 -1];
Lx((np-1):np,1:2) = beta2*[-1 0;4 -1];
for i = 3:(np-2)
    Lx(i,(i-2):(i+2)) = beta2*[-1 4 -6 4 -1]+alphat*[0 0 -1 0 0];
end

Ly = zeros(np,np);    % operator delta tau * lambda Y
Ly(1:2,1:4) = gamma2*[-6 4 -1 0;4 -6 4 -1] + alphat*[-1 0 0 0;0 -1 0 0];
Ly(1:2,(np-1):np) = gamma2*[-1 4;0 -1];
Ly((np-1):np,(np-3):np) = gamma2*[-1 4 -6 4;0 -1 4 -6]+...
    alphat*[0 0 -1 0;0 0 0 -1];
Ly((np-1):np,1:2) = gamma2*[-1 0;4 -1];

for i = 3:npmm
    Ly(i,(i-2):(i+2)) = gamma2*[-1 4 -6 4 -1]+alphat*[0 0 -1 0 0];
end

%---------------------------------------------------------------%
% Here we define the matrices appearing on the following eqs:   %
% (E-dtaoLx)htil = (E+dtaoLy)hn + 2dtaoLyhn + dtaofn+1/2        %
% (E-dtaoLy)hn,m+1 = htil - dtaoLyhn                            %
%---------------------------------------------------------------%

M1 = E - Lx;
M1i = inv(M1);

M2 = E - Ly;
M2i = inv(M2);

M3 = E + Lx;

M4 = 2*Ly;     % Dinossaur

M5 = Ly;

%---------------------------------------------------------------%
%                    Initial conditions                         %
%---------------------------------------------------------------%

hn = zeros(np,np);
hm = zeros(np,np);
htil = zeros(np,np);
hm1 = zeros(np,np);

for y = 1:np
    for x = 1:np
        %hn(y,x) = 0.1*(2*rand-1);
        hn(y,x)= 0.1*rand;
        %hn(y,x) = 0.1*sin(8*2*pi()*x/(np+1));
    end
end

h0 = hn;       % keeps the initial conditions

%---------------------------------------------------------------%
%                        Finders Keepers                        %
%---------------------------------------------------------------%

tdata = 10.;    % data stored after \tdata time steps
tsave = 500.;   % data saved after \tdata time steps

maxHN = zeros(1,ntau/tdata); % Keeps highest hn values
maxIT = zeros(1,ntau/tdata); % Records the internal iterations value

movieA = [10  30  60 100 150 200 250 300 350 400 ...
          450 500 600 700 850 1000 1250 1500 2000 2500 ... 
          3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 ...
          8000 8500 9000 9500 10000 10500 11000 11500 12000 12500 ...
          13000 14000 15000 16000 17000 18000 20000 23000 26000 30000 ...
          50000 100000 150000 200000 250000 300000 350000 400000 450000 ...
          500000 550000 600000 650000 700000 800000 900000];     
 
movieH = zeros(np,np,66);

% Probing %

pb1 = 43;   % np = 128
pb2 = 86;   % np = 128

L1p1  = zeros(1,ntau/tdata);
L1p2  = zeros(1,ntau/tdata);
L1p3  = zeros(1,ntau/tdata);
L1p4  = zeros(1,ntau/tdata);
L1p5  = zeros(1,ntau/tdata);
L1p6  = zeros(1,ntau/tdata);
L1p7  = zeros(1,ntau/tdata);
L1p8  = zeros(1,ntau/tdata);
L1p9  = zeros(1,ntau/tdata);
L1p10 = zeros(1,ntau/tdata);
L1p11 = zeros(1,ntau/tdata);
L1p12 = zeros(1,ntau/tdata);
L1p13 = zeros(1,ntau/tdata);
L1p14 = zeros(1,ntau/tdata);
L1p15 = zeros(1,ntau/tdata);
L1p16 = zeros(1,ntau/tdata);

%---------------------------------------------------------------%
%                             Main                              %
%---------------------------------------------------------------%

Mfn = zeros(np,np);
Mfm = zeros(np,np);
Mf  = zeros(np,np);

Lyh=zeros(np,np);

nmax = int32(tmax/dtau+1);
L1=zeros(1,nmax/tdata);

flag = 1; % L1 control

end

saveMode = 1;

sg2at = 6*gamma2+alphat;

% Vectorization %

vlin = 1:np;
vcol = 1:np;

vlinmm = 3:npmm;
vcolmm = 3:npmm;

for t = t1:nmax
    
    pass = 1;
    k = 1;
    hm = hn;
    
%---------------------------------------------------------------%
%                       Assembly Mfn                            %
%---------------------------------------------------------------%

    % Line 1 %
        
    Mfn(1,1) = fn(deck,hn,[1,2,3,np,npm,1,2,3,np,npm]);

    Mfn(1,2) = fn(deck,hn,[1,2,3,np,npm,2,3,4,1,np]);

    Mfn(1,npm) = fn(deck,hn,[1,2,3,np,npm,npm,np,1,npmm,npmmm]);

    Mfn(1,np) = fn(deck,hn,[1,2,3,np,npm,np,1,2,npm,npmm]);

    % Line 2 %

    Mfn(2,1) = fn(deck,hn,[2,3,4,1,np,1,2,3,np,npm]);

    Mfn(2,2) = fn(deck,hn,[2,3,4,1,np,2,3,4,1,np]);

    Mfn(2,npm) = fn(deck,hn,[2,3,4,1,np,npm,np,1,npmm,npmmm]);

    Mfn(2,np) = fn(deck,hn,[2,3,4,1,np,np,1,2,npm,npmm]);

	% Superior and inferior borders %

	for x = 3:npmm

    Mfn(1,x) = fn(deck,hn,[1,2,3,np,npm,x,x+1,x+2,x-1,x-2]);
	
	Mfn(2,x) = fn(deck,hn,[2,3,4,1,np,x,x+1,x+2,x-1,x-2]);

    Mfn(npm,x) = fn(deck,hn,[npm,np,1,npmm,npmmm,x,x+1,x+2,x-1,x-2]);
    
    Mfn(np,x) = fn(deck,hn,[np,1,2,npm,npmm,x,x+1,x+2,x-1,x-2]);

	end

	% Middle %
   
    for y = 3:npmm
    
    yp=y+1; ypp=y+2; ym=y-1; ymm=y-2;    
        
    Mfn(y,1) = fn(deck,hn,[y,yp,ypp,ym,ymm,1,2,3,np,npm]);
    
    Mfn(y,2) = fn(deck,hn,[y,yp,ypp,ym,ymm,2,3,4,1,np]);
    
        for i = 3:npmm
        Mfn(y,i) = fn(deck,hn,[y,yp,ypp,ym,ymm,i,i+1,i+2,i-1,i-2]);
        end
    
    Mfn(y,npm) = fn(deck,hn,[y,yp,ypp,ym,ymm,npm,np,1,npmm,npmmm]);
    
    Mfn(y,np) = fn(deck,hn,[y,yp,ypp,ym,ymm,np,1,2,npm,npmm]);
    
    end

    % Line np-1 %

    Mfn(npm,1) = fn(deck,hn,[npm,np,1,npmm,npmmm,1,2,3,np,npm]);

    Mfn(npm,2) = fn(deck,hn,[npm,np,1,npmm,npmmm,2,3,4,1,np]);

    Mfn(npm,npm) = fn(deck,hn,[npm,np,1,npmm,npmmm,npm,np,1,npmm,npmmm]);

    Mfn(npm,np) = fn(deck,hn,[npm,np,1,npmm,npmmm,np,1,2,npm,npmm]);

    % Line np %

    Mfn(np,1) = fn(deck,hn,[np,1,2,npm,npmm,1,2,3,np,npm]);
    
    Mfn(np,2) = fn(deck,hn,[np,1,2,npm,npmm,2,3,4,1,np]);
    
    Mfn(np,npm) = fn(deck,hn,[np,1,2,npm,npmm,npm,np,1,npmm,npmmm]);
    
    Mfn(np,np) = fn(deck,hn,[np,1,2,npm,npmm,np,1,2,npm,npmm]);

%---------------------------------------------------------------%
%                    Operator Y - rows                          %
%---------------------------------------------------------------%

    Lyh(1,vcol) = -gamma2*hn(npm,vcol)+4*gamma2*hn(np,vcol)...
       -sg2at*hn(1,vcol)+4*gamma2*hn(2,vcol)-gamma2*hn(3,vcol);
   
    Lyh(2,vcol) = -gamma2*hn(np,vcol)+4*gamma2*hn(1,vcol)...
       -sg2at*hn(2,vcol)+4*gamma2*hn(3,vcol)-gamma2*hn(4,vcol);
   
    Lyh(npm,vcol) = -gamma2*hn(npmmm,vcol)+4*gamma2*hn(npmm,vcol)...
       -sg2at*hn(npm,vcol)+4*gamma2*hn(np,vcol)-gamma2*hn(1,vcol);
   
    Lyh(np,vcol) = -gamma2*hn(npmm,vcol)+4*gamma2*hn(npm,vcol)...
       -sg2at*hn(np,vcol)+4*gamma2*hn(1,vcol)-gamma2*hn(2,vcol);
 
    Lyh(vlinmm,vcol) = -gamma2*hn(vlinmm-2,vcol)+ ...
        4*gamma2*hn(vlinmm-1,vcol)-sg2at*hn(vlinmm,vcol)+ ...
        4*gamma2*hn(vlinmm+1,vcol)-gamma2*hn(vlinmm+2,vcol);    
     
%---------------------------------------------------------------%
    
    while pass == 1
        
%---------------------------------------------------------------%
%                       Assembly Mfm                            %
%---------------------------------------------------------------%

    % Line 1 %
        
    Mfm(1,1) = fm(deck,hm,hn,[1,2,3,np,npm,1,2,3,np,npm]);

    Mfm(1,2) = fm(deck,hm,hn,[1,2,3,np,npm,2,3,4,1,np]);

    Mfm(1,npm) = fm(deck,hm,hn,[1,2,3,np,npm,npm,np,1,npmm,npmmm]);

    Mfm(1,np) = fm(deck,hm,hn,[1,2,3,np,npm,np,1,2,npm,npmm]);

    % Line 2 %

    Mfm(2,1) = fm(deck,hm,hn,[2,3,4,1,np,1,2,3,np,npm]);

    Mfm(2,2) = fm(deck,hm,hn,[2,3,4,1,np,2,3,4,1,np]);

    Mfm(2,npm) = fm(deck,hm,hn,[2,3,4,1,np,npm,np,1,npmm,npmmm]);

    Mfm(2,np) = fm(deck,hm,hn,[2,3,4,1,np,np,1,2,npm,npmm]);

	% Superior and inferior borders %

	for x = 3:npmm

    Mfm(1,x) = fm(deck,hm,hn,[1,2,3,np,npm,x,x+1,x+2,x-1,x-2]);
	
	Mfm(2,x) = fm(deck,hm,hn,[2,3,4,1,np,x,x+1,x+2,x-1,x-2]);

    Mfm(npm,x) = fm(deck,hm,hn,[npm,np,1,npmm,npmmm,x,x+1,x+2,x-1,x-2]);
    
    Mfm(np,x) = fm(deck,hm,hn,[np,1,2,npm,npmm,x,x+1,x+2,x-1,x-2]);

	end

	% Middle %

    for y = 3:npmm
    
    yp=y+1; ypp=y+2; ym=y-1; ymm=y-2; 
        
    Mfm(y,1) = fm(deck,hm,hn,[y,yp,ypp,ym,ymm,1,2,3,np,npm]);
    
    Mfm(y,2) = fm(deck,hm,hn,[y,yp,ypp,ym,ymm,2,3,4,1,np]);

        for x = 3:npmm            
        Mfm(y,x) = fm(deck,hm,hn,[y,yp,ypp,ym,ymm,x,x+1,x+2,x-1,x-2]);
        end

    Mfm(y,npm) = fm(deck,hm,hn,[y,yp,ypp,ym,ymm,npm,np,1,npmm,npmmm]);
    
    Mfm(y,np) = fm(deck,hm,hn,[y,yp,ypp,ym,ymm,np,1,2,npm,npmm]);
    
    end

    % Line np-1 %

    Mfm(npm,1) = fm(deck,hm,hn,[npm,np,1,npmm,npmmm,1,2,3,np,npm]);

    Mfm(npm,2) = fm(deck,hm,hn,[npm,np,1,npmm,npmmm,2,3,4,1,np]);

    Mfm(npm,npm) = fm(deck,hm,hn,[npm,np,1,npmm,npmmm,npm,np,1,npmm,npmmm]);

    Mfm(npm,np) = fm(deck,hm,hn,[npm,np,1,npmm,npmmm,np,1,2,npm,npmm]);

    % Line np %

    Mfm(np,1) = fm(deck,hm,hn,[np,1,2,npm,npmm,1,2,3,np,npm]);
    
    Mfm(np,2) = fm(deck,hm,hn,[np,1,2,npm,npmm,2,3,4,1,np]);
    
    Mfm(np,npm) = fm(deck,hm,hn,[np,1,2,npm,npmm,npm,np,1,npmm,npmmm]);
    
    Mfm(np,np) = fm(deck,hm,hn,[np,1,2,npm,npmm,np,1,2,npm,npmm]);
      
%---------------------------------------------------------------%
%                      Splitting Scheme                         %
%---------------------------------------------------------------%

       Mf = dtau*(Mfn + Mfm);
       
        for lin = 1:np
             htil(lin,:) = M1i*(M3*hn(lin,:)'+2*Lyh(lin,:)' + Mf(lin,:)');
        end
       
       hm1(:,vcol) = M2i*(htil(:,vcol)-M5*hn(:,vcol));
       
        if (max(max(abs(hm1-hm)))/max(max(abs(hm1)))) < 1e-8
            pass = 2;
        elseif k == 50;
            disp('   Iterations number exceeded!')
            flag = 2;
            break;
        else
            hm = hm1;
            k = k+1;   
        end
    
    end
    
    
%---------------------------------------------------------------%
%              Routine after each 10 time steps                 %
%---------------------------------------------------------------%
    
    if mod(t,tdata) == 0
    
    tr = t/tdata;
        
    %ghn(:,:,(t/10))=hn;
    maxHN(tr) = max(max(abs(hn)));  % Keeps max values
    
    maxIT(tr) = k;
    
    if ismember(t,movieA) == 1
        movieH(:,:,find(movieA==t)) = hm1;
    end
    
    % Probing %
    
    L1p1(tr)  = abs(hm1(1,1)-hn(1,1))/abs(hm1(1,1));
    L1p2(tr)  = abs(hm1(1,pb1)-hn(1,pb1))/abs(hm1(1,pb1));
    L1p3(tr)  = abs(hm1(1,pb2)-hn(1,pb2))/abs(hm1(1,pb2));
    L1p4(tr)  = abs(hm1(1,np)-hn(1,np))/abs(hm1(1,np));
    L1p5(tr)  = abs(hm1(pb1,1)-hn(pb1,1))/abs(hm1(pb1,1));
    L1p6(tr)  = abs(hm1(pb1,pb1)-hn(pb1,pb1))/abs(hm1(pb1,pb1));
    L1p7(tr)  = abs(hm1(pb1,pb2)-hn(pb1,pb2))/abs(hm1(pb1,pb2));
    L1p8(tr)  = abs(hm1(pb1,np)-hn(pb1,np))/abs(hm1(pb1,np));    
    L1p9(tr)  = abs(hm1(pb2,1)-hn(pb2,1))/abs(hm1(pb2,1));
    L1p10(tr) = abs(hm1(pb2,pb1)-hn(pb2,pb1))/abs(hm1(pb2,pb1));
    L1p11(tr) = abs(hm1(pb2,pb2)-hn(pb2,pb2))/abs(hm1(pb2,pb2));
    L1p12(tr) = abs(hm1(pb2,np)-hn(pb2,np))/abs(hm1(pb2,np));    
    L1p13(tr) = abs(hm1(np,1)-hn(np,1))/abs(hm1(np,1));
    L1p14(tr) = abs(hm1(np,pb1)-hn(np,pb1))/abs(hm1(np,pb1));
    L1p15(tr) = abs(hm1(np,pb2)-hn(np,pb2))/abs(hm1(np,pb2));
    L1p16(tr) = abs(hm1(np,np)-hn(np,np))/abs(hm1(np,np));    
    
    % Calculating L1 %
    
    Sup=0;
    Sin=0;
    for j = 1:np
        for i = i:np
            Sup = Sup + abs(hm1(j,i)-hn(j,i));
            Sin = Sin + abs(hm1(j,i));
        end
    end
    
    L1(tr) = Sup/(dtau*Sin);
    fprintf(fid,'%f \n',log10(L1(tr)));
    
    % L1 as a LimitBreak %
    
    if L1(tr) <1E-20  % OLHAAA AQUI OHHH
        flag = 2;
    end
    
    if mod(t,tsave) == 0
        hn = hm1;
        if saveMode == 1
         	save(condNameA)
            saveMode = 2;
        else
            save(condNameB)
            saveMode = 1;
        end
    end
        
    end
    
    if flag == 2;
        break;
    end
    
    % New matrices are old matrices %
    
    hn = hm1;

end

save(condName)

disp(['   End:   ',datestr(now)])
disp('********************************')