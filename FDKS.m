
%---------------------------------------------------------------%
% Finite Differences Numerical Scheme for Solving a Generalized %
% Kuramoto-Sivashinsky Equation Modelling Nano-Patterning of    %
% Surfaces by Ion Sputtering                                    %
%---------------------------------------------------------------%

clear all
clear variables

fid = fopen('log10L1.txt','w');
fprintf(fid, ['log10L1' '\n']);

ld = 2;   % ld = 1: load file

if ld == 1
    
    load arquivo.mat
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

L = 510;

np = 256;

dX = L/(np-1);
dY = dX;

%---------------------------------------------------------------%
%                    Transformed parameters                     %
%---------------------------------------------------------------%

%dtau = 2*F*amu^2*dt/a;
dtau = 0.1;
ntau = 1000000;
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

for i = 3:(np-2)
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
%ghn = zeros(np,np,ntau/10);
hm = zeros(np,np);
htil = zeros(np,np);
hm1 = zeros(np,np);

for j = 1:np
    for i = 1:np
        %hn(j,i) = 0.1*(2*rand-1);
        %hn(j,i)= 0.1*rand;
        hn(j,i) = 0.1*sin(14*2*pi()*j/(np+1));
    end
end

h0 = hn;       % keeps the initial conditions

%---------------------------------------------------------------%
%                        Finders Keepers                        %
%---------------------------------------------------------------%

maxHN = zeros(1,ntau/10); % Keeps highest hn values
maxIT = zeros(1,ntau/10); % Records the internal iterations value

movieA = [10  30  60 100 150 200 250 300 350 400 ...
          450 500 600 700 850 1000 1250 1500 2000 2500 ... 
          3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 ...
          8000 8500 9000 9500 10000 10500 11000 11500 12000 12500 ...
          13000 14000 15000 16000 17000 18000 20000 23000 26000 30000 ...
          50000 100000 150000 200000 250000 300000 350000 400000 450000 ...
          500000 550000 600000 650000 700000 800000 900000];     
 
movieH = zeros(np,np,66);

%---------------------------------------------------------------%
%                             Main                              %
%---------------------------------------------------------------%

Mfn = zeros(np,np);
Mfm = zeros(np,np);
Mf  = zeros(np,np);

Lyh=zeros(np,np);

nmax = int32(tmax/dtau+1);
L1=zeros(1,nmax/10);

flag = 1; % L1 control

end

for t = t1:nmax
    pass = 1;
    k = 1;
    hm = hn;
    
%---------------------------------------------------------------%
%                       Assembly Mfn                            %
%---------------------------------------------------------------%

    % Line 1 %
    
    Vji = [1,2,3,np,npm,1,2,3,np,npm];
    Mfn(1,1) = fn(deck,hn,Vji);

    Vji(6:10) = [2,3,4,1,np];
    Mfn(1,2) = fn(deck,hn,Vji);

    for i = 3:npmm
        Vji(6:10) = [i,i+1,i+2,i-1,i-2];
        Mfn(1,i) = fn(deck,hn,Vji);
    end

    Vji(6:10) = [npm,np,1,npmm,npmmm];
    Mfn(1,npm) = fn(deck,hn,Vji);

    Vji(6:10) = [np,1,2,npm,npmm];
    Mfn(1,np) = fn(deck,hn,Vji);

    % Line 2 %

    Vji = [2,3,4,1,np,1,2,3,np,npm];
    Mfn(2,1) = fn(deck,hn,Vji);

    Vji(6:10) = [2,3,4,1,np];
    Mfn(2,2) = fn(deck,hn,Vji);

    for i = 3:npmm
        Vji(6:10) = [i,i+1,i+2,i-1,i-2];
        Mfn(2,i) = fn(deck,hn,Vji); 
    end

    Vji(6:10) = [npm,np,1,npmm,npmmm];
    Mfn(2,npm) = fn(deck,hn,Vji);

    Vji(6:10) = [np,1,2,npm,npmm];
    Mfn(2,np) = fn(deck,hn,Vji);

    % Main %

    for j = 3:npmm
        Vji = [j,j+1,j+2,j-1,j-2,1,2,3,np,npm];
        Mfn(j,1) = fn(deck,hn,Vji);

        Vji(6:10) = [2,3,4,1,np];
        Mfn(j,2) = fn(deck,hn,Vji);

        for i = 3:npmm
            Vji(6:10) = [i,i+1,i+2,i-1,i-2];
            Mfn(j,i) = fn(deck,hn,Vji);
        end
    
        Vji(6:10) = [npm,np,1,npmm,npmmm];
        Mfn(j,(np-1)) = fn(deck,hn,Vji);

        Vji(6:10) = [np,1,2,npm,npmm];
        Mfn(j,np) = fn(deck,hn,Vji);
    end

    % Line np-1 %

    Vji = [npm,np,1,npmm,npmmm,1,2,3,np,npm];
    Mfn(npm,1) = fn(deck,hn,Vji);

    Vji(6:10) = [2,3,4,1,np];
    Mfn(npm,2) = fn(deck,hn,Vji);

    for i = 3:npmm
        Vji(6:10) = [i,i+1,i+2,i-1,i-2];
        Mfn(npm,i) = fn(deck,hn,Vji);
    end

    Vji(6:10) = [npm,np,1,npmm,npmmm];
    Mfn(npm,npm) = fn(deck,hn,Vji);

    Vji(6:10) = [np,1,2,npm,npmm];
    Mfn(npm,np) = fn(deck,hn,Vji);

    % Line np %

    Vji = [np,1,2,npm,npmm,1,2,3,np,npm];
    Mfn(np,1) = fn(deck,hn,Vji);

    Vji(6:10) = [2,3,4,1,np];
    Mfn(np,2) = fn(deck,hn,Vji);

    for i = 3:npmm
        Vji(6:10) = [i,i+1,i+2,i-1,i-2];
        Mfn(np,i) = fn(deck,hn,Vji);
    end

    Vji(6:10) = [npm,np,1,npmm,npmmm];
    Mfn(np,npm) = fn(deck,hn,Vji);

    Vji(6:10) = [np,1,2,npm,npmm];
    Mfn(np,np) = fn(deck,hn,Vji);

%---------------------------------------------------------------%
%         Beta Test                                             %
%---------------------------------------------------------------%
        
    for col = 1:np
       Lyh(1,col) = -gamma2*hn(npm,col)+4*gamma2*hn(np,col)...
       -(6*gamma2+alphat)*hn(1,col)+4*gamma2*hn(2,col)-gamma2*hn(3,col);
   
       Lyh(2,col) = -gamma2*hn(np,col)+4*gamma2*hn(1,col)...
       -(6*gamma2+alphat)*hn(2,col)+4*gamma2*hn(3,col)-gamma2*hn(4,col);
   
       Lyh(npm,col) = -gamma2*hn(npmmm,col)+4*gamma2*hn(npmm,col)...
       -(6*gamma2+alphat)*hn(npm,col)+4*gamma2*hn(np,col)-gamma2*hn(1,col);
   
       Lyh(np,col) = -gamma2*hn(npmm,col)+4*gamma2*hn(npm,col)...
       -(6*gamma2+alphat)*hn(np,col)+4*gamma2*hn(1,col)-gamma2*hn(2,col);
    end
    
    for lin = 3:npmm
        for col = 1:np
        Lyh(lin,col) = -gamma2*hn(lin-2,col)+4*gamma2*hn(lin-1,col)...
        -(6*gamma2+alphat)*hn(lin,col)+4*gamma2*hn(lin+1,col)...
        -gamma2*hn(lin+2,col);
        end        
    end
    
%---------------------------------------------------------------%
    
    while pass == 1
        
%---------------------------------------------------------------%
%                       Assembly Mfm                            %
%---------------------------------------------------------------%

        % Line 1 %
        
        Vji = [1,2,3,np,npm,1,2,3,np,npm];
        Mfm(1,1) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [2,3,4,1,np];
        Mfm(1,2) = fm(deck,hm,hn,Vji);

        for i = 3:npmm
            Vji(6:10) = [i,i+1,i+2,i-1,i-2];
            Mfm(1,i) = fm(deck,hm,hn,Vji);
        end

        Vji(6:10) = [npm,np,1,npmm,npmmm];
        Mfm(1,npm) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [np,1,2,npm,npmm];
        Mfm(1,np) = fm(deck,hm,hn,Vji);

        % Line 2 %

        Vji = [2,3,4,1,np,1,2,3,np,npm];
        Mfm(2,1) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [2,3,4,1,np];
        Mfm(2,2) = fm(deck,hm,hn,Vji);

        for i = 3:npmm
            Vji(6:10) = [i,i+1,i+2,i-1,i-2];
            Mfm(2,i) = fm(deck,hm,hn,Vji);
        end

        Vji(6:10) = [npm,np,1,npmm,npmmm];
        Mfm(2,npm) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [np,1,2,npm,npmm];
        Mfm(2,np) = fm(deck,hm,hn,Vji);

        % Main %

        for j = 3:npmm
            Vji = [j,j+1,j+2,j-1,j-2,1,2,3,np,npm];
            Mfm(j,1) = fm(deck,hm,hn,Vji);

            Vji(6:10) = [2,3,4,1,np];
            Mfm(j,2) = fm(deck,hm,hn,Vji);

            for i = 3:npmm
                Vji(6:10) = [i,i+1,i+2,i-1,i-2];
                Mfm(j,i) = fm(deck,hm,hn,Vji);
            end
    
            Vji(6:10) = [npm,np,1,npmm,npmmm];
            Mfm(j,(np-1)) = fm(deck,hm,hn,Vji);

            Vji(6:10) = [np,1,2,npm,npmm];
            Mfm(j,np) = fm(deck,hm,hn,Vji);
    
        end

        % Line np-1 %

        Vji = [npm,np,1,npmm,npmmm,1,2,3,np,npm];
        Mfm(npm,1) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [2,3,4,1,np];
        Mfm(npm,2) = fm(deck,hm,hn,Vji);

        for i = 3:npmm
            Vji(6:10) = [i,i+1,i+2,i-1,i-2];
            Mfm(npm,i) = fm(deck,hm,hn,Vji);
        end

        Vji(6:10) = [npm,np,1,npmm,npmmm];
        Mfm(npm,npm) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [np,1,2,npm,npmm];
        Mfm(npm,np) = fm(deck,hm,hn,Vji);

        % Line np %

        Vji = [np,1,2,npm,npmm,1,2,3,np,npm];
        Mfm(np,1) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [2,3,4,1,np];
        Mfm(np,2) = fm(deck,hm,hn,Vji);

        for i = 3:npmm
            Vji(6:10) = [i,i+1,i+2,i-1,i-2];
            Mfm(np,i) = fm(deck,hm,hn,Vji);
        end

        Vji(6:10) = [npm,np,1,npmm,npmmm];
        Mfm(np,npm) = fm(deck,hm,hn,Vji);

        Vji(6:10) = [np,1,2,npm,npmm];
        Mfm(np,np) = fm(deck,hm,hn,Vji);
      
%---------------------------------------------------------------%
%                      Splitting Scheme                         %
%---------------------------------------------------------------%

       Mf = dtau*(Mfn + Mfm);
       
        for lin = 1:np
            htil(lin,:) = M1i*(M3*hn(lin,:)'+2*Lyh(lin,:)' + Mf(lin,:)');
        end
        
        for col = 1:np
            hm1(:,col) = M2i*(htil(:,col)-M5*hn(:,col));
        end
        
        if (max(max(abs(hm1-hm)))/max(max(abs(hm1)))) < 1e-8
            pass = 2;
        elseif k == 50;
            disp('Iterations number exceeded!')
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
    
    if mod(t,10.) == 0
    
    tr = t/10;
        
    %ghn(:,:,(t/10))=hn;
    maxHN(tr) = max(max(abs(hn)));  % Keeps max values
    
    maxIT(tr) = k;
    
    if ismember(t,movieA) == 1
        movieH(:,:,find(movieA==t)) = hm1;
    end
    
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
    
    if L1(tr) <1E-6
        flag = 2;
    end
    
     if mod(t,500.) == 0
         hn = hm1;
         save testmat    
     end
        
    end
    
    if flag == 2;
        break;
    end
    
    % New matrices are old matrices %
    
    hn = hm1;
    
end

save testmat
