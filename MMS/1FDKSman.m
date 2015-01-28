
%---------------------------------------------------------------%
%              Method of Manufactured Solutions                 %
%  h = ho+hx*sin(ax*pi()*x/L)+hy*cos(ay*pi()*y/L)               %
%   +hxy*sin(ax*pi()*x/L)*cos(ay*pi()*y/L)*exp(b*t)             %
%---------------------------------------------------------------%

clear all
clear variables
clc

%---------------------------------------------------------------%
%                Parameters of the dynamics                     %
%---------------------------------------------------------------%

alpha = 0;         % unclear physical meaning (B&H, alpha = 0)
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

dX = 1;         
dY = 1;         

np = 128;         % number of points per side

L = np-1;

%---------------------------------------------------------------%
%                    Transformed parameters                     %
%---------------------------------------------------------------%

%dtau = 2*F*amu^2*dt/a;
dtau = 0.1;
ntau = 3000;
tmax = ntau*dtau;

Kt = (amu^2*0.54*10^16/T)*exp(-dE*1.16*10^4/T)*exp(-amu^2*s^2/2);
Kt = 5;
mut = 2*s^2-c^2-(amu*s*c)^2;

nuxt = c*(3*s^2-c^2-(amu*s*c)^2);



%---------------------------------------------------------------%
%                      Useful constants                         %
%---------------------------------------------------------------%

Dxx = (c^2-4*s^2+2*amu^2*s^2*(c^2-(2/3)*s^2)+(amu^4/3)*s^4*c^2);

Dxy = 2*(c^2-2*s^2+(amu*s*c)^2); % 2 times Dxy actually

MMS = zeros(np,np);

j=0;

t=100;
b = -1;
ebt = exp(b*t);
ho = 0.01;
hx = 0.01;
hy = 0.01;
hxy = 0.01;
ax = 2;
ay = 2;


for y = 0:dY:L
    
    i = 0  ;
    j = j+1;
    
    coy = cos(ay*pi()*y/L);
    siy = sin(ay*pi()*y/L);
    
    for x = 0:dX:L
        
        i = i+1;
        
        cox = cos(ax*pi()*x/L);
        six = sin(ax*pi()*x/L);
        
        MMS(j,i) = ebt*b*hxy*coy*six+alpha*(ho+hy*coy+hx*six+ebt*coy*six) ...
        -mut*(-ax^2*hx*pi()^2*six/L^2 -ax^2*ebt*hxy*pi()^2*coy*six/L^2) ...
        +c^2*(-ay^2*hy*pi()^2*coy/L^2 -ay^2*ebt*hxy*pi()^2*coy*six/L^2)...
        -nuxt*(ax*hx*pi()*cox/L+ax*ebt*hxy*pi()*cox*coy/L)^2 ...
        +c^3*(-ay*hy*pi()*siy/L-ay*ebt*hxy*pi()*six*siy/L)^2 ...
        +Dxx*(ax^4*hx*pi()^4*six/L^4+ax^4*ebt*hxy*pi()^4*coy*six/L^4)...
        -Dxy*ax^2*ay^2*ebt*hxy*pi()^4*coy*six/L^4 ...
        -c^2*(ay^4*hy*pi()^4*coy/L^4+ay^4*ebt*hxy*pi()^4*coy*six/L^4)...
        +Kt*(ay^4*hy*pi()^4*coy/L^4+ax^4*hx*pi()^4*six/L^4 ...
        +ax^4*ebt*hxy*pi()^4*coy*six/L^4 ...
        +2*ax^2*ay^2*ebt*hxy*pi()^4*coy*six/L^4 ...
        +ay^4*ebt*hxy*pi()^4*coy*six/L^4);
    
    end
end
        


