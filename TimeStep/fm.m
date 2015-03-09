
%---------------------------------------------------------------%
%                      Function fm for FDKS                     %
%---------------------------------------------------------------%

% Vji(1) = y
% Vji(2) = yp
% Vji(3) = ypp
% Vji(4) = ym
% Vji(5) = ymm
% Vji(6) = x
% Vji(7) = xp
% Vji(8) = xpp
% Vji(9) = xm
% Vji(10)= xmm

function fm = fm(deck,hm,hn,Vyx)

fm = 0.5*(deck(1)*(hm(Vyx(1),Vyx(7))-2*hm(Vyx(1),Vyx(6))+hm(Vyx(1),Vyx(9)))... %delta
    -deck(2)*(hm(Vyx(2),Vyx(6))-2*hm(Vyx(1),Vyx(6))+hm(Vyx(4),Vyx(6))) ... %zeta
    +deck(3)*(hm(Vyx(1),Vyx(7))-hm(Vyx(1),Vyx(9))+hn(Vyx(1),Vyx(7))... %eta
        -hn(Vyx(1),Vyx(9)))^2 ...
    -deck(4)*(hm(Vyx(2),Vyx(6))-hm(Vyx(4),Vyx(6))+hn(Vyx(2),Vyx(6))... %psi
    -hn(Vyx(4),Vyx(6)))^2 ...  
    +deck(5)*(hm(Vyx(2),Vyx(7))+hm(Vyx(4),Vyx(7))... %omega
        -2*(hm(Vyx(1),Vyx(7))+hm(Vyx(2),Vyx(6))-2*hm(Vyx(1),Vyx(6))...
        +hm(Vyx(4),Vyx(6))+hm(Vyx(1),Vyx(9)))+hm(Vyx(2),Vyx(9))...
        +hm(Vyx(4),Vyx(9)))...
    +deck(6)*(hm(Vyx(3),Vyx(6))-4*(hm(Vyx(2),Vyx(6))+hm(Vyx(4),Vyx(6)))... %phi
        +6*hm(Vyx(1),Vyx(6))+hm(Vyx(5),Vyx(6))));
    
end
