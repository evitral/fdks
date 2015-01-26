
%---------------------------------------------------------------%
%                      Function fm for FDKS                     %
%---------------------------------------------------------------%

% Vji(1) = j
% Vji(2) = jp
% Vji(3) = jpp
% Vji(4) = jm
% Vji(5) = jmm
% Vji(6) = i
% Vji(7) = ip
% Vji(8) = ipp
% Vji(9) = im
% Vji(10)= imm

function fm = fm(deck,hm,hn,Vji)

fm = 0.5*(deck(1)*(hm(Vji(1),Vji(7))-2*hm(Vji(1),Vji(6))+hm(Vji(1),Vji(9)))... %delta
    -deck(2)*(hm(Vji(2),Vji(6))-2*hm(Vji(1),Vji(6))+hm(Vji(4),Vji(6))) ... %zeta
    +deck(3)*(hm(Vji(1),Vji(7))-hm(Vji(1),Vji(9))+hn(Vji(1),Vji(7))... %eta
        -hn(Vji(1),Vji(9)))^2 ...
    -deck(4)*(hm(Vji(2),Vji(6))-hm(Vji(4),Vji(6))+hn(Vji(2),Vji(6))... %psi
    -hn(Vji(4),Vji(6)))^2 ...
    +deck(5)*(hm(Vji(2),Vji(7))+hm(Vji(4),Vji(7))... %omega
        -2*(hm(Vji(1),Vji(7))+hm(Vji(2),Vji(6))-2*hm(Vji(1),Vji(6))...
        +hm(Vji(4),Vji(6))+hm(Vji(1),Vji(9)))+hm(Vji(2),Vji(9))...
        +hm(Vji(4),Vji(9)))...
    +deck(6)*(hm(Vji(3),Vji(6))-4*(hm(Vji(2),Vji(6))+hm(Vji(4),Vji(6)))... %phi
        +6*hm(Vji(1),Vji(6))+hm(Vji(5),Vji(6))));
    
end
