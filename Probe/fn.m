
%---------------------------------------------------------------%
%                      Function fn for FDKS                     %
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

function fn = fn(deck,hn,Vyx)

fn = 0.5*(deck(1)*(hn(Vyx(1),Vyx(7))-2*hn(Vyx(1),Vyx(6))+hn(Vyx(1),Vyx(9)))...    
    -deck(2)*(hn(Vyx(2),Vyx(6))-2*hn(Vyx(1),Vyx(6))+hn(Vyx(4),Vyx(6))) ...
    +deck(5)*(hn(Vyx(2),Vyx(7))+hn(Vyx(4),Vyx(7))...
        -2*(hn(Vyx(1),Vyx(7))+hn(Vyx(2),Vyx(6))-2*hn(Vyx(1),Vyx(6))...
        +hn(Vyx(4),Vyx(6))+hn(Vyx(1),Vyx(9)))+hn(Vyx(2),Vyx(9))...
        +hn(Vyx(4),Vyx(9)))...
    +deck(6)*(hn(Vyx(3),Vyx(6))-4*(hn(Vyx(2),Vyx(6))+hn(Vyx(4),Vyx(6)))...
        +6*hn(Vyx(1),Vyx(6))+hn(Vyx(5),Vyx(6))));
    
end