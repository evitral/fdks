
%---------------------------------------------------------------%
%                      Function fn for FDKS                     %
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

function fn = fn(deck,hn,Vji)

fn = 0.5*(deck(1)*(hn(Vji(1),Vji(7))-2*hn(Vji(1),Vji(6))+hn(Vji(1),Vji(9)))...
    -deck(2)*(hn(Vji(2),Vji(6))-2*hn(Vji(1),Vji(6))+hn(Vji(4),Vji(6))) ...
    +deck(5)*(hn(Vji(2),Vji(7))+hn(Vji(4),Vji(7))...
        -2*(hn(Vji(1),Vji(7))+hn(Vji(2),Vji(6))-2*hn(Vji(1),Vji(6))...
        +hn(Vji(4),Vji(6))+hn(Vji(1),Vji(9)))+hn(Vji(2),Vji(9))...
        +hn(Vji(4),Vji(9)))...
    +deck(6)*(hn(Vji(3),Vji(6))-4*(hn(Vji(2),Vji(6))+hn(Vji(4),Vji(6)))...
        +6*hn(Vji(1),Vji(6))+hn(Vji(5),Vji(6))));
    
end