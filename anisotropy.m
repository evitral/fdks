clear all
clear variables

Kt = 5.0;
j = 1;
amu = 4;

for theta = 0:0.01:1.21
    
    s = sin(theta);
    c = cos(theta);

    mut = 2*s^2-c^2-(amu*s*c)^2;
    nut = -c^2;

    qc = (abs(mut)/(2*Kt))^0.5;

    delta = 0.75*(abs(mut)-abs(nut))*qc^2;

    tt(j) = theta;
    A(j) = delta/(Kt*qc^4);

    j = j+1;
    
end