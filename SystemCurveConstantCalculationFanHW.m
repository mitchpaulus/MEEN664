clc
clear all

syms eta N Phi Psi P

%Coefficients for PLR calculation.
c1 = 0.00153;
c2 = 0.0052;
c3 = 1.1086;
c4 = -0.1164;

%Coefficients for efficiency calculation.
b1 = 2.957;
b2 = -2.896;
b3 = 0.29;

%Coefficients for Psi calculation.
a0 = 3.46;
a1 = 2.286;
a2=-4.104;
a3=-1.53;

Vdesign = 10000; % cfm
VdesignM = Vdesign / (35.315*60);

rhoAir = 1.2; %kg/m^3
Wrated = 2384.06678; %W

Dia = 0.76; %m

Vcfm = 10000:-500:500;

for i = 1:length(Vcfm)
    V = Vcfm(i) / (35.315*60); %m^3/s
    PLR = V / VdesignM;
    Wfraction = c1+c2*PLR+c3*PLR^2+c4*PLR^3;   
    Ws = Wrated*Wfraction;     %W
    
    eqn1 =  Phi == V / (N*Dia^3);
    eqn2 =  Ws == V * P / eta;
    eqn3 =  eta == b1*Phi + b2*Phi^2 + b3*Phi^3;
    eqn4 =  Psi == P / (rhoAir*N^2*Dia^2);
    eqn5 =  Psi == a0 + a1*Phi + a2*Phi^2 + a3*Phi^3;

    eqns = [eqn1, eqn2, eqn3, eqn4, eqn5];

    %Set constraints on solution. 
    assume(N > 0);
    assume(eta >= 0);
    assume(eta <= 1);
    assume(P >= 0);
    assume(Psi >= 0);
    assume(Phi >= 0);
    S = solve(eqns, 'Real', true)


    TotalAnswer(i,1) = S.N;
    TotalAnswer(i,2) = S.eta;
    TotalAnswer(i,3) = S.Phi;
    TotalAnswer(i,4) = S.Psi;
    TotalAnswer(i,5) = S.P;

end


TotalAnswer = double(TotalAnswer)





