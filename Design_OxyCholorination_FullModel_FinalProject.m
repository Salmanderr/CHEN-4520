% Cp and Hrxn Data imported from Excel Sheet:

close all 
clear
clc 

%{
For this example, the heat capacities' coefficients are not used in 
the MatLab code since the calculations are already done in Excel. 
To arrive to the same computed values in Cp533, integerate: 

Cp = a + b*( (c/Tk1)/sinh(c/TK1) ) + d*( (e/TK1)/cosh(e/TK1) )  ; 

Temperatures are in absolute units. (K in this example). 
Values imported from Perry's Handbook. 
%}
CpCoeff = [33380,94790,1596,55100,740.800000000000;
    29157,9048,2093.80000000000,-107,120;
    29103,10040,2526.50000000000,9356,1153.80000000000;
    29142,9176,949,10030,425;
    29370,34540,1428,26400,588;
    33363,26790,2610.50000000000,8896,1169;
    65271,112540,1737.60000000000,87800,795.450000000000;
    66554,112570,1545.40000000000,97196,717.040000000000];


Cp533 = [ 224.26; % kJ/mol-K
          53.46 ; 
          44.94 ; 
          26.370;
          86.80 ;
          137.12;
          293.51;
          253.93]; % values computed - assuming constant Heat Capacity

 Hrxn533 = [-239.230000000000;  % reaction 1
     -161.570000000000;         % reaction 2
     -1321.10000000000;         % reaction 3
     -116.160000000000];        % reaction 4

 Hr1 = Hrxn533(1,1); % kJ/mol
 Hr2 = Hrxn533(2,1); 
 Hr3 = Hrxn533(3,1); 
 Hr4 = Hrxn533(4,1); 
 
 
% Define Gas Constant:
global Rgas Rgas2 

Rgas = 8.314; % kJ/kgmol/K
Rgas2 = 0.082057; % atm*m3/kgmol/K

% Molecular Weights of Species:

MwHCl     = 36.458; % kg/kgmol 
MwC2H4    = 28.05; 
MwO2      = 16; 
MwCl2     = 70; 
MwH2O     = 18; 
MwC2H4Cl2 = 98.96; 
MwC2H3Cl3 = 133.4; 
MwCO2     = 44.01; 

% Assigning Heat Capacities for each species at 533 K:

CpC2H4 = Cp533(1,1)*1000; %kJ/mol-K
CpHCl  = Cp533(2,1)*1000; % converted to kJ/kmol-k
CpO2   = Cp533(3,1)*1000; 
CpCl2  = Cp533(4,1)*1000; 
CpCO2  = Cp533(5,1)*1000; 
CpH2O  = Cp533(6,1)*1000; 
CpC2H4Cl2  = Cp533(7,1)*1000; 
CpC2H3Cl3 = Cp533(8,1)*1000; 
CpCoolent = 290.3;  % 290.3 kJ/kmol-K (or 3.221 kJ/kg-k)
% Found from HYSYS simulation for Dowtherm (C4H10O)

% Reactor Parameters: 

Dr = 0.767; % diameter, m ---> GUESSED (Or from HYSYS?) 1.6765
Ntube = 750; % number of tubes 
Lr = 9; % length, m   ---> GUESSED (Or from HYSYS?)
Ar = pi*Dr^2/4; % x-sectional area, m2
Vr = Ar*Lr; % volume, m3
U  = 300; % W/m^2-K
a  = 4/Dr; % (1/m) heat transfer area 
% Feed Conditions: 

FeedC2H4 = 340.87  ; % kgmole/hr ---> All feed conditions are given 
FeedHCl  = 603 ; 
FeedO2   = 170 ; 
FeedCl2  = 2.233  ; 
FeedH2O  = 0      ; 
FeedCO2  = 0      ; 
FeedC2H4Cl2 = 0   ; 
FeedC2H3Cl3 = 0   ; 

FeedCoolent = 22.19; %  (22.19 kmole/hr converted from 2000 kg/hr in HYSYS)
Feedin_tot =  sum( [  FeedC2H4; FeedHCl; FeedO2; FeedCl2; FeedH2O; FeedCO2; ...
    FeedC2H4Cl2; FeedC2H3Cl3] ); 
%{
% Previos Condtions
    T = 520, Tc = 180C, Dr = 1.67m,Lr=3m,Nt=5000
       500         180        0.967   3    5000
%}
FeedP = 2000/101; % kPa --> given (atm /101.325)
% Tinput = input('Guess temperature in degC:');
FeedT = 523             ; % K   --> Guessed
Tc0    = 180 + 273.15    ; % K 
%{
% Feed Partial Pressures: 

% pC2H4 = FeedC2H4/Feedin_tot * FeedP  ; % atm
% pHCl  = FeedHCl/Feedin_tot * FeedP   ; 
% pO2   = FeedO2/Feedin_tot * FeedP    ; 
% pCl2  = FeedCl2/Feedin_tot * FeedP   ; 
% pH2O  = FeedH2O/Feedin_tot * FeedP   ; 
% pCO2  = FeedCO2/Feedin_tot * FeedP   ; 
% pC2H4Cl2 = FeedC2H4CL2/Feedin_tot * FeedP; 
% pC2H3Cl3 = FeedC2H3CL3/Feedin_tot * FeedP; 
%}

% Catalyst Particles and Packing: 


Dp = 1/8*Dr; % particle diameter, m ---> Heuristic?? 
Eps = 0.5 ; % void fraction       ---> GIVEN
rhoCat_Bulk    = 1975; % kg cat/m^3  ---> GIVEN
rhoCat_Int = rhoCat_Bulk/eps; % kg cat/m^3 --> intrinsic density
% Coolant kinematic viscosity estimate from Hysys case
nu = 1.534e-7; % m2/s

% Mass flux
M0 = FeedC2H4 * MwC2H4 + FeedHCl*MwHCl + FeedO2 * MwO2 +FeedCl2*MwCl2 + MwH2O * MwH2O + ...
    FeedCO2*MwCO2 + FeedC2H4Cl2 * MwC2H4Cl2 +  FeedC2H3Cl3*MwC2H3Cl3 ; 

G0 = M0/Ar/3600; %kg/m^2/sec (mass flux in Ergun equation)

% % Initial Enthalpy. 
% H0 = FeedC2H4 * CpC2H4 + FeedHCl*CpHCl + FeedO2 * CpO2 + FeedCl2*CpCl2 + FeedH2O * CpH2O + ...
%     FeedCO2*CpCO2 + FeedC2H4Cl2 * CpC2H4Cl2 + FeedC2H3Cl3 * CpC2H3Cl3 ; 

% Initial conditions: 

y0 = [ FeedC2H4; FeedHCl; FeedO2; FeedCl2; FeedH2O; FeedCO2; ...
    FeedC2H4Cl2; FeedC2H3Cl3; FeedT; Tc0; FeedP];
FeedP0 = FeedP;
FeedT0 = FeedT; 

% Define Solution Span: 

vspan = linspace(0, Vr, 200); 

% Define Anonymous Function

PBRanon = @(v,y) PBRfullModel( v, y, FeedC2H4, FeedHCl, FeedO2, FeedCl2,...
    FeedH2O, FeedCO2,FeedC2H4Cl2, FeedC2H3Cl3, FeedCoolent, Eps,...
    CpC2H4, CpHCl, CpO2, CpCl2, CpCO2, ...
    CpH2O, CpC2H4Cl2, CpC2H3Cl3, CpCoolent, Hr1, Hr2, Hr3, Hr4,U, a, ...
    MwC2H4, MwHCl,MwO2, MwCl2, MwH2O, MwCO2, MwC2H4Cl2, MwC2H3Cl3, ...
    Dp, Ar, nu, G0, Feedin_tot, FeedP0,FeedT0, Dr, Ntube);

%{
 PBRsimplified( v, y, FeedC2H4, FeedHCl, FeedO2, FeedCl2,...
    P, Eps, CpC2H4, CpHCl, CpO2, CpCl2, CpCO2, ...
    CpH2O, CpC2H4Cl2, CpC2H3Cl3, Hr1, Hr2, Hr3, Hr4)
%}

% options = odeset('Mass',M);
    
[ v , ysoln ] = ode23tb(PBRanon,vspan,y0);


% Unpack Soltuions: 
for i =1:35
    for j =1:11
   if  ysoln(i,j) < 0
       ysoln(i,j) = 0;
   end 
    end
end


MflowC2H4 = ysoln(:,1); 
MflowHCl  = ysoln(:,2); 
MflowO2   = ysoln(:,3);
MflowCl2  = ysoln(:,4); 
MflowH2O  = ysoln(:,5); 
MflowCO2  = ysoln(:,6); 
MflowC2H4Cl2 = ysoln(:,7);
MflowC2H3Cl3 = ysoln(:,8);
T            = ysoln(:,9);
Tc           = ysoln(:,10); 
P            = ysoln(:,11)*Ntube; 
% H            = ysoln(:,12); 

Mflow_tot = ( MflowC2H4 + MflowCl2 +  MflowO2 + MflowH2O + MflowCO2 + MflowC2H4Cl2 + MflowC2H3Cl3); 

Conv = (FeedC2H4 - MflowC2H4)/FeedC2H4 *100;% Conversion, c2h4 basis, in %
% create plots
figure(1) % Temperature Profile
plot(v,T,'k-')
grid off
xlabel('Reactor Volume - L','interpreter','latex')
ylabel('Reactor Temperature - K','interpreter','latex')
% title('Temperature Profile of the Reactor','interpreter','latex')

figure(2) % Molar Flows Profile
plot(v,MflowC2H4,'b--',v,MflowHCl,'m--',v,MflowO2,'g--',v, MflowCl2,'c--',v,MflowH2O,'--', v,MflowC2H4Cl2,'y--', v,MflowC2H4Cl2, 'k--')
grid off
xlabel('Reactor Volume - L','interpreter','latex')
ylabel('Molar Flow Rate - kgmol/h','interpreter','latex')
% title('Molar Flow Profile','interpreter','latex')
legend('$C_{2}H_{4}$','$HCl$','$O_{2}$','$Cl_{2}$','$H_{2}O$','$C_{2}H_{4}Cl_{2}$','$C_{2}H_{4}Cl_{3}$','Location','northeastoutside','interpreter','latex')

figure(3)
plot(v,P,'*');
xlabel('Reactor Volume - L','interpreter','latex');
ylabel('Pressure drop - kPa','interpreter','latex'); 
% title('Pressure Drop in PBR','interpreter','latex');
figure(4) 
plot(v,Tc)
ylabel('Dowtherm Temperature - K','interpreter','latex');
xlabel('Reactor Volume - L','interpreter','latex'); 
% title('Temperature Profile of Dowtherm','interpreter','latex'); 

Vflow = 36.4220548595289; % m^3/hr --> hysys
disp(Conv(end));
disp(v(end)); 
disp(v(end)*Eps*rhoCat_Bulk/1000); % 1000 to convert form L to m^3 
Residence = v(end)/Vflow*3600;
disp(Residence); 
disp(Lr/Dr); 

%%%%%%%%%%%%%%%%%%%
% YIELD OF PRODUCT%
%%%%%%%%%%%%%%%%%%%

pC2H4 = MflowC2H4/Feedin_tot .* P   ; 
pHCl  = MflowHCl/Feedin_tot  .* P   ; 
pO2   = MflowO2/Feedin_tot   .* P   ; 
pCl2  = MflowCl2/Feedin_tot  .* P   ; 
pH2O  = MflowH2O/Feedin_tot  .* P   ; 
pCO2  = MflowCO2/Feedin_tot  .* P   ; 
pC2H4Cl2 = MflowC2H4Cl2/Feedin_tot .* P; 

r1 = 10^(4.2) * exp( -40.1e03/(Rgas*(T) ) ) .* pC2H4 .* pCl2.^(0.5)    * 1000     ; 
r2 = 10^(13.21) * exp( -128.04e03/(Rgas*(T)) ) .* pC2H4Cl2 .* pCl2.^(0.5)  * 1000    ; 
r3 = 10^(6.78) * exp( -112e03/(Rgas*(T)) )   .* pO2 .* pCl2.^(0.5) .* pC2H4  * 1000   ; 
r4 = ( 1000 * exp( 17.13 - 13000/(1.987*(T)) ) ) ...
    /    ( 1 * exp( 5.4 + 16000/(1.987*(T)) ) ) .* pO2 .* pCl2.^(-1)    * 1000  ;

% Yield = (r1(end) - r2(end)) / (2*r1(end) + 4*r3(end) + 2 * r4(end) );
Yield = (r1(end) - r2(end)) / ( 3.5 * r1(end) + 1.5 * r2(end) + 4 * r3(end) + 5 * r4(end) );
fprintf('Yield of EDC is %0.2f',Yield)