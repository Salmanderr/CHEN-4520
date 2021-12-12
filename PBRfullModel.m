function [ outvar ] = PBRfullModel( v, y, FeedC2H4, FeedHCl, FeedO2, FeedCl2,...
    FeedH2O, FeedCO2,FeedC2H4Cl2, FeedC2H3Cl3, FeedCoolent, Eps,...
    CpC2H4, CpHCl, CpO2, CpCl2, CpCO2, ...
    CpH2O, CpC2H4Cl2, CpC2H3Cl3, CpCoolent, Hr1, Hr2, Hr3, Hr4,U, a, ...
    MwC2H4, MwHCl,MwO2, MwCl2, MwH2O, MwCO2, MwC2H4Cl2, MwC2H3Cl3, ...
    Dp, Ar, nu, G0, Feedin_tot, FeedP0,FeedT0,Dr,Ntube)
global Rgas 
%{
    The function takes initial values defined for Oxycholorination reactor
    and solves the associated mass balances and energy balance; temperature
    profile is non-isothermal. The function is based on Prof. David Clough
    MatLab code provided in the course website by Prof. Alan Wymeir and 
    Prof. Wendy Young for Design course in Fall 2022. The code is modified
    to fit the reactor given in HW 8. 
    
    Team Members: 
    Salman Salman, Megan Smith, Maddie, Erast, Abdularhman Mansouri.
%}
    
% Unpack depedent variables: 
% y0 = [ FeedC2H4; FeedHCl; FeedO2; FeedCl2; FeedH2O; FeedCO2; ...
  %   FeedC2H4Cl2; FeedC2H3Cl3; FeedT];
MflowC2H4 = y(1); 
MflowHCl  = y(2); 
MflowO2   = y(3);
MflowCl2  = y(4); 
MflowH2O  = y(5); 
MflowCO2  = y(6); 
MflowC2H4Cl2 = y(7);
MflowC2H3Cl3 = y(8);
T            = y(9); 
Tc           = y(10); 
P            = y(11);


% Define total flow rate:

Feed_tot =  sum( [  MflowC2H4; MflowHCl; MflowO2; MflowCl2; MflowH2O; MflowCO2; ...
    MflowC2H4Cl2; MflowC2H3Cl3] ); 

% Define Partial Pressures in kPa : 
pC2H4 = MflowC2H4/Feed_tot * P   ; 
pHCl  = MflowHCl/Feed_tot  * P   ; 
pO2   = MflowO2/Feed_tot   * P   ; 
pCl2  = MflowCl2/Feed_tot  * P   ; 
pH2O  = MflowH2O/Feed_tot  * P   ; 
pCO2  = MflowCO2/Feed_tot  * P   ; 
pC2H4Cl2 = MflowC2H4Cl2/Feed_tot * P; 
pC2H3Cl3 = MflowC2H3Cl3/Feed_tot * P; 

% Define reaction rates: [source: Lakshmanan et al, 1999]( mol/Liter cat-h)
% converted units to mol/m^3 cat-hr
r1 = 10^(4.2) * exp( -40.1e03/(Rgas*(T) ) ) * pC2H4 * pCl2^(0.5)    * 1000     ; 
r2 = 10^(13.21) * exp( -128.04e03/(Rgas*(T)) ) * pC2H4Cl2 * pCl2^(0.5)  * 1000    ; 
r3 = 10^(6.78) * exp( -112e03/(Rgas*(T)) )   * pO2 * pCl2^(0.5) * pC2H4  * 1000   ; 
r4 = ( 1000 * exp( 17.13 - 13000/(1.987*(T)) ) ) ...
    /    ( 1 * exp( 5.4 + 16000/(1.987*(T)) ) ) * pO2 * pCl2^(-1)    * 1000  ;

% Define differential mole balances on species: 
% (mol/m^3-hr) 
outvar(1) = (-r1-r3) * (1-Eps);                        % C2H4
outvar(2) = (-2*r1-r2-4*r4) * (1-Eps);                 % HCl
outvar(3) = (-3*r3-r4-1/2*r2-1/2*r1) * (1-Eps);        % O2
outvar(4) = (2*r4) * (1-Eps);                          % Cl2
outvar(5) = (2*r4+2*r3+r2+r1) * (1-Eps);               % H2O
outvar(6) = (2*r3) * (1-Eps);                          % CO2
outvar(7) = (r1-r2) * (1-Eps);                         % C2H4Cl2
outvar(8) = (r2) * (1-Eps);                            % C2H3Cl3

% Define Average molecular weight: kg/kgmol

MwAvg =( MflowC2H4 * MwC2H4 + MflowHCl*MwHCl + MflowO2 * MwO2 + MflowCl2*MwCl2 + MflowH2O * MwH2O + ...
    MflowCO2*MwCO2 + MflowC2H4Cl2 * MwC2H4Cl2 +  MflowC2H3Cl3*MwC2H3Cl3 ) / Feed_tot ;




% Define Energy Balance for the reactions: 
% The negative signs accounting for consumption. 
% K/m^3 --> K/L
outvar(9) = ( ( -r1*Hr1 + -r2 * Hr2 + -r3 * Hr3  + -2*r4 * Hr4 ) - U*a*(T-Tc)/ (1-Eps) * (3600) )  / ... 
    ( MflowC2H4 * CpC2H4 + MflowHCl * CpHCl + MflowO2 * CpO2 + MflowCl2 * ...
    CpCl2 + MflowH2O * CpH2O + MflowCO2 * CpCO2 + MflowC2H4Cl2 * CpC2H4Cl2 ...
    + MflowC2H3Cl3 * CpC2H3Cl3 ); 

% Coolant Energy Balance K/m^3 --> K/L
outvar(10) = U*a*(T-Tc)* (3600) /FeedCoolent/CpCoolent/(Dr*Ntube*(1-Eps)); 

% Pressure drop: kPa/L 

rho = MwAvg * P/Rgas/T; % density as function o
% gas viscosity from kinematic viscosity and density
mu = nu * rho; % Pa-s

outvar(11) = -(20*(1-Eps)/Dp/G0*mu+7/4)*(1-Eps)/Eps^3/Dp/rho*G0^2 ...
    /Ar/101325* FeedP0/P * Feed_tot/Feedin_tot * T/FeedT0;

outvar = outvar'; 
