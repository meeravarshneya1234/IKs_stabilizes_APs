function deriv = dydt_Hund(t,statevar,Id,p,c)

statevarcell = num2cell(statevar) ;

[V,Cai,Cass,CaNSR,CaJSR,Nai,Ki,Cli,m,h,j,hL, ...
  d,dp,f,f2,fCa,fCa2,xKr,xs1,xs2,ato,ito1,ito2,aa,ro,ri,CAMK_trap] = ...
  deal(statevarcell{:}) ;

% % CaMK
CAMK_bound = p.CaMK0*(1-CAMK_trap)/(1+p.KmCaM/Cass) ;
CAMK_a = CAMK_bound + CAMK_trap ;

dCAMK_trap = 0.05*CAMK_a*CAMK_bound - 6.8e-4*CAMK_trap ;

% Reversal potentials
ENa = p.RTF*log(p.Nao/Nai) ;
EK = p.RTF*log(p.Ko/Ki) ;
EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)/(Ki + p.pKNa*Nai)) ;
ECa = 0.5*p.RTF*log(p.Cao/Cai) ;
ECl = -p.RTF*log(p.Clo/Cli) ;

%% compute ionic currents

% Na currents
INa = c.GNa_*m^3*h*j*(V - ENa) ;
INaL = c.GNaL_*m^3*hL*(V - ENa) ;

% Ca currents

ICa_ = c.PCa_*4*p.F/p.RTF*(V-15)* ...
  (p.gamma_Cai*Cass*exp(2*(V-15)/p.RTF) - p.gamma_Cao*p.Cao)/ ...
  (exp(2*(V-15)/p.RTF) - 1) ;
ICa = ICa_*d^(dp)*f*f2*fCa*fCa2 ;

ICab = c.PCab*4*V*p.F/p.RTF*(p.gamma_Cai*Cai*exp(2*V/p.RTF)- ...
  p.gamma_Cao*p.Cao)/(exp(2*V/p.RTF)-1) ;

IpCa = c.IpCa_*Cai/(Cai + p.KmpCa) ;


% % IK1

ak1 = 1.02/(1+exp(0.2385*(V-EK-59.215))) ;
bk1 = (0.49124*exp(0.08032*(V-EK+5.476))+exp(0.06175*(V-EK-594.31)))/...
    (1+exp(-0.5143*(V-EK+4.753))) ;
IK1 = c.GK1_*sqrt(p.Ko/5.4)*ak1/(ak1+bk1)*(V-EK) ;

% % IK1 updated

% Move this to parameters in main program

RKr = 1/(exp((V+10)/15.4) + 1) ;

IKr = c.GKr_*sqrt(p.Ko/5.4)*xKr*RKr*(V - EK) ;

 % % IKr updated

% Move this to parameters in main program

IKs = c.GKs_*(1+0.6/((3.8e-5/Cai)^(1.4)+1))*xs1*xs2*(V - EKs) ;

 % % IKs updated

Rto = exp(V/300) ;
Ito = c.Gto_*ato^3*ito1*ito2*Rto*(V - EK) ;

 % % Ito updated

 % Chloride currents

CTKCl = p.CTKCl_*(EK-ECl)/((EK-ECl)+87.8251) ;
CTNaCl = p.CTNaCl_*(ENa-ECl)^4/( (ENa-ECl)^4 + 87.8251^4) ;

IClb = c.GClb_*(V - ECl) ;

Ito2 = aa*c.PCl*V*p.F/p.RTF*(Cli-p.Clo*exp(V/p.RTF))/(1-exp(V/p.RTF)) ;


Kp = 1/(1+exp((7.488-V)/5.98)) ;
IKp = c.GKp_*Kp*(V - EK) ;

% Pumps and transporters
sigma_NaK = (exp(p.Nao/67.3) - 1)/7 ; 
fNaK = 1/(1 + 0.1245*exp(-0.1*V/p.RTF) + 0.0365*sigma_NaK*exp(-V/p.RTF)) ;

INaK = c.INaK_*fNaK*p.Ko/(p.Ko + p.KmK_NaK)*(1 + (p.KmNa_NaK/Nai)^2)^-1 ;

% % INaK updated


num = c.kNaCa*(Nai^3*p.Cao*exp(p.eta*V/p.RTF)-p.Nao^3*Cai*exp((p.eta-1)*V/p.RTF)) ;
denom1 = 1+p.ksat*exp((p.eta-1)*V/p.RTF) ;
denom2 = p.KmCao*Nai^3+p.KmNao^3*Cai+p.KmNai^3*p.Cao*(1+Cai/p.KmCai) ;
denom3 = p.KmCai*p.Nao^3*(1+(Nai/p.KmNai)^3)+Nai^3*p.Cao+p.Nao^3*Cai ;
allo = 1/(1 + (p.KmCa_allo/Cai)^2) ;
INCX = allo*num/(denom1*(denom2+denom3)) ;

Iion = INa + INaL + ICa + ICab + IpCa + IK1 + IKr + IKs + Ito + Ito2 + ...
  IKp + INaK + INCX  ;

%% compute rate constants to update gates

% % 1) Na current
if V >= -40
    ah = 0.0;
    aj = 0.0;
    bh = 1/(0.13*(1+exp((V+10.66)/(-11.1)))) ;
    bj = 0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32))) ;
else
    ah = 0.135*exp((80+V)/(-6.8)) ;
    aj = (-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))...
        *(V+37.78)/(1+exp(0.311*(V+79.23))) ;
    bh = 3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V) ;
    bj = 0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))) ;
end

am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))) ;
bm = 0.08*exp(-V/11) ;

dmdt = am - (am+bm)*m ;
dhdt = ah - (ah+bh)*h ;
djdt = aj - (aj+bj)*j ;

% m, h, j updated

% % 2) late Na current
hLinf = 1/(exp((V+91)/6.1) + 1) ;
tauhL = 600 ;

% m gate for late Na current the same
dhLdt = (hLinf - hL)/tauhL ;

% hL updated

% % % 2) Ca current
dinf=1/(1+exp(-(V-4)/6.74)) ;
taud=0.59+0.8*exp(0.052*(V+13))/(1+exp(0.132*(V+13))) ;

dpinf=9-8/(1+exp(-(V+65)/3.4)) ;               % power to which 'd' is raised
taudp = 10 ;                                  % ms
finf=0.7/(1+exp((V+17.12)/7)) + 0.3;
f2inf=0.77/(1+exp((V+17.12)/7)) + 0.23;

tauf= 1/(0.2411*exp(-(0.045*(V-9.6914))^2)+0.0529) ;
tauf2=1/(0.0423*exp(-(0.059*(V-18.5726))^2)+0.0054) ;

fCainf = 0.3/(1-ICa/0.05)+0.55/(1+Cass/0.003)+0.15;%
fCa2inf =1/(1-ICa/0.01) ;

dtaufCa =p.dtaufCa_max*CAMK_a/(p.KmCaMK+CAMK_a) ;
taufCa = dtaufCa + 0.5 + 1/(1+Cass/0.003) ;
taufCa2 = 300/(1+exp((-ICa-0.175)/0.04)) + 125 ;

dddt = (dinf - d)/taud ;
ddpdt = (dpinf - dp)/taudp ;

dfdt = (finf - f)/tauf ;
df2dt = (f2inf - f2)/tauf2 ;

dfCadt = (fCainf - fCa)/taufCa ;
dfCa2dt = (fCa2inf - fCa2)/taufCa2 ;

% d, dp, f, f2, fCa, fCa2 updated
% Still need to figure things out wrt CaMK-dependent variables

xKrinf = 1/(1+exp(-(V+10.085)/4.25)) ;
tauxKr = 1/(6e-4*(V-1.7384)/(1-exp(-0.136*(V-1.7384)))...
    +3e-4*(V+38.3608)/(exp(0.1522*(V+38.3608))-1)) ;

dxKrdt = (xKrinf - xKr)/tauxKr ;

% xKr updated

xsinf = 1/(1+exp(-(V-10.5)/24.7)) ;

tauxs1 = 1/ ...
  ( 7.61e-5*(V+44.6)/(1-exp(-9.97*(V+44.6)))+ ...
  3.6e-4*(V-0.55)/(exp(0.128*(V-0.55))-1) ) ;
tauxs2 = 2*tauxs1 ;

dxs1dt = (xsinf - xs1)/tauxs1 ;
dxs2dt = (xsinf - xs2)/tauxs2 ;

% xs1 xs2 updated

aato = 25*exp((V-40)/25)/(1+exp((V-40)/25)) ;
bato = 25*exp(-(V+90)/25)/(1+exp(-(V+90)/25)) ;

aito1 = 0.03/(1+exp((V+60)/5)) ;
bito1 = 0.2*exp((V+25)/5)/(1+exp((V+25)/5)) ;

aito2 = 0.00225/(1+exp((V+60)/5)) ;
bito2 = 0.1*exp((V+25)/5)/(1+exp((V+25)/5)) ;

datodt = aato - (aato+bato)*ato ;
dito1dt = aito1 - (aito1+bito1)*ito1 ;
dito2dt = aito2 - (aito2+bito2)*ito2 ;

% ato ito1 ito2 updated

% second transient outward current
aainf = 1/(1.0+p.Kmto2/Cass) ;
tauaa = 1 ;

daadt = (aainf - aa)/tauaa ;

% aa updated

% % Intracellular Ca fluxes
% % SR Ca release, uptake, and leak
% 
% 

% % This is correct.  cafac depends on ICa; vg depends on ICa_
cafac=1/(exp((ICa+0.05)/0.015) + 1) ;

vg=1/(exp((ICa_+13)/5) + 1) ;
Jrel = c.Krel_*vg*ro*ri*(CaJSR-Cass) ;

dKmPLB=p.dKmPLBmax*CAMK_a/(p.KmCaMK+CAMK_a) ;
dJup=p.dJupmax*CAMK_a/(p.KmCaMK+CAMK_a) ;
Jup = (dJup+1.0)*c.Vup*Cai/(Cai+p.Kmup-dKmPLB) ;
Jleak = c.Vup*CaNSR/p.CaNSR_max ;

Jdiff = (Cass-Cai)/p.tau_diff ;
Jtr = (CaNSR-CaJSR)/p.tau_transfer ;

% SR release gating variables
% move to different section
droinf = (CaJSR^1.9)/(CaJSR^1.9+(49.28*Cass/(Cass+0.0028))^1.9) ;

dtau_rel = p.dtau_rel_max*CAMK_a/(p.KmCaMK+CAMK_a) ;

roinf = droinf*ICa^2/(ICa^2 + 1) ;
tauro = 3 ;

riinf = 1/(exp((Cass - 4e-4 + 2e-3*cafac)/2.5e-5) + 1) ;
tauri = 3+dtau_rel+(350.0-dtau_rel)/(1.0+exp((Cass-0.003+0.003*cafac)/2e-4)) ;

drodt = (roinf - ro)/tauro ;
dridt = (riinf - ri)/tauri ;

% % % Buffering factors for rapid buffering approximation

BJSR = (1 + p.CSQNtot*p.KmCSQN/(p.KmCSQN + CaJSR)^2)^-1 ;

Bi = 1/(1+ p.CMDNtot*p.KmCMDN/(Cai+p.KmCMDN)^2+ ...
  p.TRPNtot*p.KmTRPN/(Cai+p.KmTRPN)^2) ;

Bss = 1/(1+p.BSRtot*p.KmBSR/(Cass+p.KmBSR)^2+p.BSLtot*p.KmBSL/(Cass+p.KmBSL)^2) ;

% % % Derivatives for voltage and ionic concentrations
dVdt = -(Id + Iion) ;

dNai = -(INa + INaL + 3*INCX + 3*INaK)*p.Acap/(p.Vmyo*p.F) + CTNaCl ;
dKi = -(IKr + IKs + IK1 + IKp - 2*INaK + Ito + 0.5*Id)*p.Acap/(p.Vmyo*p.F) + CTKCl ;
dCli = (IClb + Ito2 + 0.5*Id)*p.Acap/(p.Vmyo*p.F) + CTNaCl + CTKCl ;

dCass = Bss*( Jrel*p.VJSR/p.Vss - Jdiff - ICa*p.Acap/(2*p.Vss*p.F) ) ;
% Bss
% Jrel*VJSR/Vss
% Jdiff
% ICa*Acap/(2*Vss*F)

dCai = Bi*( (Jleak - Jup)*p.VNSR/p.Vmyo + Jdiff*p.Vss/p.Vmyo - ...
  (ICab + IpCa - 2*INCX)*p.Acap/(2*p.Vmyo*p.F) ) ;
% Bi 
% Jleak*VNSR/Vmyo
% Jup 
% Jdiff 
% 
% guyer = guyer
% 


dCaJSR = BJSR*(Jtr - Jrel) ;
dCaNSR = Jup - Jleak - Jtr*p.VJSR/p.VNSR ;

deriv = [dVdt;dCai;dCass;dCaNSR;dCaJSR;dNai;dKi;dCli; ...
  dmdt;dhdt;djdt;dhLdt;dddt;ddpdt;dfdt;df2dt;dfCadt;dfCa2dt; ... 
  dxKrdt;dxs1dt;dxs2dt;datodt;dito1dt;dito2dt;daadt; ...
  drodt;dridt;dCAMK_trap] ;

return


