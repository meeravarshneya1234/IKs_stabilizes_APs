function deriv = dydt_TT06(t,statevar,Id,p,c)

statevarcell = num2cell(statevar) ;

[V,m,h,j,d,f,f2,fCass,r,s,xs,xr1,xr2,Rbar_ryr,Cai,Cass,CaSR,Nai,Ki] = ...
  deal(statevarcell{:}) ;
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %% 
% %% Compute ionic currents
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%

ENa = p.RTF*log(p.Nao/Nai) ;
EK = p.RTF*log(p.Ko/Ki) ;
ECa = 0.5*p.RTF*log(p.Cao/Cai) ;

EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)/(Ki + p.pKNa*Nai)) ;

INa = c.GNa_*m^3*h*j*(V - ENa) ;

ICa = c.PCa_*d*f*f2*fCass*4*p.F/p.RTF*(V-15)* ...
  (0.25*Cass*exp(2*(V-15)/p.RTF) - p.Cao)/ ...
  (exp(2*(V-15)/p.RTF) - 1) ;

Ito = c.Gto_*r*s*(V-EK) ;

IKs = c.GKs_*xs^2*(V-EKs) ;

IKr = c.GKr_*sqrt(p.Ko/5400)*xr1*xr2*(V-EK) ;

aK1 = 0.1/(exp(0.06*(V-EK-200)) + 1) ;
bK1 = ( 3*exp(2e-4*(V-EK+100)) + exp(0.1*(V-EK-10)) )/ ...
  (exp(-0.5*(V-EK)) + 1) ;
xK1 = aK1/(aK1+bK1) ;
IK1 = c.GK1_*sqrt(p.Ko/5400)*xK1*(V-EK) ;

INCX = c.kNaCa*(p.KmNa^3 + p.Nao^3)^-1*(p.KmCa + p.Cao)^-1* ...
  (p.ksat*exp((p.eta-1)*V/p.RTF) + 1)^-1* ...
  (exp(p.eta*V/p.RTF)*Nai^3*p.Cao - exp((p.eta-1)*V/p.RTF)*p.Nao^3*Cai*p.alpha_ncx) ;

fNaK = 1/(1 + 0.1245*exp(-0.1*V/p.RTF) + 0.0353*exp(-V/p.RTF)) ;
INaK = c.INaK_*p.Ko*Nai*fNaK/( (p.Ko + p.KmKo)*(Nai + p.KmNai) ) ;

IpCa = c.GpCa_*Cai/(Cai + p.KpCa) ;

IpK = p.GpK_*(V-EK)/(exp(-(V-25)/5.98) + 1) ;

INab = c.GNab_*(V-ENa) ;
ICab = c.GCab_*(V-ECa) ;

Iion = INa + ICa + Ito + IKs + IKr + IK1 + INCX + INaK + ...
  IpCa + IpK + INab + ICab ;
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %% 
% %% compute rate constants to update gates
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
% % Na current
minf = 1/(exp((-56.86-V)/9.03) + 1)^2 ;
am = 1/(exp((-60-V)/5) + 1) ;
bm = 0.1/(exp((V+35)/5) + 1) + 0.1/(exp((V-50)/200) + 1) ;
taum = am*bm ;

dmdt = (minf-m)/taum ;

hinf = 1/(exp((V+71.55)/7.43) + 1)^2 ;
jinf = 1/(exp((V+71.55)/7.43) + 1)^2 ;

if V >= -40
  ah = 0 ;
  bh = 0.77/(0.13*(exp((-V-10.66)/11.1) + 1)) ;
  aj = 0 ;
  bj = 0.6*exp(0.057*V)/(exp(-0.1*(V+32)) + 1) ;
else
  ah = 0.057*exp((-V-80)/6.8) ;
  bh = 2.7*exp(0.079*V)+3.1e5*exp(0.3485*V) ;
  aj = (-2.5428e4*exp(0.2444*V)-6.948e-6*exp(-0.04391*V))*(V+37.78)/ ...
    (exp(0.311*(V+79.23))+1) ;
  bj = 0.02424*exp(-0.01052*V)/(exp(-0.1378*(V+40.14)) + 1);
end

tauh = 1/(ah+bh) ;
tauj = 1/(aj+bj) ;

dhdt = (hinf-h)/tauh ;
djdt = (jinf-j)/tauj ;

%%%%%%%%%%%%%%%%%
% % L-type Ca current
dinf = 1/(exp(-(V+8)/7.5) + 1) ;
ad = 1.4/(exp(-(V+35)/13) + 1) + 0.25 ;
bd = 1.4/(exp((V+5)/5) + 1) ;
gd = 1/(exp((50-V)/20) + 1) ;      % gammad in paper

taud = (ad*bd + gd) ;

finf = 1/(exp((V+20)/7) + 1) ;
af = 1102.5*exp(-(V+27)^2/225) ;
bf = 200/(1 + exp((13-V)/10)) ;
gf = 180/(1 + exp((V+30)/10)) + 20 ;

tauf = (af + bf + gf) ;

f2inf = 0.67/(exp((V+35)/7) + 1) + 0.33 ;
af2 = 600*exp(-(V+25)^2/170) ;
bf2 = 31/(1 + exp((25-V)/10)) ;
gf2 = 16/(1 + exp((V+30)/10)) ;

tauf2 = (af2 + bf2 + gf2) ;

%%%%%%%%%%%%%%%%%
% % % % all units in uM

fCassinf = 0.6/(1 + (Cass/50)^2) + 0.4 ;

taufCass = 80/(1 + (Cass/50)^2) + 2 ;

dddt = (dinf-d)/taud ;
dfdt = (finf-f)/tauf ;
df2dt = (f2inf-f2)/tauf2 ;
dfCassdt = (fCassinf-fCass)/taufCass ;

% %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%
% % % % % % Need to make sure this is working correctly
% kfCa = (fCainf < fCa || V < -60) ; 
% dfCadt = kfCa*(fCainf-fCa)/taufCa ;

%%%%%%%%%%%%%%%%%
% % transient outward current
rinf = 1/(exp((20-V)/6) + 1) ;
taur = (9.5*exp(-(V+40)^2/1800) + 0.8) ;

switch (p.celltype)
  case 'endo'
    sinf = 1/(exp((V+28)/5) + 1) ;
    taus = (1000*exp(-(V+67)^2/1000) + 8) ;
  case 'mid'
    sinf = 1/(exp((V+20)/5) + 1) ;
    taus = (85*exp(-(V+45)^2/320) + 5/(exp((V-20)/5) + 1) + 3) ;
  otherwise    % % % epicardial
    sinf = 1/(exp((V+20)/5) + 1) ;
    taus = (85*exp(-(V+45)^2/320) + 5/(exp((V-20)/5) + 1) + 3) ;
end

drdt = (rinf-r)/taur ;
dsdt = (sinf-s)/taus ;

%%%%%%%%%%%%%%%%%
% % slow delayed rectifier current
xsinf = 1/(exp(-(V+5)/14) + 1) ;
axs = 1400/sqrt(exp(-(V-5)/6) + 1) ;
bxs = 1/(exp((V-35)/15) + 1) ;
tauxs = (axs*bxs + 80) ;

dxsdt = (xsinf-xs)/tauxs ;

%%%%%%%%%%%%%%%%%
% % rapid delayed rectifier current
xr1inf = 1/(exp(-(V+26)/7) + 1) ;
axr1 = 450/(exp(-(V+45)/10) + 1) ;
bxr1 = 6/(exp((V+30)/11.5) + 1) ;
tauxr1 = (axr1*bxr1) ;

xr2inf = 1/(exp((V+88)/24) + 1) ;
axr2 = 3/(exp(-(V+60)/20) + 1) ;
bxr2 = 1.12/(exp((V-60)/20) + 1) ;
tauxr2 = (axr2*bxr2) ;

dxr1dt = (xr1inf-xr1)/tauxr1 ;
dxr2dt = (xr2inf-xr2)/tauxr2 ;

% %%%%%%%%%%%%%%%%%%%%%%%%%
% %% 
% %% CICR and other ion transport processes
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%

% % % Units:  uM/ms
Ileak = c.Vleak*(CaSR - Cai) ;
Iup = c.Iup_/((p.Kmup/Cai)^2 + 1) ;

kcasr = p.maxsr - (p.maxsr - p.minsr)/(1 + (p.EC_ryr/CaSR)^2) ;
k1_ryr = p.k1_ryr_prime/kcasr ;
k2_ryr = p.k2_ryr_prime*kcasr ;

O_ryr = k1_ryr*Cass^2*Rbar_ryr/(p.k3_ryr + k1_ryr*Cass^2) ;
Irel = c.Vrel*O_ryr*(CaSR - Cass) ;

Ixfer = p.Vxfer*(Cass - Cai) ;

dRbar_ryr = -k2_ryr*Cass*Rbar_ryr + p.k4_ryr*(1 - Rbar_ryr) ;

% % % % These not used but should be equivalent
% % % % Rapid buffer scaling factors
Bi = (1 + p.Bufc*p.Kbufc/(p.Kbufc + Cai)^2)^-1 ;
Bss = (1 + p.Bufss*p.Kbufss/(p.Kbufss + Cass)^2)^-1 ;
BSR = (1 + p.BufSR*p.KbufSR/(p.KbufSR + CaSR)^2)^-1 ;

dVdt = -(Id + Iion) ;

dCai = Bi*(-(IpCa + ICab - 2*INCX)*1e6*p.Acap/(2*p.F*p.Vmyo) + ...
  (p.VSR/p.Vmyo)*(Ileak - Iup) + Ixfer) ;
dCass = Bss*(-ICa*1e6*p.Acap/(2*p.F*p.Vss) + p.VSR/p.Vss*Irel - p.Vmyo/p.Vss*Ixfer) ;
dCaSR = BSR*(Iup - Ileak - Irel) ;

% dCai = (-(IpCa + ICab - 2*INCX)*1e6*p.Acap/(2*p.F*p.Vmyo) + ...
%   (p.VSR/p.Vmyo)*(Ileak - Iup) + Ixfer) ;
% dCass = (-ICa*1e6*p.Acap/(2*p.F*p.Vss) + p.VSR/p.Vss*Irel - p.Vmyo/p.Vss*Ixfer) ;
% dCaSR = (Iup - Ileak - Irel) ;

dNai = -(INa + 3*INCX + 3*INaK + INab)*1e6*p.Acap/(p.F*p.Vmyo) ;
dKi = -(Ito + IKs + IKr + IK1 - 2*INaK + IpK)*1e5*p.Acap/(p.F*p.Vmyo) ;

deriv = [dVdt;dmdt;dhdt;djdt;dddt;dfdt;df2dt;dfCassdt;drdt;dsdt; ...
  dxsdt;dxr1dt;dxr2dt;dRbar_ryr;dCai;dCass;dCaSR;dNai;dKi] ;

return

