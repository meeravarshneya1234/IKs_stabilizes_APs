function deriv = dydt_Fox(t,statevar,Id,p,c)

statevarcell = num2cell(statevar) ;

[V,Cai,CaSR,f,d,m,h,j,fCa,xKr,xKs,xto,yto] = ...
  deal(statevarcell{:}) ;

% Reversal potentials
ENa = p.RTF*log(p.Nao/p.Nai) ;
EK = p.RTF*log(p.Ko/p.Ki) ;
EKs = p.RTF*log((p.Ko+0.01833*p.Nao)/(p.Ki+0.01833*p.Nai)) ;
ECa = 0.5*p.RTF*log(p.Cao/Cai) ;

%% compute ionic currents

% Na currents
INa = c.GNa_*m^3*h*j*(V - ENa);

INab = c.GNab_*(V - ENa) ;

% Ca currents
ICa_ = c.PCa_/p.Cm*4*p.F/p.RTF*V*(Cai*exp(2*V/p.RTF) - 0.341*p.Cao)/ ...
  (exp(2*V/p.RTF) - 1) ;
ICa = ICa_*d*f*fCa ;

ICab = c.GCab_*(V - ECa) ;

IpCa = c.IpCa_*Cai^2/(Cai^2 + p.KmpCa^2) ;

% K currents
K1inf = 1/(exp(1.62/p.RTF*(V - EK)) + 2) ;
IK1 = c.GK1_*K1inf*p.Ko/(p.Ko + p.KmK1)*(V - EK) ;

R_V = 1/(2.5*exp(0.1*(V+28)) + 1) ;
IKr = c.GKr_*R_V*xKr*sqrt(p.Ko/4000)*(V - EK) ;
  
IKs = c.GKs_*xKs^2*(V - EKs) ;

Ito = c.Gto_*xto*yto*(V - EK) ;
  
Kp = 1/(1+exp((7.488-V)/5.98));
IKp = c.GKp_*Kp*(V - EK) ;

ICaK = c.PCaK/p.Cm*f*d*fCa/(1 + ICa_/p.ICahalf)* ...
  p.F/p.RTF*V*(p.Ki*exp(V/p.RTF) - p.Ko)/(exp(V/p.RTF) - 1) ;

% Pumps and transporters
sigma = (exp(p.Nao/67300) - 1)/7 ; 
fNaK = 1/(1 + 0.1245*exp(-0.1*V/p.RTF) + 0.0365*sigma*exp(-V/p.RTF)) ;

INaK = c.INaK_*fNaK*p.Ko/(p.Ko + p.KmKo)*(1 + (p.KmNai/p.Nai)^1.5)^-1 ;

INCX = c.kNaCa*(p.KmNa^3 + p.Nao^3)^-1*(p.KmCa + p.Cao)^-1*(p.ksat*exp((p.eta-1)*V/p.RTF) + 1)^-1* ...
  (exp(p.eta*V/p.RTF)*p.Nai^3*p.Cao - exp((p.eta-1)*V/p.RTF)*p.Nao^3*Cai) ;

Iion = INa + INab + ICa + ICab + IpCa + IK1 + IKr + IKs + Ito + ...
  IKp + ICaK + INaK + INCX  ;

%% compute rate constants to update gates
am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))) ;
bm = 0.08*exp(-V/11); 

ah = 0.135*exp(-(V+80)/6.8) ;
bh = 7.5/(exp(-0.1*(V+11)) + 1) ;
  
aj = 0.175*exp(-(V+100)/23)/(exp(0.15*(V+79)) + 1) ;
bj = 0.3/(exp(-0.1*(V+32)) + 1) ;
  
dmdt = am - (am+bm)*m ;
dhdt = ah - (ah+bh)*h ;
djdt = aj - (aj+bj)*j ;

finf = 1/(exp((V+12.5)/5) + 1) ;
tauf = 30 + 200/(exp((V+20)/9.5) + 1) ;
  
dinf = 1/(exp(-(V+10)/6.24) + 1) ;
taud = 1/(0.25*exp(-0.01*V)/(exp(-0.07*V) + 1) + ...
  0.07*exp(-0.05*(V+40))/(exp(0.05*(V+40)) + 1) ) ; 
  
fCainf = 1/((Cai/p.KmfCa)^3 + 1) ;
taufCa = 30 ;
  
dddt = (dinf - d)/taud ;
dfdt = (finf - f)/tauf ;
dfCadt = (fCainf - fCa)/taufCa ;

xKrinf = 1/(exp(-2.182-0.1819*V) + 1) ;
tauxKr = 43 + 1/(exp(-5.495 + 0.1691*V) + exp(-7.677 - 0.0128*V)) ;

dxKrdt = (xKrinf - xKr)/tauxKr ;
  
xKsinf = 1/(exp(-(V-16)/13.6) + 1) ;
tauxKs = 1/ ...
  (7.19e-5*(V-10)/(1 - exp(-0.148*(V-10))) + ...
  1.31e-4*(V-10)/(exp(0.0687*(V-10)) - 1)) ;

dxKsdt = (xKsinf - xKs)/tauxKs ;

axto = 0.04516*exp(0.03577*V) ;
bxto = 0.0989*exp(-0.06237*V) ;
  
ayto = 0.005415*exp(-(V+33.5)/5)/(0.051335*exp(-(V+33.5)/5) + 1) ;
byto = 0.005415*exp((V+33.5)/5)/(0.051335*exp((V+33.5)/5) + 1) ;
  
dxtodt = axto - (axto+bxto)*xto ;
dytodt = ayto - (ayto+byto)*yto ;

%% Fluxes for ionic balance
% % Rapid buffer scaling factors
Bi = (1 + p.CMDNtot*p.KmCMDN/(p.KmCMDN + Cai)^2)^-1 ;
BSR = (1 + p.CSQNtot*p.KmCSQN/(p.KmCSQN + CaSR)^2)^-1 ;

% % SR Ca release, uptake, and leak
gamma = 1/((2000/CaSR)^3 + 1) ;
Jrel = c.Prel*f*d*fCa*(gamma*CaSR - Cai)/(1.65*exp(V/20) + 1) ;

Jup = c.Vup/((p.Kmup/Cai)^2 + 1) ;
Jleak = c.Pleak*(CaSR - Cai) ;

%% Compute additional derivatives

dVdt = -(Id + Iion) ;

dCai = Bi*(Jrel + Jleak - Jup - p.Acap*p.Cm/(2*p.Vmyo*p.F)* ...
  (ICa + ICab + IpCa - 2*INCX) ) ;

dCaSR = BSR*(Jup - Jleak - Jrel)*p.Vmyo/p.VSR ;

deriv = [dVdt;dCai;dCaSR;dfdt;dddt;dmdt;dhdt;djdt;dfCadt; ...
  dxKrdt;dxKsdt;dxtodt;dytodt] ;


return


