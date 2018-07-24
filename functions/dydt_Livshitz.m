function deriv = dydt_Livshitz(t,statevar,Id,p,c)

statevarcell = num2cell(statevar) ;

[V,Cai,CaNSR,CaJSR,Nai,Ki,m,h,j,d,f,b,g,xKr,xs1,xs2,Jrel] = ...
  deal(statevarcell{:}) ;

% Reversal potentials

ENa = p.RTF*log(p.Nao/Nai) ;
EK = p.RTF*log(p.Ko/Ki) ;
EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)/(Ki + p.pKNa*Nai)) ;
ECa = 0.5*p.RTF*log(p.Cao/Cai) ;

% %% compute ionic currents

% Na currents
INa = c.GNa_*m^3*h*j*(V - ENa) ;
INab = c.GNab*(V - ENa) ;

% % L-type Ca current
ICa_ = c.PCa_*4*p.F*p.FRT*V* ...
  (p.gamma_Cai*Cai*exp(2*V*p.FRT) -p.gamma_Cao*p.Cao)/ ...
  (exp(2*V*p.FRT) - 1) ;
ICaK_ = c.PCa_K*p.F*p.FRT*V* ...
  (p.gamma_Ki*Ki*exp(V*p.FRT) - p.gamma_Ko*p.Ko)/ ...
  (exp(V*p.FRT) - 1) ;
ICaNa_ = c.PCa_Na*p.F*p.FRT*V* ...
  (p.gamma_Nai*Nai*exp(V*p.FRT) - p.gamma_Nao*p.Nao)/ ...
  (exp(V*p.FRT) - 1) ;
fCa = 1/(Cai/p.KmCa + 1) ;

ICaL = ICa_*d*f*fCa ;
ICaL_K = ICaK_*d*f*fCa ;
ICaL_Na = ICaNa_*d*f*fCa ;

ICab = c.GCab*(V - ECa) ;

IpCa = c.IpCa_*Cai/(Cai + p.KmpCa) ;

% % T-type Ca current 
ICaT = c.GCaT*b^2*g*(V-ECa) ;

% K currents
xK1 = 0.004*(1+exp(0.6987*(V-EK+11.724)))/ ...
  (1+exp(0.6168*(V-EK+4.872))) ;
IK1 = c.GK1_*sqrt(p.Ko/5.4)*(V-EK)/(1+xK1) ;

RKr = 1/(exp((V+9)/22.4) + 1) ;

IKr = c.GKr_*sqrt(p.Ko/5.4)*xKr*RKr*(V - EK) ;

IKs = c.GKs_*(1+0.6/((3.8e-5/Cai)^(1.4)+1))*xs1*xs2*(V - EKs) ;

Kp = 1/(1+exp((7.488-V)/5.98)) ;
IKp = c.GKp_*Kp*(V - EK) ;

% Pumps and transporters
sigma_NaK = (exp(p.Nao/67.3) - 1)/7 ; 
fNaK = 1/(1 + 0.1245*exp(-0.1*V*p.FRT) + 0.0365*sigma_NaK*exp(-V*p.FRT)) ;
INaK = c.INaK_*fNaK*p.Ko/( (p.Ko + p.KmK_NaK)*(1 + (p.KmNa_NaK/Nai)^2) ) ;

INCX = c.kNCX*exp((p.eta-1)*V*p.FRT)*(Nai^3*p.Cao*exp(V*p.FRT)-p.Nao^3*Cai) / ...
  ( 1+p.ksat*exp((p.eta-1)*V*p.FRT)*(Nai^3*p.Cao*exp(V*p.FRT)+p.Nao^3*Cai) ) ;

Iion = INa + INab + ICaL + ICaL_Na + ICaL_K + ICab + ICaT + ...
  IpCa + IKr + IKs + IK1 + IKp + INCX + INaK ;

% %% Compute rate constants to update gates
lambda_na=1-1/(1+exp(-(V+40)/0.024));
ah= lambda_na*0.135*exp(-(80+V)/6.8);
bh= (1-lambda_na)/(0.13*(1+exp((V+10.66)/(-11.1)))) + ...
  lambda_na*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));
hinf = ah/(ah + bh) ;

lambda_na=1-1/(1+exp(-(V+40)/0.024));
aj =  lambda_na*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))* ...
  (V+37.78)/ ...
  (1+exp(0.311*(V+79.23)));

bj= (1-lambda_na)*(0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)))) + ...
 lambda_na*(0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))));
jinf = aj/(aj + bj) ;

am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))) ;
bm = 0.08*exp(-V/11) ;
minf = am/(am + bm) ;

lambda_na=1-1/(1+exp(-(V+40)/0.024));
ah= lambda_na*0.135*exp(-(80+V)/6.8);
bh= (1-lambda_na)/(0.13*(1+exp((V+10.66)/(-11.1)))) + ...
  lambda_na*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));
tauh = 1/( (ah + bh)) ;

aj =  lambda_na*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))* ...
  (V+37.78)/ ...
  (1+exp(0.311*(V+79.23)));
bj= (1-lambda_na)*(0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)))) + ...
 lambda_na*(0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))));
tauj = 1/( (aj + bj)) ;

am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))) ;
bm = 0.08*exp(-V/11) ;
taum = 1/( (am + bm)) ;

dmdt = (minf-m)/taum ;
dhdt = (hinf-h)/tauh ;
djdt = (jinf-j)/tauj ;

% % % 2) Ca current
dinf_0 = 1/(1+exp(-(V+10)/6.24)) ;
dinf_1 = 1/(1+exp(-(V+60)/0.024)) ;
dinf = dinf_0*dinf_1 ;

taud = (1)/(1+exp(-(V+10)/6.24))* ...
  (1-exp(-(V+10)/6.24))/(0.035*(V+10)) ;

finf = 1/(1+exp((V+32)/8))+(0.6)/(1+exp((50-V)/20)) ;

tauf = (1)/(0.0197*exp(-(0.0337*(V+10))^2)+0.02) ;

dddt = (dinf - d)/taud ;
dfdt = (finf - f)/tauf ;

% % T Type current
binf = 1./(1+exp(-(V+14.0)/10.8));

taub = (1)*(3.7+6.1/(1+exp((V+25.0)/4.5))) ;

ginf = 1/(1+exp((V+60.0)/5.6));
lambda_g=1-1/(1+exp(-V/0.0024));
taug = (1)*(lambda_g*(-0.875*V+12.0)+12.0*(1-lambda_g)) ;

dbdt = (binf-b)/taub ;
dgdt = (ginf-g)/taug ;

% % IKr
xKrinf = 1/(1+exp(-(V+21.5)/7.5)) ;
tauxKr = (1)*(1/(0.00138*(V+14.2)/(1-exp(-0.123*(V+14.2)))+ ...
  0.00061*(V+38.9)/(exp(0.145*(V+38.9))-1)) ) ;
dxKrdt = (xKrinf - xKr)/tauxKr ;

% % IKr
% % % Have to decide about this
% % % Two gating variables, same infinity values, taus same by factor of 4
% % % Vary these indepdendently, or dependently???
% % % Start with same variation of both
xsinf = 1/(1+exp(-(V-1.5)/16.7));

tauxs1 = (10000)/(0.719*(V+30)/(1-exp(-0.148*(V+30)))+1.31*(V+30)/ ...
  (exp(0.0687*(V+30))-1));
tauxs2 = 4*tauxs1 ;

dxs1dt = (xsinf - xs1)/tauxs1 ;
dxs2dt = (xsinf - xs2)/tauxs2 ;

% % Intracellular Ca fluxes
% % SR Ca release, uptake, and leak
% 
Jrelinf = c.alpha_rel*p.beta_tau*ICaL/((p.Krel_inf/CaJSR)^p.hrel + 1) ;
if (abs(Jrelinf) < 1e-5)
  Jrelinf = 0 ;
end
tau_rel = (1)*p.beta_tau/(p.Krel_tau/CaJSR + 1) ;
if (tau_rel < 0.1)
  tau_rel = 0.1 ;
end
  
dJreldt = - (Jrelinf + Jrel)/tau_rel ;
% if (abs(dJreldt) > 100)
%   dJreldt = 0 ;
%   Jrel = -Jrelinf ;
% end

Jserca = c.Vserca*(Cai/(Cai+p.Kmserca) - CaNSR/p.CaNSR_max ) ;

Jtr = (CaNSR-CaJSR)/p.tau_transfer ;

% % % Buffering factors for rapid buffering approximation

BJSR = (1 + p.CSQNtot*p.KmCSQN/(p.KmCSQN + CaJSR)^2)^-1 ;

Bi = 1/(1+ p.CMDNtot*p.KmCMDN/(Cai+p.KmCMDN)^2+ ...
  p.TRPNtot*p.KmTRPN/(Cai+p.KmTRPN)^2) ;

% % % Derivatives for voltage and ionic concentrations
dVdt = -(Id + Iion) ;

dNai = -(INa + INab + ICaL_Na + 3*INCX + 3*INaK)*p.Acap/(p.Vmyo*p.F) ;
dKi = -(IKr + IKs + IK1 + IKp + ICaL_K - 2*INaK)*p.Acap/(p.Vmyo*p.F) ;

dCai = Bi*( -Jserca*p.VNSR/p.Vmyo + Jrel*p.VJSR/p.Vmyo - ...
  (ICaL + ICaT + ICab + IpCa - 2*INCX)*p.Acap/(2*p.Vmyo*p.F) ) ;
dCaJSR = BJSR*(Jtr - Jrel) ;
dCaNSR = Jserca - Jtr*p.VJSR/p.VNSR ;

deriv = [dVdt;dCai;dCaNSR;dCaJSR;dNai;dKi; ...
  dmdt;dhdt;djdt;dddt;dfdt;dbdt;dgdt; ... 
  dxKrdt;dxs1dt;dxs2dt;dJreldt] ;

return


