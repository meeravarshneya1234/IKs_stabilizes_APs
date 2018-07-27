function deriv = dydt_Grandi(t,y,Id,p,c)

ydot = zeros(size(y));                                    

V = y(39);

% Nernst Potentials
ena_junc = (1/p.FoRT)*log(p.Nao/y(32));     % [mV]
ena_sl = (1/p.FoRT)*log(p.Nao/y(33));       % [mV]
% ena_junc = (1/FoRT)*log(Nao/7.65);     % [mV]
% ena_sl = (1/FoRT)*log(Nao/7.65);       % [mV]
ek = (1/p.FoRT)*log(p.Ko/y(35));	        % [mV]
eca_junc = (1/p.FoRT/2)*log(p.Cao/y(36));   % [mV]
eca_sl = (1/p.FoRT/2)*log(p.Cao/y(37));     % [mV]
ecl = (1/p.FoRT)*log(p.Cli/p.Clo);            % [mV]
% %% K current parameters
pNaK = 0.01833;      
eks = (1/p.FoRT)*log((p.Ko+pNaK*p.Nao)/(y(35)+pNaK*y(34)));

% %% Membrane Currents

I_Na_junc = p.Fjunc*c.GNa*y(1)^3*y(2)*y(3)*(V-ena_junc);
I_Na_sl = p.Fsl*c.GNa*y(1)^3*y(2)*y(3)*(V-ena_sl);
I_Na = I_Na_junc+I_Na_sl;

I_nabk_junc = p.Fjunc*c.GNaB*(V-ena_junc);
I_nabk_sl = p.Fsl*c.GNaB*(V-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

% I_nak: Na/K Pump Current
sigma = (exp(p.Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*V*p.FoRT)+0.0365*sigma*exp(-V*p.FoRT));
I_nak_junc = 1*p.Fjunc*c.IbarNaK*fnak*p.Ko /(1+(p.KmNaip/y(32))^4) /(p.Ko+p.KmKo);
I_nak_sl = 1*p.Fsl*c.IbarNaK*fnak*p.Ko /(1+(p.KmNaip/y(33))^4) /(p.Ko+p.KmKo);
I_nak = I_nak_junc+I_nak_sl;


% %% I_kr: Rapidly Activating K Current
gkr = c.GKr_*sqrt(p.Ko/5.4);
xrss = 1/(1+exp(-(V+10)/5));
rkr = 1/(1+exp((V+74)/24));
I_kr = gkr*y(12)*rkr*(V-ek);

gks_junc = c.GKs_ ;
gks_sl = c.GKs_ ;
I_ks_junc = p.Fjunc*gks_junc*y(13)^2*(V-eks);
I_ks_sl = p.Fsl*gks_sl*y(13)^2*(V-eks);                                                                                                                                   
I_ks = I_ks_junc+I_ks_sl;

%I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-V/5.98));
I_kp_junc = p.Fjunc*c.gkp*kp_kp*(V-ek);
I_kp_sl = p.Fsl*c.gkp*kp_kp*(V-ek);
I_kp = I_kp_junc+I_kp_sl;

I_tos = c.GtoSlow*y(8)*y(9)*(V-ek);    % [uA/uF]
I_tof = c.GtoFast*y(10)*y(11)*(V-ek);
I_to = I_tos + I_tof;

% %% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(V-ek-59.215)));
bki =(0.49124*exp(0.08032*(V+5.476-ek)) + exp(0.06175*(V-ek-594.31))) /(1 + exp(-0.5143*(V-ek+4.753)));
% Ak1=0.1/(1+exp(0.06*(y(39)-ek-200)));
% Bk1=(3*exp(0.0002*(y(39)-ek+100))+ exp(0.1*(y(39)-ek-10)))/(1+exp(-0.5*(y(39)-ek)));
% kiss=Ak1/(Ak1+Bk1);
% I_ki = GK1*sqrt(Ko/5.4)*kiss*(y(39)-ek);
kiss = aki/(aki+bki);
I_ki =c.GK1_*sqrt(p.Ko/5.4)*kiss*(V-ek);

% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = p.Fjunc*c.GClCa/(1+p.KdClCa/y(36))*(V-ecl);
I_ClCa_sl = p.Fsl*c.GClCa/(1+p.KdClCa/y(37))*(V-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = c.GClB*(V-ecl);

% % % L-type Ca current
fcaCaMSL= 0.1/(1+(0.01/y(37)));
fcaCaj= 0.1/(1+(0.01/y(36)));
fcaCaMSL=0;
fcaCaj= 0;
ibarca_j = c.PCa_*4*(V*p.Frdy*p.FoRT) * (0.341*y(36)*exp(2*V*p.FoRT)-0.341*p.Cao) /(exp(2*V*p.FoRT)-1);
ibarca_sl = c.PCa_*4*(V*p.Frdy*p.FoRT) * (0.341*y(37)*exp(2*V*p.FoRT)-0.341*p.Cao) /(exp(2*V*p.FoRT)-1);
ibark = p.pK*(V*p.Frdy*p.FoRT)*(0.75*y(35)*exp(V*p.FoRT)-0.75*p.Ko) /(exp(V*p.FoRT)-1);
ibarna_j = p.pNa*(V*p.Frdy*p.FoRT) *(0.75*y(32)*exp(V*p.FoRT)-0.75*p.Nao)  /(exp(V*p.FoRT)-1);
ibarna_sl = p.pNa*(V*p.Frdy*p.FoRT) *(0.75*y(33)*exp(V*p.FoRT)-0.75*p.Nao)  /(exp(V*p.FoRT)-1);
I_Ca_junc = (p.Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*p.Q10CaL^p.Qpow)*0.45*1;
I_Ca_sl = (p.Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*p.Q10CaL^p.Qpow)*0.45*1;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(p.Fjunc_CaL*(fcaCaj+(1-y(6)))+p.Fsl_CaL*(fcaCaMSL+(1-y(7))))*p.Q10CaL^p.Qpow)*0.45*1;
I_CaNa_junc = (p.Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*p.Q10CaL^p.Qpow)*0.45*1;
I_CaNa_sl = (p.Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*p.Q10CaL^p.Qpow)*.45*1;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(p.Kdact/y(36))^2);
Ka_sl = 1/(1+(p.Kdact/y(37))^2);
s1_junc = exp(p.nu*V*p.FoRT)*y(32)^3*p.Cao;
s1_sl = exp(p.nu*V*p.FoRT)*y(33)^3*p.Cao;
s2_junc = exp((p.nu-1)*V*p.FoRT)*p.Nao^3*y(36);
s3_junc = p.KmCai*p.Nao^3*(1+(y(32)/p.KmNai)^3) + p.KmNao^3*y(36)*(1+y(36)/p.KmCai)+p.KmCao*y(32)^3+y(32)^3*p.Cao+p.Nao^3*y(36);
s2_sl = exp((p.nu-1)*V*p.FoRT)*p.Nao^3*y(37);
s3_sl = p.KmCai*p.Nao^3*(1+(y(33)/p.KmNai)^3) + p.KmNao^3*y(37)*(1+y(37)/p.KmCai)+p.KmCao*y(33)^3+y(33)^3*p.Cao+p.Nao^3*y(37);

I_ncx_junc = p.Fjunc*c.IbarNCX*p.Q10NCX^p.Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+p.ksat*exp((p.nu-1)*V*p.FoRT));
I_ncx_sl = p.Fsl*c.IbarNCX*p.Q10NCX^p.Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+p.ksat*exp((p.nu-1)*V*p.FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = p.Fjunc*p.Q10SLCaP^p.Qpow*c.IbarSLCaP*y(36)^1.6/(p.KmPCa^1.6+y(36)^1.6);
I_pca_sl = p.Fsl*p.Q10SLCaP^p.Qpow*c.IbarSLCaP*y(37)^1.6/(p.KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

% I_cabk: Ca Background Current
I_cabk_junc = p.Fjunc*c.GCaB*(V-eca_junc);
I_cabk_sl = p.Fsl*c.GCaB*(V-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Update gating variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mss = 1 / ((1 + exp( -(56.86 + V) / 9.03 ))^2);
taum = (1)*(0.1292 * exp(-((V+45.79)/15.54)^2) + 0.06487 * exp(-((V-4.823)/51.12)^2));                 

ah = (V >= -40) * (0)... 
   + (V < -40) * (0.057 * exp( -(V + 80) / 6.8 )); 
bh = (V >= -40) * (0.77 / (0.13*(1 + exp( -(V + 10.66) / 11.1 )))) ...
   + (V < -40) * ((2.7 * exp( 0.079 * V) + 3.1*10^5 * exp(0.3485 * V))); 
hss = 1 / ((1 + exp( (V + 71.55)/7.43 ))^2);

ah = (V >= -40) * (0)... 
   + (V < -40) * (0.057 * exp( -(V + 80) / 6.8 )); 
bh = (V >= -40) * (0.77 / (0.13*(1 + exp( -(V + 10.66) / 11.1 )))) ...
   + (V < -40) * ((2.7 * exp( 0.079 * V) + 3.1*10^5 * exp(0.3485 * V))); 
tauh = 1 / ((ah + bh)); 

aj = (V >= -40) * (0) ...
    +(V < -40) * (((-2.5428 * 10^4*exp(0.2444*V) - 6.948*10^-6 * exp(-0.04391*V)) * (V + 37.78)) / ...
                     (1 + exp( 0.311 * (V + 79.23) )));
bj = (V >= -40) * ((0.6 * exp( 0.057 * V)) / (1 + exp( -0.1 * (V + 32) ))) ...
   + (V < -40) * ((0.02424 * exp( -0.01052 * V )) / (1 + exp( -0.1378 * (V + 40.14) ))); 
jss = 1 / ((1 + exp( (V + 71.55)/7.43 ))^2);   

aj = (V >= -40) * (0) ...
    +(V < -40) * (((-2.5428 * 10^4*exp(0.2444*V) - 6.948*10^-6 * exp(-0.04391*V)) * (V + 37.78)) / ...
                     (1 + exp( 0.311 * (V + 79.23) )));
bj = (V >= -40) * ((0.6 * exp( 0.057 * V)) / (1 + exp( -0.1 * (V + 32) ))) ...
   + (V < -40) * ((0.02424 * exp( -0.01052 * V )) / (1 + exp( -0.1378 * (V + 40.14) ))); 
tauj = 1 / ((aj + bj));

ydot(1) = (mss - y(1)) / taum;
ydot(2) = (hss - y(2)) / tauh;
ydot(3) = (jss - y(3)) / tauj;
    
% %% I_kr: Rapidly Activating K Current
xrss = 1/(1+exp(-(V+10)/5));
tauxr = (1)*(550/(1+exp((-22-V)/9))*6/(1+exp((V-(-11))/9))+230/(1+exp((V-(-40))/20)));
ydot(12) = (xrss-y(12))/tauxr;

xsss = 1 / (1+exp(-(V + 3.8)/14.25)); % fitting Fra
tauxs=(1)*(990.1/(1+exp(-(V+2.436)/14.12)));
ydot(13) = (xsss-y(13))/tauxs;


% %% I_to: Transient Outward K Current (slow and fast components)
% modified for human myocytes
xtoss = 1/(1+exp(-(V-19.0)/13));
ytoss = 1/(1+exp((V+19.5)/5));
% rtoss = 1/(1+exp((y(39)+33.5)/10));
tauxtos = (1)*(9/(1+exp((V+3.0)/15))+0.5);
tauytos = (1)*(800/(1+exp((V+60.0)/10))+30);
% taurtos = 2.8e3/(1+exp((y(39)+60.0)/10))+220; %Fei changed here!! time-dependent gating variable
ydot(8) = (xtoss-y(8))/tauxtos;
ydot(9) = (ytoss-y(9))/tauytos;

tauxtof = (1)*(8.5*exp(-((V+45)/50)^2)+0.5);
%tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5;
tauytof = (1)*(85*exp((-(V+40)^2/220))+7);
%tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;

fss = 1/(1+exp((V+35)/9))+0.6/(1+exp((50-V)/20));
tauf = (1)*(1/(0.0197*exp( -(0.0337*(V+14.5))^2 )+0.02));

dss = 1/(1+exp(-(V+5)/6.0));
taud = (1)*(dss*(1-exp(-(V+5)/6.0))/(0.035*(V+5)));

ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl

% %% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(p.ec50SR/y(31))^2.5);
koSRCa = p.koCa/kCaSR;
kiSRCa = p.kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (p.kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-p.kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-p.kom*y(15))-(kiSRCa*y(36)*y(15)-p.kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-p.kim*y(16))-(p.kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = c.ks*y(15)*(y(31)-y(36));          % [mM/ms]
% if t<2000
% J_SRCarel = c.ks*(y(31)-y(36));          % [mM/ms]
% end
J_serca = 1*p.Q10SRCaP^p.Qpow*c.Vmax_SRCaP*((y(38)/p.Kmf)^p.hillSRCaP-(y(31)/p.Kmr)^p.hillSRCaP)...
    /(1+(y(38)/p.Kmf)^p.hillSRCaP+(y(31)/p.Kmr)^p.hillSRCaP);
J_SRleak = c.Kleak_*(y(31)-y(36));           %   [mM/ms]


% %% Sodium and Calciurm Buffering
ydot(17) = p.kon_na*y(32)*(p.Bmax_Naj-y(17))-p.koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = p.kon_na*y(33)*(p.Bmax_Nasl-y(18))-p.koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = p.kon_tncl*y(38)*(p.Bmax_TnClow-y(19))-p.koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = p.kon_tnchca*y(38)*(p.Bmax_TnChigh-y(20)-y(21))-p.koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = p.kon_tnchmg*p.Mgi*(p.Bmax_TnChigh-y(20)-y(21))-p.koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = p.kon_cam*y(38)*(p.Bmax_CaM-y(22))-p.koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = p.kon_myoca*y(38)*(p.Bmax_myosin-y(23)-y(24))-p.koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = p.kon_myomg*p.Mgi*(p.Bmax_myosin-y(23)-y(24))-p.koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = p.kon_sr*y(38)*(p.Bmax_SR-y(25))-p.koff_sr*y(25);                    % SRB       [mM/ms]
J_CaB_cytosol = sum(ydot(19:25));

% Junctional and SL Ca Buffers
ydot(26) = p.kon_sll*y(36)*(p.Bmax_SLlowj-y(26))-p.koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = p.kon_sll*y(37)*(p.Bmax_SLlowsl-y(27))-p.koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = p.kon_slh*y(36)*(p.Bmax_SLhighj-y(28))-p.koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = p.kon_slh*y(37)*(p.Bmax_SLhighsl-y(29))-p.koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

% %% Ion concentrations
% SR Ca Concentrations
ydot(30) = p.kon_csqn*y(31)*(p.Bmax_Csqn-y(30))-p.koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*p.Vmyo/p.Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current
% ydot(31)=0;

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*p.Cmem/(p.Vjunc*p.Frdy)+p.J_na_juncsl/p.Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*p.Cmem/(p.Vsl*p.Frdy)+p.J_na_juncsl/p.Vsl*(y(32)-y(33))...
   +p.J_na_slmyo/p.Vsl*(y(34)-y(33))-ydot(18);
%FluxNaSL=ydot(33);
% ydot(32) = 0;
% ydot(33) = 0;
ydot(34) = p.J_na_slmyo/p.Vmyo*(y(33)-y(34));             % [mM/msec] 
% ydot(34)=0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp;     % [uA/uF]
% ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*p.Frdy);           % [mM/msec]
ydot(35) =0; % -I_K_tot*Cmem/(Vmyo*p.Frdy);

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*p.Cmem/(p.Vjunc*2*p.Frdy)+p.J_ca_juncsl/p.Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*p.Vsr/p.Vjunc+J_SRleak*p.Vmyo/p.Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*p.Cmem/(p.Vsl*2*p.Frdy)+p.J_ca_juncsl/p.Vsl*(y(36)-y(37))...
    + p.J_ca_slmyo/p.Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
% ydot(36)=0;
% ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+p.J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*p.Vsr/p.Vmyo-J_CaB_cytosol +p.J_ca_slmyo/p.Vmyo*(y(37)-y(38));
% ydot(38)=0;
%if (t<15000)
%    ydot(41) = 0;
%    ydot(42) = 0;
%else
%ydot(41) = -I_Na*Cmem/(p.Vmyo*p.Frdy);
%ydot(42) = -I_Ca_sl*Cmem/(Vjunc*2*p.Frdy)*Vjunc/p.Vmyo;
%end


% %% Membrane Potential
% %%
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
%ydot(39) = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);
ydot(39) = -(I_tot+Id);
deriv = ydot ;

return
