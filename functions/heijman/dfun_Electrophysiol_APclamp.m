% *************************************************************************
% *   Copyright (C) 2011 by Jordi Heijman, Paul G. A. Volders,            *
% *   Ronald L. Westra and Yoram Rudy                                     *
% *   Email: jordi.heijman@maastrichtuniversity.nl  / rudy@wustl.edu      *
% *   Web:   http://rudylab.wustl.edu                                     *
% *                                                                       *
% *   MATLAB Implementation of the model described in the article:        *
% *   "Local control of beta-adrenergic stimulation: Effects on           *
% *   ventricular myocyte electrophysiology and Ca2+ transient"           *
% *                                                                       *
% *   Journal of Molecular and Cellular Cardiology. 2011.                 *
% *                                                                       *
% *   Usage:                                                              *
% *   When called with the current state (y) and settings structure, this *
% *   function will calculate dy/dt for the electrophysiological          *
% *   components of the model.                                            *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   constants_Electrophysiol.m                                          *
% *                                                                       *
% *   This program is free software; you can redistribute it and/or modify*
% *   it under the terms of the GNU General Public License as published by*
% *   the Free Software Foundation; either version 2 of the License, or   *
% *   (at your option) any later version.                                 *
% *                                                                       *
% *   This program is distributed in the hope that it will be useful,     *
% *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
% *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
% *   GNU General Public License for more details.                        *
% *                                                                       *
% *   You should have received a copy of the GNU General Public License   *
% *   along with this program; if not, write to the                       *
% *   Free Software Foundation, Inc.,                                     *
% *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.           *
% *************************************************************************

function dy = dfun_Electrophysiol_APclamp(t, y, V,flags, settings)
    % The new Electrophysiological state vector looks as follows:
    %   1  Vm           2  Ca_i
    %   3  Ca_ss        4  Ca_ssCaL
    %   5  Ca_JSR       6  Ca_NSR
    %   7  Ca_trap      8  Na_i
    %   9  Na_ss       10  K_i
    %  11  Cl_i        12  Cl_ss
    %  13  Irel        14  H
    %  15  m           16  J
    %  17  mL          18  hL
    %  19  to_a        20  to_if
    %  21  to_is       22  xr
    %  23  AA
    %
    %  24 - 31 ICaL Markov states (NP)
    %  32 - 48 IKs Markov states (NP)
    %  49 - 56 ICaL Markov states (P)
    %  57 - 73 IKs Markov states (P)
    %  74  H (P)         75  m (P)
    %  76  J (P)         77  m (P, CAMKII)
    %  78  h (P, CAMKII) 79  J (P, CAMKII)
    %  80  IRel (P)
    %  81  to_if_CaMKII  82  to_is_CaMKII
    %  83  fPLBP_CaMKII  84  fRyRP_CaMKII
    %  85  fIToP_CaMKII  86  fINaP_CaMKII
    %  87  fIK1P_CaMKII  88  fICaLP_CaMKII
    y = [V;y];
    %y(1) = V;
    data = settings.dataElectrophysiol;
    dy = zeros(size(y,1), size(y,2));
            
    %Calcium Buffers
    y(2) = Buff_Ca_i(y(2), data, settings);
    if isfield(settings, 'Cai_clamp')
        if settings.Cai_clamp >= 0, y(2) = settings.Cai_clamp; end
    end
    y(3) = Buff_Ca_ss(y(3), data);
    y(4) = Buff_Ca_ss(y(4), data);
    y(5) = Buff_Ca_jsr(y(5), data);
    
    y(13) = (settings.IrelB) * y(13);
    
    %Determine stimulus current
	if settings.applyVoltageClamp == 1
        if isfield(settings, 'vcpFileData') && ~isempty(settings.vcpFileData)
            [val, minind] = min(abs(t - settings.vcpFileData(:,1)));            
            Istim = settings.vcRate * (settings.vcpFileData(minind,2) - y(1));
            y(1) = settings.vcpFileData(minind,2);
        else
            Istim = settings.vcRate * (settings.Vhold - y(1));
        end
	else
        Istim = settings.Istim * (t > 0 && t <= 0+settings.stimdur);
	end   
    
    Cactive = Ca_CaMKII(y, settings, data);    
    [dyt, INa] = Ion_HH_Na(y, Cactive, settings, data); dy = dy + dyt;
    [dyt, INaCa, INaCass] = Ion_HH_NaCa(y, settings, data); dy = dy + dyt;
    [dyt, INaK] = Ion_HH_NaK(y, settings, data); dy = dy + dyt;
    [dyt, INaL] = Ion_HH_NaL(y, Cactive, settings, data); dy = dy + dyt;
    [dyt, INab] = Ion_HH_Nab(y, settings, data); dy = dy + dyt;
    [dyt, IK1] = Ion_HH_K1(y, settings, data); dy = dy + dyt;
    [dyt, IKp] = Ion_HH_Kp(y, settings, data); dy = dy + dyt;
    [dyt, IKr] = Ion_HH_Kr(y, settings, data); dy = dy + dyt;
    [dyt, ITo1] = Ion_HH_To1(y, settings, data, Cactive); dy = dy + dyt;
    [dyt, CTNaCl,CTKCl,IClb,ITo2] = Ion_HH_Cl(y, settings, data); dy = dy + dyt;
    
    [dyt, IKs] = Ion_Markov_Ks(y, settings, data); dy = dy + dyt;
    [dyt, ICaL, ICaLNP, ICaLP] = Ion_Markov_CaL(y, Cactive, settings, data); dy = dy + dyt;
    if isfield(settings, 'applyJSRClamp') && settings.Vhold == 0
        ICaL = (1 ./ (1 + exp(-(t - 0.1) ./ 0.05))) .* -1.311 .* (1 ./ (1 + exp(-(t - 0.7966) ./ 0.7354))) .* (( (1 ./ 0.1803) ./ (1 + exp( (t - 0.7966) ./ 16.03)) ) + 0.1803);
    end
    
    % Calcium buffers and time-independent calcium currents
    [dyt, Ileak, Iup, IpCa, ICab, Itr, Idiff_ssCaL, Idiff_Ca, Idiff_Na, Idiff_Cl] = Ion_HH_Ca(y, Cactive, settings, data); dy = dy + dyt;
    
    [dyt, Ileak_JSR] = Ca_CICR(y, ICaL, settings, data); dy = dy + dyt;
    
    caiont = ICaL+ICab+IpCa-2*(INaCa+INaCass);
    naiont = INa+3*(INaCa+INaCass)+3*INaK+INaL+INab;
    kiont = IKr+IKs+IK1+IKp-2*INaK+ITo1+Istim;
    clont = IClb+ITo2;
    
    dy(1,:) =[];
    %dy(1) = -(naiont+kiont+caiont+clont);  
    dyt = UpdateConcentrations(y, ICaL, ICab, IpCa, INa, INaL, INaCa, INaCass, INaK, INab, IClb, ITo2, Iup, Ileak, Itr, Idiff_Ca, Idiff_ssCaL, Idiff_Na, Idiff_Cl, CTNaCl, CTKCl, kiont, Ileak_JSR, data, settings); 
    dyt(1,:)=[];   
    dy = dy + (dyt);
end

function [dyt, INa] = Ion_HH_Na(y, Cactive, settings, data)       
    dyt = zeros(size(y,1), size(y,2));        
    
    V = y(1);   
    Na_i = y(8);
    
    ENa =log(settings.Na_o./Na_i)/data.frt;       % Nernst potential of Na, mV    
    
    %Dual phosphorylated (both PKA and CaMK, kinetics are combined)
    %============================================
	params = settings.INaNPParams;
    gNa = params(1) * y(75) * y(75) * y(75) * y(78) * y(79);
	INaBothP = gNa * (V - ENa);
    
    %Non Phosphorylated
    %============================================
    params = settings.INaNPParams;
    H = y(14);
    m = y(15);
    J = y(16);
    
    gNa = params(1)*m*m*m*H*J;
    INaNP = gNa*(V-ENa);

    if V >= -40
        ah = 0;
        bh = 1 / (params(2)*(1+exp((V+params(3))/params(4))));
        aj = 0;
        bj = (params(5)*exp(params(6)*V)./(1+exp(params(7)*(V+params(8)))));
    else
        ah = params(9).*exp((params(10)+V)./params(11));
        bh = (params(12)*exp(params(13)*(V-params(29)))+params(14)*exp(params(15)*(V-params(29))));
        aj = (params(16)*exp(params(17)*V)+params(18)*exp(params(19)*V)).*(V+params(20))./(1+exp(params(21)*(V+params(22))));
        bj = (params(23)*exp(params(24)*V)./(1+exp(params(25)*(V+params(26)))));
    end
    
    am = 0.32*(V-params(27))/(1-exp(-0.1*(V-params(27))));
    bm = 0.08*exp(-(V+params(28))/11);
    
    dyt(14) = ( ah*(1-H)-bh*H );
    dyt(15) = ( am*(1-m)-bm*m );
    dyt(16) = ( aj*(1-J)-bj*J );
    
    % PKA Phosphorylated
    %============================================
    dVIn = -4.9;
    dVAc = -3.7;
    params = settings.INaNPParams;
    H = y(74);
    m = y(75);
    J = y(76);
    
    gNa = 1.25*params(1)*m*m*m*H*J;
    INaPKAP = gNa*(V-ENa);

    if V >= -40
        ah = 0;
        bh = 1 / (params(2)*(1+exp(((V-dVIn)+params(3))/params(4))));
        aj = 0;
        bj = (params(5)*exp(params(6)*(V-dVIn))./(1+exp(params(7)*((V-dVIn)+params(8)))));
    else
        ah = params(9).*exp((params(10)+(V-dVIn))./params(11));
        bh = (params(12)*exp(params(13)*(V-dVIn-params(29)))+params(14)*exp(params(15)*(V-dVIn-params(29))));
        aj = (params(16)*exp(params(17)*(V-dVIn))+params(18)*exp(params(19)*(V-dVIn))).*((V-dVIn)+params(20))./(1+exp(params(21)*((V-dVIn)+params(22))));
        bj = (params(23)*exp(params(24)*(V-dVIn))./(1+exp(params(25)*((V-dVIn)+params(26)))));
    end
    
    am = 0.32*(V-dVAc-params(27))/(1-exp(-0.1*(V-dVAc-params(27))));
    bm = 0.08*exp(-(V-dVAc+params(28))/11);
    
    dyt(74) = ( ah*(1-H)-bh*H );
    dyt(75) = ( am*(1-m)-bm*m );
    dyt(76) = ( aj*(1-J)-bj*J );
    
     % CaMKII Phosphorylated
    %============================================
    dVIn = -3.25;
    dVAc = 0.0;
    params = settings.INaNPParams;
    H = y(78);
    m = y(77);
    J = y(79);
    
    gNa = params(1)*m*m*m*H*J;
    INaCaMKP = gNa*(V-ENa);

    if V - dVIn >= -40
        ah = 0;
        bh = 1 / (params(2)*(1+exp(((V-dVIn)+params(3))/params(4))));
        aj = 0;
        bj = (params(5)*exp(params(6)*(V-dVIn))./(1+exp(params(7)*((V-dVIn)+params(8)))));
    else
        ah = params(9).*exp((params(10)+(V-dVIn))./params(11));
        bh = (params(12)*exp(params(13)*(V-dVIn-params(29)))+params(14)*exp(params(15)*(V-dVIn-params(29))));
        aj = (params(16)*exp(params(17)*(V-dVIn))+params(18)*exp(params(19)*(V-dVIn))).*((V-dVIn)+params(20))./(1+exp(params(21)*((V-dVIn)+params(22))));
        bj = (params(23)*exp(params(24)*(V-dVIn))./(1+exp(params(25)*((V-dVIn)+params(26)))));
    end
    
    am = 0.32*(V-dVAc-params(27))/(1-exp(-0.1*(V-dVAc-params(27))));
    bm = 0.08*exp(-(V-dVAc+params(28))/11);
    
    dyt(78) = ( ah*(1-H)-bh*H );
    dyt(77) = ( am*(1-m)-bm*m );
    dyt(79) = ( aj*(1-J)-bj*J );
    
    % Combine
	f_INa_CaMKP = y(86);
	f_INa_BothP = f_INa_CaMKP * settings.fINaP;
	f_INa_CaMK_Only = f_INa_CaMKP - f_INa_BothP;
	f_INa_PKA_Only = settings.fINaP - f_INa_BothP;
	f_INa_NP = 1 - f_INa_CaMK_Only - f_INa_PKA_Only - f_INa_BothP;

	INa = (settings.INaB) * (f_INa_NP * INaNP + f_INa_PKA_Only * INaPKAP + f_INa_CaMK_Only * INaCaMKP + f_INa_BothP * INaBothP);
end

function [dyt, INaL] = Ion_HH_NaL(y, Cactive, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    V = y(1);
    HL = y(18);
    mL = y(17);
    Na_i = y(8);
    
    ENa =log(settings.Na_o./Na_i)/data.frt;       % Nernst potential of Na, mV
    dGNaLCaMK = 0.0095 * y(86);
    INaL = (settings.INaLB)*(data.GNaL+dGNaLCaMK)*mL^3*HL*(V-ENa);

    amL=0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
    bmL=0.08*exp(-V/11.0);

    hLss=1/(1+exp((V+91)/6.1));
    
    dyt(17) = amL*(1-mL)-bmL*mL;
    dyt(18) = (hLss-HL)/600;
end

function [dyt, INaCa, INaCass] = Ion_HH_NaCa(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    V = y(1);
    Na_i = y(8);
    Ca_i = y(2);
    Na_ss = y(9);
    Ca_ss = y(3);
    
    num = 0.8 * data.NCXmax * (Na_i^3 * settings.Ca_o * exp(data.eta * V * data.frt) - settings.Na_o^3 * Ca_i * exp((data.eta-1) * V * data.frt));
    denom1 = 1 + (data.KmCa / Ca_i).^2;
    denom2 = 1 + data.ksat * exp((data.eta - 1) * V * data.frt);
    denom3 = data.KmCao * Na_i^3 + data.KmNao^3 * Ca_i + data.KmNai^3 * settings.Ca_o * (1 + Ca_i/data.KmCai);
    denom4 = data.KmCai * settings.Na_o^3 * (1 + (Na_i/data.KmNai)^3) + Na_i^3 * settings.Ca_o + settings.Na_o^3 * Ca_i;
    
    INaCa = (settings.INaCaB) * ( (num) / (denom1 * denom2 * (denom3 + denom4)) );
    
    num = 0.2 * data.NCXmax * (Na_ss^3 * settings.Ca_o * exp(data.eta * V * data.frt) - settings.Na_o^3 * Ca_ss * exp((data.eta-1) * V * data.frt));
    denom1 = 1 + (data.KmCa / Ca_ss).^2;
    denom2 = 1 + data.ksat * exp((data.eta - 1) * V * data.frt);
    denom3 = data.KmCao * Na_ss^3 + data.KmNao^3 * Ca_ss + data.KmNai^3 * settings.Ca_o * (1 + Ca_ss/data.KmCai);
    denom4 = data.KmCai * settings.Na_o^3 * (1 + (Na_ss/data.KmNai)^3) + Na_ss^3 * settings.Ca_o + settings.Na_o^3 * Ca_ss;
    
    INaCass = (settings.INaCaB) * ((num) / (denom1 * denom2 * (denom3 + denom4)));    
end

function [dyt, INaK] = Ion_HH_NaK(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    V = y(1);
    Na_i = y(8);

    V_half = -92;
    phi = 1 * (V - V_half) * data.frt;
    fv = 1 / (1 + exp(-phi));
    
    Pk = settings.K_o / (settings.K_o + data.kmko);       
    
    INaK_NP = data.ibarnak .* fv .* Pk .* ((Na_i ./ (Na_i + data.kmnai_NP)).^3);
    INaK_P = data.ibarnak .* fv .* Pk .* ((Na_i ./ (Na_i + data.kmnai_P)).^3);            
    INaK = (settings.INaKB) .* (1 - settings.fINaKP) .* INaK_NP + settings.fINaKP .* INaK_P;
end

function [dyt, INab] = Ion_HH_Nab(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    V = y(1);
    Na_i = y(8);
    
    phi = V * data.frt;
    PNab = 0.32E-8;

    INab = (settings.INabB) * data.F * PNab * phi * ( (Na_i * exp(phi) - settings.Na_o) / (exp(phi) - 1) );
end

function [dyt, IK1] = Ion_HH_K1(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    V = y(1);
    K_i = y(10);
    fIK1P_CaMKII = y(87);
    
    GK1 = data.GK1max * sqrt(settings.K_o/5.4);
    EK1 = log(settings.K_o/K_i)/data.frt;

    ak1 = 1.02./(1+exp(0.2385*(V-EK1-59.215)));
    bk1 = (0.49124*exp(0.08032*(V-EK1+5.476))+exp(0.06175*(V-EK1-594.31)))/(1+exp(-0.5143*(V-EK1+4.753)));

    IK1_NP = GK1*(ak1./(ak1+bk1))*(V-EK1);
    IK1_CaMKII = 1.2 * GK1*(ak1./(ak1+bk1))*(V-EK1);
    
    IK1 = (settings.IK1B) .* ((1 - fIK1P_CaMKII) .* IK1_NP + fIK1P_CaMKII * IK1_CaMKII);
end

function [dyt, IKp] = Ion_HH_Kp(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    V = y(1);
    K_i = y(10);    
    
    EK1 = log(settings.K_o/K_i)/data.frt;
    IKp_NP = data.GKpmax*(V-EK1)/(1+exp((15.0-V)/17));
    IKp_P = 3.62 * data.GKpmax * (V - EK1) / (1 + exp(36.0 - V)/17);
    
    IKp = (settings.IKpB) * ((1 - settings.fIKurP) * IKp_NP + settings.fIKurP * IKp_P);
end

function [dyt, IKr] = Ion_HH_Kr(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    V = y(1);
    K_i = y(10);
    xr = y(22);
    
    gkr = data.GKrmax*sqrt(settings.K_o/5.4);
    ekr = log(settings.K_o/K_i)/data.frt;
    r = 1./(1+exp((V+10)/15.4));

    IKr = (settings.IKrB) *gkr*xr*r*(V-ekr);
    xrss = 1./(1+exp(-(V+10.085)/4.25));
    tauxr = 1./(6e-4*(V-1.7384)./(1-exp(-0.136*(V-1.7384)))+3e-4*(V+38.3608)./(exp(0.1522*(V+38.3608))-1));
    
    dyt(22) = (xrss-xr)/tauxr;
end

function [dyt, ITo1] = Ion_HH_To1(y, settings, data, Cactive)
    dyt = zeros(size(y,1), size(y,2));

    V = y (1);
    K_i = y(10);
    i_f = y(20);
    i_s = y(21);
    a = y(19);
    i_f_CaMKII = y(81);
    i_s_CaMKII = y(82);
    
    EK = log(settings.K_o/K_i)/data.frt;
    Rto1 = exp(V / 550);
    
    alpha_a = 1 / (1.2089 * (1 + exp( (V - 18.4099) / -29.3814 )));
    beta_a = 3.5 / (1 + exp( (V + 100) / 29.3814 ));    
    tau_a = 1 / (alpha_a + beta_a);
    a_ss = 1 / (1 + exp( (V + 9.437) / -7.133 ));
    dyt(19) = (a_ss - a) / tau_a; 
        
    alpha_if = 0.02144 / (1 + exp( (V + 58) / 5 ));
    beta_if = 1.0 / (0.5 * 9.7953 * (1 + exp( (V + 19) / -9 )));
    
    alpha_is = 0.56034 / (250 * (1 + exp( (V + 60) / 5 )));
    beta_is = beta_if;

    dyt(20) = alpha_if * (1 - i_f) - beta_if * i_f;
    dyt(21) = alpha_is * (1 - i_s) - beta_is * i_s;    
    
    ITo1_NP = data.gitodv * a^3 * (0.7356 * i_f + 0.2644 * i_s) * Rto1 * (V - EK);

    alpha_if = 0.04796 / (1 + exp( (V + 58) / 5 ));
    beta_if = 1.0 / (0.5 * 9.7953 * (1 + exp( (V + 19) / -9 )));
    
    alpha_is = 2.4600 / (250 * (1 + exp( (V + 60) / 5 )));
    beta_is = beta_if;

    dyt(81) = alpha_if * (1 - i_f_CaMKII) - beta_if * i_f_CaMKII;
    dyt(82) = alpha_is * (1 - i_s_CaMKII) - beta_is * i_s_CaMKII;    
    
    ITo1_CaMKII = data.gitodv * a^3 * (0.7356 * i_f_CaMKII + 0.2644 * i_s_CaMKII) * Rto1 * (V - EK);

    fITo_CaMKII_P = y(85);
    if isfield(settings, 'fIToP_CaMKII'), fITo_CaMKII_P = settings.fIToP_CaMKII; end    
    ITo1 = (settings.ITo1B) .* ((1 - fITo_CaMKII_P) .* ITo1_NP + fITo_CaMKII_P .* ITo1_CaMKII);
end

function [dyt, CTNaCl,CTKCl,IClb,Ito2] = Ion_HH_Cl(y, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    V = y(1);
    K_i = y(10);
    Na_i = y(8);
    AA = y(23);
    Cl_i = y(11);    
    Irel = (settings.IrelB) * ((1 - settings.fRyRP) * y(13) + settings.fRyRP * y(80));
    
    ENa=log(settings.Na_o/Na_i)/data.frt;
    EK=log(settings.K_o/K_i)/data.frt;
    ECl=-log(settings.Cl_o/Cl_i)/data.frt;
    
    CTKCl=(1 - settings.CTKClB) * data.CTKClmax*(EK-ECl)/((EK-ECl)+87.8251);
    CTNaCl=(1 - settings.CTNaClB) * data.CTNaClmax * ((ENa-ECl)^4.0) / ((ENa-ECl)^4.0 + 87.8251^4.0);

    kCaIto2 = 0.4;
    KCaIto2 = 1 - (1 / (1 + (Irel / kCaIto2)^2));
    
    tau_AA = 8;
    alpha_AA = 0.025 / (1 + exp( (V + 58) / 5 ));
    beta_AA = 1 / (5 * (1 + exp( (V + 19) / -9 )));
    AA_ss = alpha_AA / (alpha_AA + beta_AA);
        
    Ito2_max=data.PCl*V*data.F*data.frt*(Cl_i-settings.Cl_o*exp(V*data.frt))/(1-exp(V*data.frt));
    Ito2 = (settings.ITo2B) * Ito2_max * AA * KCaIto2;
    
    IClb=(settings.IClB) * data.GClb*(V-ECl);
    
    dyt(23) = (AA_ss-AA) / tau_AA;
end

function [dyt, IKs] = Ion_Markov_Ks(Y, settings, data)
    dyt = zeros(size(Y,1), size(Y,2));    

    v = Y(1);
    Ca_i = Y(2);
    Na_i = Y(8);
	K_i = Y(10);    

    EKs = log((settings.K_o+data.prnak*settings.Na_o)/(K_i+data.prnak*Na_i))/data.frt;
    
	fKsP = settings.fIKsP;	
    
	%Non phosphorylated channel rates
	startindex = 1;
    params = settings.IKsNPParams;
    
    alphaNP = params(1) ./ (1 + exp((-(v - params(2))./params(3)) * data.F/(data.R*data.T)));
    betaNP = params(4) ./ (1 + exp(((v - params(5))./params(6)) * data.F/(data.R*data.T)));
    gammaNP = params(7) ./ (1 + exp((-(v - params(8))./params(9)) * data.F/(data.R*data.T)));
    deltaNP = params(10)*exp((v*data.F/(data.R*data.T))*params(11));
    etaNP = (params(12)-params(13))/(1+exp((v-params(14))/params(15) * data.F/(data.R*data.T))) + params(13);
    thetaNP = params(16);
    omegaNP = params(17)*exp((v*data.F/(data.R*data.T))*params(18));
    psiNP = params(19)*exp((v*data.F/(data.R*data.T))*params(20));
    muNP = params(21);

  	%Phosphorylated channel rates
    params = settings.IKsPParams;
    alphaP = params(1) ./ (1 + exp((-(v - params(2))./params(3)) * data.F/(data.R*data.T)));
    betaP = params(4) ./ (1 + exp(((v - params(5))./params(6)) * data.F/(data.R*data.T)));
    gammaP = params(7) ./ (1 + exp((-(v - params(8))./params(9)) * data.F/(data.R*data.T)));
    deltaP = params(10)*exp((v*data.F/(data.R*data.T))*params(11));
    etaP = (params(12)-params(13))/(1+exp((v-params(14))/params(15) * data.F/(data.R*data.T))) + params(13);
    thetaP = params(16);
    omegaP = params(17)*exp((v*data.F/(data.R*data.T))*params(18));
    psiP = params(19)*exp((v*data.F/(data.R*data.T))*params(20));
    muP = params(21);
    
    %Calculate all state derivatives for Non-phosphorylated states
    states = Y(32:48);
    if abs(sum(states)-1) > 1E-3, disp(['IKs NP State Error: Non unity sum of states! (sum = ', num2str(sum(states)), ')']); end
    dstates = zeros(17,1);
	O = 15; C = 0;
	dstates(O+2) = psiNP*states(O+1)															-omegaNP*states(O+2);
	dstates(O+1) = thetaNP*states(C+15)+omegaNP*states(O+2)										-(etaNP+psiNP)*states(O+1);
	dstates(C+15)=									   gammaNP*states(C+14)+etaNP*states(O+1)		-(4*deltaNP+thetaNP)*states(C+15);

	dstates(C+14)=   alphaNP*states(C+13)+				 2*gammaNP*states(C+12)+4*deltaNP*states(C+15)-(betaNP+3*deltaNP+gammaNP)*states(C+14);
	dstates(C+13)=					  betaNP*states(C+14)+  gammaNP*states(C+11)					-(alphaNP+3*deltaNP)*states(C+13);
	dstates(C+12)=   alphaNP*states(C+11)+				 3*gammaNP*states(C+9)+ 3*deltaNP*states(C+14)-(2*betaNP+2*deltaNP+2*gammaNP)*states(C+12);
	dstates(C+11)= 2*alphaNP*states(C+10)+2*betaNP*states(C+12)+2*gammaNP*states(C+8)+ 3*deltaNP*states(C+13)-(betaNP+alphaNP+gammaNP+2*deltaNP)*states(C+11);
	dstates(C+10)=					  betaNP*states(C+11)+  gammaNP*states(C+7)					-(2*alphaNP+2*deltaNP)*states(C+10);
	dstates(C+9) =   alphaNP*states(C+8)+					 4*gammaNP*states(C+5)+ 2*deltaNP*states(C+12)-(3*betaNP+deltaNP+3*gammaNP)*states(C+9);
	dstates(C+8) = 2*alphaNP*states(C+7)+	3*betaNP*states(C+9)+ 3*gammaNP*states(C+4)+ 2*deltaNP*states(C+11)- (2*betaNP+alphaNP+2*gammaNP+deltaNP)*states(C+8);
	dstates(C+7) = 3*alphaNP*states(C+6)+	2*betaNP*states(C+8)+ 2*gammaNP*states(C+3)+ 2*deltaNP*states(C+10)- (betaNP+2*alphaNP+deltaNP+gammaNP)*states(C+7);
	dstates(C+6) =					  betaNP*states(C+7)+   gammaNP*states(C+2)					- (3*alphaNP+deltaNP)*states(C+6);
	dstates(C+5) =   alphaNP*states(C+4)+										 deltaNP*states(C+9)	- (4*betaNP+4*gammaNP)*states(C+5);
	dstates(C+4) = 2*alphaNP*states(C+3)+	4*betaNP*states(C+5)+					 deltaNP*states(C+8)	- (3*betaNP+alphaNP+3*gammaNP)*states(C+4);
	dstates(C+3) = 3*alphaNP*states(C+2)+	3*betaNP*states(C+4)+					 deltaNP*states(C+7)	- (2*betaNP+2*alphaNP+2*gammaNP)*states(C+3);
	dstates(C+2) = 4*alphaNP*states(C+1)+	2*betaNP*states(C+3)+					 deltaNP*states(C+6) - (betaNP+3*alphaNP+gammaNP)*states(C+2);
	dstates(C+1) =					  betaNP*states(C+2)										- 4*alphaNP*states(C+1); 
   			
	IKsNP = muNP*(1+0.6/(1+(((3.8E-05)/Ca_i)^1.4)))*(states(O+1)+states(O+2))*(v-EKs);
    dyt(32:48) = dstates;
    if abs(sum(dstates)) > 1E-5, disp(['IKs NP Rates Error: Non zero sum of derivatives! (sum = ', num2str(abs(sum(dstates))), ')']); end
    
	%Calculate all state derivatives for PhosphorYlated states
    states = Y(57:73);
    if abs(sum(states)-1) > 1E-3, disp(['IKs P State Error: Non unity sum of states! (sum = ', num2str(sum(states)), ')']); end
    dstates = zeros(17,1);
	O = 15; C = 0;
	dstates(O+2) = psiP*states(O+1)															-omegaP*states(O+2);
	dstates(O+1) = thetaP*states(C+15)+omegaP*states(O+2)										-(etaP+psiP)*states(O+1);
	dstates(C+15)=									   gammaP*states(C+14)+etaP*states(O+1)		-(4*deltaP+thetaP)*states(C+15);

	dstates(C+14)=   alphaP*states(C+13)+				 2*gammaP*states(C+12)+4*deltaP*states(C+15)-(betaP+3*deltaP+gammaP)*states(C+14);
	dstates(C+13)=					  betaP*states(C+14)+  gammaP*states(C+11)					-(alphaP+3*deltaP)*states(C+13);
	dstates(C+12)=   alphaP*states(C+11)+				 3*gammaP*states(C+9)+ 3*deltaP*states(C+14)-(2*betaP+2*deltaP+2*gammaP)*states(C+12);
	dstates(C+11)= 2*alphaP*states(C+10)+2*betaP*states(C+12)+2*gammaP*states(C+8)+ 3*deltaP*states(C+13)-(betaP+alphaP+gammaP+2*deltaP)*states(C+11);
	dstates(C+10)=					  betaP*states(C+11)+  gammaP*states(C+7)					-(2*alphaP+2*deltaP)*states(C+10);
	dstates(C+9) =   alphaP*states(C+8)+					 4*gammaP*states(C+5)+ 2*deltaP*states(C+12)-(3*betaP+deltaP+3*gammaP)*states(C+9);
	dstates(C+8) = 2*alphaP*states(C+7)+	3*betaP*states(C+9)+ 3*gammaP*states(C+4)+ 2*deltaP*states(C+11)- (2*betaP+alphaP+2*gammaP+deltaP)*states(C+8);
	dstates(C+7) = 3*alphaP*states(C+6)+	2*betaP*states(C+8)+ 2*gammaP*states(C+3)+ 2*deltaP*states(C+10)- (betaP+2*alphaP+deltaP+gammaP)*states(C+7);
	dstates(C+6) =					  betaP*states(C+7)+   gammaP*states(C+2)					- (3*alphaP+deltaP)*states(C+6);
	dstates(C+5) =   alphaP*states(C+4)+										 deltaP*states(C+9)	- (4*betaP+4*gammaP)*states(C+5);
	dstates(C+4) = 2*alphaP*states(C+3)+	4*betaP*states(C+5)+					 deltaP*states(C+8)	- (3*betaP+alphaP+3*gammaP)*states(C+4);
	dstates(C+3) = 3*alphaP*states(C+2)+	3*betaP*states(C+4)+					 deltaP*states(C+7)	- (2*betaP+2*alphaP+2*gammaP)*states(C+3);
	dstates(C+2) = 4*alphaP*states(C+1)+	2*betaP*states(C+3)+					 deltaP*states(C+6) - (betaP+3*alphaP+gammaP)*states(C+2);
	dstates(C+1) =					  betaP*states(C+2)										- 4*alphaP*states(C+1); 
   			
	IKsP = muP*(1+0.6/(1+(((3.8E-05)/Ca_i)^1.4)))*(states(O+1)+states(O+2))*(v-EKs);
    dyt(57:73) = dstates;
    if abs(sum(dstates)) > 1E-5, disp(['IKs P Rates Error: Non zero sum of derivatives! (sum = ', num2str(abs(sum(dstates))), ')']); end
    
	IKs = (settings.IKsB) * (fKsP * IKsP + (1- fKsP) * IKsNP);
end

function [dyt, ICaL, ICaLNP, ICaLP] = Ion_Markov_CaL(Y, Cactive, settings, data)
    dyt = zeros(size(Y,1), size(Y,2));
        
    V = Y(1);           
    Ca_ssCaL = Y(3);
    if settings.observeICaLSS == 1, Ca_ssCaL = Y(4); end
        
    %% Nonphosphorylated parameters and rates
    params = settings.LCCNPParams;
    
    ACT_tau = params(1) + (params(2) * exp(params(3) * (V + params(4)))) / (1 + exp(params(5) * (V + params(4))));
    ACT_ss = 1 / ((1 + exp( - (V - params(6)) / (params(7)))) * (1 + exp( - (V + 25) / 5)));
        
    f_ICaL_P_CaMK = Y(88);
    if isfield(settings, 'fICaLP_CaMKII') && settings.fICaLP_CaMKII >= 0, f_ICaL_P_CaMK = settings.fICaLP_CaMKII; end
    
    I_Valpha = 1 / (((1 - f_ICaL_P_CaMK) * params(8) + f_ICaL_P_CaMK * 0.5 * params(8)) * (1 + exp( ( V + params(9) ) / params(10) )));
    I_Vbeta = 1 / (params(11) * (1 + exp( -(V + params(12)) / params(13) )));
    I_Vss = (1 / (1 + params(16))) * (params(16) + (1 / (1 + exp( (V + params(14)) / params(15) ))));
    I_Vtau = 1 / (I_Valpha + I_Vbeta);
    
    dtau_CaMKII = 5 * f_ICaL_P_CaMK;
    Is_Valpha = I_Valpha;
    Is_Vbeta_ss = params(17) - ( (params(18) - dtau_CaMKII) / (1 + (params(23) / Ca_ssCaL)^4) + params(24) / (1 + (params(25) / Ca_ssCaL) ^ 10) );
    Is_Vbeta = 1 / ((Is_Vbeta_ss) * (1 + exp( -(V + params(12)) / params(13) )));
    Is_Vss = (1 / (1 + params(19))) * (params(19) + (1 / (1 + exp( (V + params(14)) / params(15) ))));
    Is_Vtau = 1 / (Is_Valpha + Is_Vbeta);
    
    alpha = ACT_ss / ACT_tau;
    beta = (1 - ACT_ss) / ACT_tau;
    x = I_Vss / I_Vtau;
    y = (1 - I_Vss) / I_Vtau;
    xs = Is_Vss / Is_Vtau;
    ys = (1 - Is_Vss) / Is_Vtau;
    delta = params(20) / (1 + (params(23) / Ca_ssCaL)^4);
    theta = 1;
    theta_I = params(21);
    
    if abs(y) < 1E-12 || abs(xs) < 1E-12
        delta_I = theta_I * ( (x * ys * delta) / ((y + 0.0001) * (xs + 0.0001) * theta) );
    else
        delta_I = theta_I * ( (x * ys * delta) / (y * xs * theta) );
    end

	PCa = params(22) * (1 + 0.4*f_ICaL_P_CaMK);

    %% Nonphosphorylated states
	dyt(24) = beta * Y(25) + theta * Y(26) + x * Y(28) - (alpha + delta + y) * Y(24);
	dyt(25) = alpha * Y(24) + theta * Y(27) + x * Y(29) - (beta + delta + y) * Y(25);
	dyt(26) = delta * Y(24) + beta * Y(27) + xs * Y(30) - (theta + alpha + ys) * Y(26);
	dyt(27) = delta * Y(25) + alpha * Y(26) + xs * Y(31) - (theta + beta + ys) * Y(27);
	dyt(28) = y * Y(24) + beta * Y(29) + theta_I * Y(30) - (x + alpha + delta_I) * Y(28);
	dyt(29) = y * Y(25) + alpha * Y(28) + theta_I * Y(31) - (x + beta + delta_I) * Y(29);
	dyt(30) = ys * Y(26) + delta_I * Y(28) + beta * Y(31) - (xs + theta_I + alpha) * Y(30);
	dyt(31) = ys * Y(27) + delta_I * Y(29) + alpha * Y(30) - (xs + theta_I + beta) * Y(31);

    gamma_cai = 1;
    gamma_cao = 0.341;
    ICaLBar = PCa * 4 * V * data.F * data.frt * ((gamma_cai * Ca_ssCaL * exp(2*V*data.frt) - gamma_cao * settings.Ca_o) / (exp(2*V*data.frt) - 1));
    
    ICaLNP = ICaLBar * (Y(25) + Y(27));
    
    %% Phosphorylated parameters and rates
    params = settings.LCCPParams;   
    
    ACT_tau = params(1) + (params(2) * exp(params(3) * (V + params(4)))) / (1 + exp(params(5) * (V + params(4))));
    ACT_ss = 1 / ((1 + exp( - (V - params(6)) / (params(7)))) * (1 + exp( - (V + 25) / 5)));

    I_Valpha = 1 / (((1 - f_ICaL_P_CaMK) * params(8) + f_ICaL_P_CaMK * 0.5 * params(8)) * (1 + exp( ( V + params(9) ) / params(10) )));
    I_Vbeta = 1 / (params(11) * (1 + exp( -(V + params(12)) / params(13) )));
    I_Vss = (1 / (1 + params(16))) * (params(16) + (1 / (1 + exp( (V + params(14)) / params(15) ))));
    I_Vtau = 1 / (I_Valpha + I_Vbeta);
    
    dtau_CaMKII = 0.1 * f_ICaL_P_CaMK;
    Is_Valpha = I_Valpha;
    Is_Vbeta_ss = params(17) - ( (params(18) - dtau_CaMKII) / (1 + (params(23) / Ca_ssCaL)^4) + params(24) / (1 + (params(25) / Ca_ssCaL) ^ 10) );
    Is_Vbeta = 1 / ((Is_Vbeta_ss) * (1 + exp( -(V + params(12)) / params(13) )));
    Is_Vss = (1 / (1 + params(19))) * (params(19) + (1 / (1 + exp( (V + params(14)) / params(15) ))));
    Is_Vtau = 1 / (Is_Valpha + Is_Vbeta);
    
    alpha = ACT_ss / ACT_tau;
    beta = (1 - ACT_ss) / ACT_tau;
    x = I_Vss / I_Vtau;
    y = (1 - I_Vss) / I_Vtau;
    xs = Is_Vss / Is_Vtau;
    ys = (1 - Is_Vss) / Is_Vtau;
    delta = params(20) / (1 + (params(23) / Ca_ssCaL)^4);
    theta = 1;
    theta_I = params(21);
    
    if abs(y) < 1E-12 || abs(xs) < 1E-12
        delta_I = theta_I * ( (x * ys * delta) / ((y + 0.0001) * (xs + 0.0001) * theta) );
    else
        delta_I = theta_I * ( (x * ys * delta) / (y * xs * theta) );
    end

	PCa = params(22) * (1 + 0.1*f_ICaL_P_CaMK);
      
	dyt(49) = beta * Y(50) + theta * Y(51) + x * Y(53) - (alpha + delta + y) * Y(49);
	dyt(50) = alpha * Y(49) + theta * Y(52) + x * Y(54) - (beta + delta + y) * Y(50);
	dyt(51) = delta * Y(49) + beta * Y(52) + xs * Y(55) - (theta + alpha + ys) * Y(51);
	dyt(52) = delta * Y(50) + alpha * Y(51) + xs * Y(56) - (theta + beta + ys) * Y(52);
	dyt(53) = y * Y(49) + beta * Y(54) + theta_I * Y(55) - (x + alpha + delta_I) * Y(53);
	dyt(54) = y * Y(50) + alpha * Y(53) + theta_I * Y(56) - (x + beta + delta_I) * Y(54);
	dyt(55) = ys * Y(51) + delta_I * Y(53) + beta * Y(56) - (xs + theta_I + alpha) * Y(55);
	dyt(56) = ys * Y(52) + delta_I * Y(54) + alpha * Y(55) - (xs + theta_I + beta) * Y(56);

    gamma_cai = 1;
    gamma_cao = 0.341;
    ICaLBar = PCa * 4 * V * data.F * data.frt * ((gamma_cai * Ca_ssCaL * exp(2*V*data.frt) - gamma_cao * settings.Ca_o) / (exp(2*V*data.frt) - 1));
    
    ICaLP = ICaLBar * (Y(50) + Y(52));
    
    ICaL = (settings.ICaLB) * ((1 - settings.fICaLP) * ICaLNP + settings.fICaLP * ICaLP);
end

function [dyt, Ileak, Iup, IpCa, ICab, Itr, Idiff_ssCaL, Idiff_Ca, Idiff_Na, Idiff_Cl] = Ion_HH_Ca(y, Cactive, settings, data)
    dyt = zeros(size(y,1), size(y,2));
    
    v = y(1); 
    Ca_i = y(2);
    cass = y(3);
    Ca_ssCaL = y(4);
    jsr = y(5);
    nsr = y(6);
    Na_i = y(8);
    Na_ss = y(9);
    Cl_i = y(11);
    Cl_ss = y(12);
    
    IpCa = (settings.IpCaB)*data.ibarpca*Ca_i/(data.kmpca+Ca_i);	 % sarcolema pump

    ICab = (settings.ICabB) * (1.995E-7)*4*v*data.F*data.frt*(Ca_i*exp(2*v*data.frt)-0.341*settings.Ca_o)/(exp(2*v*data.frt)-1);
    Itr = (1 - settings.ItrB) * (nsr-jsr)/ settings.TauTr;
    
    fPLBP = settings.fPLBP;	
    fPLBBoth = fPLBP * y(83);
	fPLBNP = 1 - fPLBP - y(83) + fPLBBoth;
    fPLBCAMK = y(83) - fPLBBoth;
    fPLBPKA = fPLBP - fPLBBoth;        
    
    iupmaxCAMK = 3.25 * data.iupmax;
    frac_SERCA2a_P = 1 / (1 + (0.03 / Cactive) ^ 2);    
    Iup = (settings.IupB) * ((1 - frac_SERCA2a_P) * data.iupmax + frac_SERCA2a_P * iupmaxCAMK) * (Ca_i / (Ca_i + (fPLBNP * data.Kmup + fPLBCAMK * (data.Kmup - data.dKmPLBCAMK) + fPLBPKA * (data.Kmup - data.dKmPLBPKA) + fPLBBoth * (data.Kmup - data.dKmPLBBoth))));
    Ileak = (settings.IleakB) .* (((1 - frac_SERCA2a_P) * data.iupmax + frac_SERCA2a_P * iupmaxCAMK) / data.nsrmax) * nsr;

    tau_diff_ss = 0.2;
    tau_diff_ssCaL = 2.0;
    if settings.observeICaLSS ~= 1, tau_diff_ssCaL = 0.02; end
        
    Idiff_ssCaL = (cass - Ca_ssCaL) / tau_diff_ssCaL;
    Idiff_Ca = (1 - settings.IdiffB) * (cass - Ca_i) / tau_diff_ss;
    Idiff_Na = (1 - settings.IdiffB) * (Na_ss - Na_i) / tau_diff_ss;
    Idiff_Cl = (1 - settings.IdiffB) * (Cl_ss - Cl_i) / tau_diff_ss;
end

function dyt = UpdateConcentrations(y, ICaL, ICab, IpCa, INa, INaL, INaCa, INaCass, INaK, INab, IClb, Ito2, Iup, Ileak, Itr, Idiff_Ca, Idiff_ssCaL, Idiff_Na, Idiff_Cl, CTNaCl, CTKCl, IKtot, Ileak_JSR, data, settings)
    dyt = zeros(size(y,1), size(y,2));
    
    Irel = (settings.IrelB) * ((1 - settings.fRyRP) * y(13) + settings.fRyRP * y(80) + Ileak_JSR);
    CaJSR = y(5); Cass = y(3);    
    Ctrap = y(7);
    
    Cactive = Ca_CaMKII(y, settings, data);
    
    dCa_i = -((ICab+IpCa-2*INaCa)*data.AF/(data.vmyo*2)  +  (Iup - Ileak)*data.vnsr/data.vmyo  -  Idiff_Ca*data.vss/data.vmyo);
    dCa_ss = -( -2*INaCass*data.AF/(data.vss*2) - Irel*data.vjsr/data.vss + Idiff_Ca + Idiff_ssCaL);
    dCa_ssCaL = -(ICaL*data.AF/(data.vssCaL*2) - Idiff_ssCaL*data.vss/data.vssCaL);

    dNa_i = -((3*INaCa + 3*INaK + INa + INaL + INab)*data.AF/data.vmyo - CTNaCl - Idiff_Na*data.vss/data.vmyo);
    dNa_ss = -(3*INaCass*data.AF/data.vss + Idiff_Na);
    
    dCl_i = -(IClb*data.AF/(-1*data.vmyo) - CTNaCl - CTKCl - Idiff_Cl*data.vss/data.vmyo);
    dCl_ss = -(Ito2*data.AF/(-1*data.vss) + Idiff_Cl);
    
    dK_i = -(IKtot*data.AF/data.vmyo - CTKCl);
    
    %Derrivatives of the state variables:
    dyt(2) = dCa_i;
    dyt(3) = dCa_ss;
    dyt(4) = dCa_ssCaL;
    dyt(5) = Itr - Irel;
    dyt(6) = Iup - Itr*data.vjsr/data.vnsr-Ileak;
    betaCaMKII = 0.1 * data.betaCamK + (0.9 * data.betaCamK / 0.1371 * settings.Whole_cell_PP1);
    dyt(7) = data.alphaCamK*Cactive*(Cactive - Ctrap) - betaCaMKII*Ctrap;
    dyt(8) = dNa_i;
    dyt(9) = dNa_ss;
    dyt(10) = dK_i;
    dyt(11) = dCl_i;
    dyt(12) = dCl_ss;
    
    dyt(83) = (Cactive / (Cactive + data.KmCaMK) - y(83)) / data.tau_PLB_CaMKII;
    dyt(84) = (1.0 / (1.0 + (data.KmCaMK / Cactive)^2) - y(84)) / data.tau_RyR_CaMKII;
    dyt(85) = (Cactive / (Cactive + data.KmCaMK) - y(85)) / data.tau_PLB_CaMKII;
    dyt(86) = (Cactive / (Cactive + data.KmCaMK) - y(86)) / data.tau_PLB_CaMKII;
    dyt(87) = (Cactive / (Cactive + data.KmCaMK) - y(87)) / data.tau_PLB_CaMKII;
    dyt(88) = (Cactive / (Cactive + data.KmCaMK) - y(88)) / data.tau_RyR_CaMKII;
end

function [dyt, Ileak_JSR] = Ca_CICR(y, ICaL, settings, data)
    dyt = zeros(size(y,1), size(y,2));
        
	jsr = y(5);
    
    fLeakP = settings.ISO / (0.05 + settings.ISO);
	Ileak_JSR = (1 - settings.fRyRP) * data.k_JSR_Leak * exp(jsr / data.Km_JSR_Leak) * (jsr - y(3)) + settings.fRyRP * data.k_JSR_Leak_P * exp(jsr / data.Km_JSR_Leak_P) * (jsr - y(3));
    
	Irelcicr_NP = y(13);
    Irelcicr_P = y(80);
   
	Cactive = Ca_CaMKII(y, settings, data);
    
	deltaBetaCAMK = 2.000 * y(84);
	betatau = 0.6667 * 4.75 * (1 + deltaBetaCAMK);
	alpha_Rel = betatau * 0.1125;	
	Irelss_NP =  (ICaL) * (alpha_Rel) / (1 + ((1.0 / jsr)^8));
    tau_Irel_NP = betatau / (1 + (0.0123 / jsr));

    deltaBetaCAMK = 0.000 * y(84);
	betatau = 0.6667 * 4.75 * (1 + deltaBetaCAMK);
	alpha_Rel = betatau * 0.1125;
    Irelss_P = settings.RyRP_Amp * (ICaL) * (alpha_Rel) / (1 + ((1.0 / jsr)^8));
    tau_Irel_P = settings.RyRP_Tau * betatau / (1 + (0.0123 / jsr));
    
    dyt(13) = -(Irelss_NP + Irelcicr_NP) / tau_Irel_NP;
    dyt(80) = -(Irelss_P + Irelcicr_P) / tau_Irel_P;
end

function Cactive = Ca_CaMKII(y, settings, data)
    Ctrap = y(7);
    cass = y(3);
    Cactive = (1 - settings.CAMKIIB) * data.CaMK0*(1-Ctrap)/(1+data.Km/cass)+Ctrap;
end

function [cai] = Buff_Ca_i(ca_t, data, settings)
    kmtrpn = (1 - settings.fTnIP) * data.kmtrpn + settings.fTnIP * 1.5 * data.kmtrpn;
    
 	bmyo = data.cmdnbar+data.trpnbar-ca_t+kmtrpn+data.kmcmdn;
	cmyo = data.kmcmdn*kmtrpn -ca_t*(kmtrpn+data.kmcmdn)+data.trpnbar*data.kmcmdn+data.cmdnbar*kmtrpn;
	dmyo = -kmtrpn*data.kmcmdn*ca_t;
      
    cai =( 2*(bmyo.*bmyo-3*cmyo).^(1/2)/3).*cos(acos((9*bmyo.*cmyo-2*bmyo.*bmyo.*bmyo-27*dmyo)./(2*(bmyo.*bmyo-3*cmyo).^1.5))/3)-(bmyo/3);
end

function [cajsr] = Buff_Ca_jsr(ca_t, data)
    b=data.csqnbar+data.kmcsqn-ca_t;
    c=ca_t*data.kmcsqn;
    
    cajsr=-b/2+(b.^2+4*c).^(1/2)/2;
end

function [cass] = Buff_Ca_ss(ca_t,data)
	bmyo = data.BSLmax+data.BSRmax-ca_t+data.KmBSR+data.KmBSL;
	cmyo = data.KmBSL*data.KmBSR -ca_t*(data.KmBSR+data.KmBSL)+data.BSRmax*data.KmBSL+data.BSLmax*data.KmBSR;
	dmyo = -data.KmBSR*data.KmBSL*ca_t;      
    
    cass =( 2*(bmyo.*bmyo-3*cmyo).^(1/2)/3).*cos(acos((9*bmyo.*cmyo-2*bmyo.*bmyo.*bmyo-27*dmyo)./(2*(bmyo.*bmyo-3*cmyo).^1.5))/3)-(bmyo/3);
end