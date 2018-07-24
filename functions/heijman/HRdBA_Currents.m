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
% *   When called with the State matrix and settings structure, this      *
% *   function will calculate all currents, free calcium concentrations   *
% *   and phosphorylated fractions, for analysis purposes.                *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   constants_Signaling.m, constants_Electrophysiol.m                   *
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

function [currents] = HRdBA_Currents(State,t,Stm,settings)

    data = constants_Electrophysiol;
    
    V = State(:,1); Ca_i = State(:,2); Ca_ss = State(:,3); Ca_ssCaL = State(:,4);
    Ca_JSR = State(:,5); Ca_NSR = State(:,6); Ca_Trap = State(:,7); Na_i = State(:,8);
    Na_ss = State(:,9); K_i = State(:, 10); Cl_i = State(:,11); Cl_ss = State(:,12);
    Irel_NP = State(:,13); H = State(:, 14); m = State(:,15); J = State(:,16);
    mL = State(:,17); hL = State(:,18); to_a = State(:,19); to_if = State(:,20);
    to_is = State(:,21); xr = State(:,22);
    AA = State(:,23);

    H_PKAP = State(:,74); m_PKAP = State(:,75); J_PKAP = State(:,76); Irel_P = State(:,80);
    m_CaMKP = State(:,77); H_CaMKP = State(:,78); J_CaMKP = State(:,79);
    
    to_if_CaMKII = State(:,81); to_is_CaMKII = State(:,82);
    
    MM_ICaLNP_O1 = State(:,25); MM_ICaLNP_O2 = State(:,27);
    MM_ICaLP_O1 = State(:,50); MM_ICaLP_O2 = State(:,52);
    MM_IKsNP_O1 = State(:,47); MM_IKsNP_O2 = State(:,48);
    MM_IKsP_O1 = State(:,72); MM_IKsP_O2 = State(:,73);

    % Phosphorylation levels
    % ===================================================================
    Cactive=(1 - settings.CAMKIIB) .* data.CaMK0.*(1-Ca_Trap)./(1+data.Km./Ca_ss)+Ca_Trap;
    
    fPLBP_CaMKII = State(:,83);    fRyRP_CaMKII = State(:,84);
    fIToP_CaMKII = State(:,85);    fINaP_CaMKII = State(:,86);
    fIK1P_CaMKII = State(:,87);    fICaLP_CaMKII = State(:,88);

    if isfield(settings, 'fIToP_CaMKII') && settings.fIToP_CaMKII >= 0, fIToP_CaMKII = settings.fIToP_CaMKII; end
	if isfield(settings, 'fICaLP_CaMKII') && settings.fICaLP_CaMKII >= 0, fICaLP_CaMKII = settings.fICaLP_CaMKII; end
    
    if settings.runSignalingPathway == 1
        sigdata = settings.dataSignaling;
        fICaLP = min(max(((State(:,128) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1);
        fIKsP = min(max(((State(:,129) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);    
        fPLBP = min(max((State(:,130) - 6.591000e-001) / (9.945000e-001 - 6.591000e-001), 0), 1);    
        fTnIP = min(max((State(:,131) - 6.735188e-001) / (9.991797e-001 - 6.735188e-001), 0), 1);    
        fINaP = min(max((State(:,132) - 2.394795e-001) / (9.501431e-001 - 2.394795e-001), 0), 1);    
        fINaKP = min(max((State(:,133) - 1.263453e-001) / (9.980137e-001 - 1.263453e-001), 0), 1);    
        fRyRP = min(max(((State(:,134) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);            
        fIKurP = min(max((State(:,135) - 5.893798e-002) / (3.937470e-001 - 5.893798e-002), 0), 1);        
    else
        fICaLP = settings.fICaLP .* ones(length(t),1);
        fIKsP = settings.fIKsP .* ones(length(t),1);
        fPLBP = settings.fPLBP .* ones(length(t),1);
        fTnIP = settings.fTnIP .* ones(length(t),1);
        fINaP = settings.fINaP .* ones(length(t),1);
        fINaKP = settings.fINaKP .* ones(length(t),1);
        fRyRP = settings.fRyRP .* ones(length(t),1);
        fIKurP = settings.fIKurP .* ones(length(t),1);
    end
    
    % Buffers
    % ===================================================================
    Ca_i = Buff_Ca_i(Ca_i, data, settings, fTnIP);
	if isfield(settings, 'Cai_clamp')
        if settings.Cai_clamp >= 0, Ca_i = settings.Cai_clamp * ones(size(Ca_i,1),size(Ca_i,2)); end
	end
    Ca_ss = Buff_Ca_ss(Ca_ss, data);
    Ca_ssCaL = Buff_Ca_ss(Ca_ssCaL, data);
    Ca_JSR = Buff_Ca_jsr(Ca_JSR, data);
    
    % Reversal potential
    % ===================================================================
    ENa=log(settings.Na_o./Na_i)/data.frt;
    EK=log(settings.K_o./K_i)/data.frt;
    ECl=-log(settings.Cl_o./Cl_i)/data.frt;
    EKs = log((settings.K_o+data.prnak*settings.Na_o)./(K_i+data.prnak*Na_i))/data.frt;
    
    % Sodium currents
    % ===================================================================
    GNaNP = settings.INaNPParams(1).*m.*m.*m.*H.*J;
    INaNP = GNaNP.*(V-ENa);
    INaPKAP = 1.25 * settings.INaNPParams(1).*m_PKAP.*m_PKAP.*m_PKAP.*H_PKAP.*J_PKAP .* (V - ENa);
    INaCaMKP = settings.INaNPParams(1).*m_CaMKP.*m_CaMKP.*m_CaMKP.*H_CaMKP.*J_CaMKP .* (V - ENa);    
    INaBothP = settings.INaNPParams(1).*m_PKAP.*m_PKAP.*m_PKAP.*H_CaMKP.*J_CaMKP .* (V - ENa);    
    
	f_INa_BothP = fINaP_CaMKII .* fINaP;
	f_INa_CaMK_Only = fINaP_CaMKII - f_INa_BothP;
	f_INa_PKA_Only = fINaP - f_INa_BothP;
	f_INa_NP = -f_INa_CaMK_Only - f_INa_PKA_Only - f_INa_BothP + 1;

	INa = (settings.INaB) * (f_INa_NP .* INaNP + f_INa_PKA_Only .* INaPKAP + f_INa_CaMK_Only .* INaCaMKP + f_INa_BothP .* INaBothP);
    
    dGNaLCaMK = 0.0095 * fINaP_CaMKII;
    INaL = (settings.INaLB) * (data.GNaL + dGNaLCaMK).*mL.^3.*hL.*(V-ENa);
    [INaCa, INaCass] = Ion_HH_NaCa(V, Na_i, Ca_i, Na_ss, Ca_ss, settings, data);       
    INaK_NP = data.ibarnak .* (1 ./ (1 + exp(-((V + 92) .* data.frt)))) .* ((Na_i ./ (Na_i + data.kmnai_NP)).^3) .* (settings.K_o / (settings.K_o + data.kmko));
    INaK_P = data.ibarnak .* (1 ./ (1 + exp(-((V + 92) .* data.frt)))) .* ((Na_i ./ (Na_i + data.kmnai_P)).^3) .* (settings.K_o / (settings.K_o + data.kmko));        
    INaK = (settings.INaKB) .* (1 - fINaKP) .* INaK_NP + fINaKP .* INaK_P;
    INab = (settings.INabB) * data.F .* 0.32E-8 .* V .* data.frt .* ( (Na_i .* exp(V .* data.frt) - settings.Na_o) ./ (exp(V .* data.frt) - 1) );

    % Potassium currents
    % ===================================================================
    IKp_NP = data.GKpmax.*(V-EK)./(1+exp((15.0-V)./17));
    IKp_P = 3.62 * data.GKpmax .* (V - EK) ./ (1 + exp((36.0 - V)./17));
    IKp = (settings.IKpB) * ((1 - fIKurP) .* IKp_NP + fIKurP .* IKp_P); 
    ak1 = 1.02./(1+exp(0.2385*(V-EK-59.215)));
    bk1 = (0.49124*exp(0.08032*(V-EK+5.476))+exp(0.06175*(V-EK-594.31)))./(1+exp(-0.5143*(V-EK+4.753)));
    IK1 = (settings.IK1B) .* ((1 - fIK1P_CaMKII) .* data.GK1max + fIK1P_CaMKII .* 1.2 .* data.GK1max) .* sqrt(settings.K_o/5.4) .* (ak1./(ak1+bk1)) .* (V-EK);
    IKr = (settings.IKrB) * (data.GKrmax*sqrt(settings.K_o/5.4)) .* xr .* (1./(1+exp((V+10)/15.4))) .* (V-EK);
    ITo1_NP = data.gitodv .* to_a.^3 .* (0.7356 .* to_if + 0.2644 .* to_is) .* exp(V ./ 550) .* (V - EK);
    ITo1_CaMKII = data.gitodv .* to_a.^3 .* (0.7356 .* to_if_CaMKII + 0.2644 .* to_is_CaMKII) .* exp(V ./ 550) .* (V - EK);
    ITo1 = (settings.ITo1B) .* ((1 - fIToP_CaMKII) .* ITo1_NP + fIToP_CaMKII .* ITo1_CaMKII);    
    IKsNP = settings.IKsNPParams(end) .* (1+0.6./(1+(((3.8E-05)./Ca_i).^1.4))) .* (MM_IKsNP_O1 + MM_IKsNP_O2) .* (V-EKs);
    IKsP = settings.IKsPParams(end) .* (1+0.6./(1+(((3.8E-05)./Ca_i).^1.4))) .* (MM_IKsP_O1 + MM_IKsP_O2) .* (V-EKs);
    IKs = (settings.IKsB) .* (fIKsP .* IKsP + (1 - fIKsP) .* IKsNP);
    
    % Chloride currents
    % ===================================================================
    CTKCl = (1 - settings.CTKClB) * data.CTKClmax .* (EK-ECl)./((EK-ECl)+87.8251);
    CTNaCl = (1 - settings.CTNaClB) * data.CTNaClmax .* (ENa-ECl).^4.0 ./ ((ENa-ECl).^4.0+87.8251^4.0);
    Ito2_max = data.PCl .* V .* data.F .* data.frt .* (Cl_i-settings.Cl_o.*exp(V*data.frt))./(1-exp(V*data.frt));
    Irel_adjust = (1 - fRyRP) .* Irel_NP + fRyRP .* Irel_P;
    Ito2 = (settings.ITo2B) * Ito2_max .* AA .* (1 - (1 ./ (1 + (Irel_adjust ./ 0.4).^2)));
    IClb = (settings.IClB) * data.GClb.*(V-ECl);
    
    % Calcium currents
    % ===================================================================    
    ICaLBarNP = settings.LCCNPParams(22) .* (1 + 0.4 .* fICaLP_CaMKII) .* 4 .* V .* data.F .* data.frt .* ((Ca_ssCaL .* exp(2*V*data.frt) - 0.341 .* settings.Ca_o) ./ (exp(2*V*data.frt) - 1));    
    ICaLBarP = settings.LCCPParams(22) .* (1 + 0.1 .* fICaLP_CaMKII) .* 4 .* V .* data.F .* data.frt .* ((Ca_ssCaL .* exp(2*V*data.frt) - 0.341 .* settings.Ca_o) ./ (exp(2*V*data.frt) - 1));
    ICaL = (settings.ICaLB) * ((1 - fICaLP) .* ICaLBarNP .* (MM_ICaLNP_O1 + MM_ICaLNP_O2) + fICaLP .* ICaLBarP .* (MM_ICaLP_O1 + MM_ICaLP_O2));
    ICab = (settings.ICabB) * (1.995E-7)*4*V*data.F*data.frt.*(Ca_i.*exp(2*V*data.frt)-0.341*settings.Ca_o)./(exp(2*V*data.frt)-1);
    IpCa = (settings.IpCaB) * data.ibarpca .* Ca_i ./ (data.kmpca+Ca_i);
    
    Itr = (1 - settings.ItrB) * (Ca_NSR-Ca_JSR)/25;
    Idiff_ssCaL = (Ca_ss - Ca_ssCaL) / 2.0;
    Idiff_Ca = (1 - settings.IdiffB) * (Ca_ss - Ca_i) / 0.2;
    Idiff_Na = (1 - settings.IdiffB) * (Na_ss - Na_i) / 0.2;
    Idiff_Cl = (1 - settings.IdiffB) * (Cl_ss - Cl_i) / 0.2;    

	fPLBNP = 1 - fPLBP - fPLBP_CaMKII + fPLBP .* fPLBP_CaMKII;
    fPLBCAMK = fPLBP_CaMKII - fPLBP .* fPLBP_CaMKII;
    fPLBPKA = fPLBP - fPLBP .* fPLBP_CaMKII;
    fPLBBoth = fPLBP .* fPLBP_CaMKII;        
    
    iupmaxCAMK = 3.25 * data.iupmax;
    frac_SERCA2a_P = 1 ./ (1 + (0.03 ./ Cactive) .^ 2);    
    iup = (settings.IupB) .* ((1 - frac_SERCA2a_P) .* data.iupmax + frac_SERCA2a_P .* iupmaxCAMK) .* (Ca_i ./ (Ca_i + (fPLBNP .* data.Kmup + fPLBCAMK .* (data.Kmup - data.dKmPLBCAMK) + fPLBPKA .* (data.Kmup - data.dKmPLBPKA) + fPLBBoth .* (data.Kmup - data.dKmPLBBoth))));
    Ileak = (settings.IleakB) .* (((1 - frac_SERCA2a_P) * data.iupmax + frac_SERCA2a_P .* iupmaxCAMK) ./ data.nsrmax) .* Ca_NSR;
        
    fLeakP = settings.ISO / (0.05 + settings.ISO);
	Ileak_JSR = (1 - fRyRP) .* data.k_JSR_Leak .* exp(Ca_JSR ./ data.Km_JSR_Leak) .* (Ca_JSR - Ca_ss) + fRyRP .* data.k_JSR_Leak_P .* exp(Ca_JSR ./ data.Km_JSR_Leak_P) .* (Ca_JSR - Ca_ss);
    Irel_adjust = (1 - fRyRP) .* Irel_NP + fRyRP .* Irel_P + Ileak_JSR;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    currents.Cafree=Ca_i;
    currents.JSRfree=Ca_JSR;
    currents.Cassfree=Ca_ss;
    currents.CassCaLfree=Ca_ssCaL;
    currents.inab=INab;
    currents.ical=ICaL;
    currents.icalbarNP = ICaLBarNP;
    currents.icalbarP = ICaLBarP;
    currents.icab=ICab;
    currents.inaca=INaCa;
    currents.inacass=INaCass;
    currents.itr=Itr;
    currents.iup=iup;
    currents.ipca=IpCa;
    currents.camka=Cactive;
    currents.iks=IKs;
    currents.ikr=IKr;
    currents.ik1=IK1;
    currents.ikp=IKp;
    currents.ina=INa;
%     currents.gna=GNa;
    currents.inal=INaL;
    currents.irel=Irel_adjust;
    currents.ito2=Ito2;
    currents.CTNaCl=CTNaCl;
    currents.CTKCl=CTKCl;
    currents.iclb=IClb;
    currents.ito=ITo1;
    currents.inak=INaK;
    currents.ileak=Ileak;
    currents.ileak_JSR = Ileak_JSR;    
    currents.stm=Stm;
    currents.caiont =ICaL+ICab+IpCa-2*INaCa;%
    currents.naiont = INa+3*INaCa+3*INaK+INaL+INab;
    currents.clont=IClb+Ito2;
    currents.fICaLP = fICaLP;
    currents.fIKsP = fIKsP;    
    currents.fPLBP = fPLBP;
    currents.fTnIP = fTnIP;
    currents.fINaP = fINaP;
    currents.fINaKP = fINaKP;
    currents.fRyRP = fRyRP;
    currents.fIKurP = fIKurP;
    currents.kiont =IKr+IKs+IK1+IKp-2*INaK+ITo1+Stm';
    currents.DVDT=currents.caiont+currents.naiont+currents.kiont+currents.clont;
end

function [INaCa, INaCass] = Ion_HH_NaCa(V, Na_i, Ca_i, Na_ss, Ca_ss, settings, data)   
    num = 0.8 .* data.NCXmax .* (Na_i.^3 .* settings.Ca_o .* exp(data.eta .* V .* data.frt) - settings.Na_o^3 .* Ca_i .* exp((data.eta-1) .* V .* data.frt));
    denom1 = 1 + (data.KmCa ./ Ca_i).^2;
    denom2 = 1 + data.ksat .* exp((data.eta - 1) .* V .* data.frt);
    denom3 = data.KmCao .* Na_i.^3 + data.KmNao^3 .* Ca_i + data.KmNai^3 .* settings.Ca_o .* (1 + Ca_i./data.KmCai);
    denom4 = data.KmCai .* settings.Na_o^3 * (1 + (Na_i./data.KmNai).^3) + Na_i.^3 .* settings.Ca_o + settings.Na_o^3 .* Ca_i;
    
    INaCa = (settings.INaCaB) * (num) ./ (denom1 .* denom2 .* (denom3 + denom4));
    
    num = 0.2 .* data.NCXmax .* (Na_ss.^3 .* settings.Ca_o .* exp(data.eta .* V .* data.frt) - settings.Na_o.^3 .* Ca_ss .* exp((data.eta-1) .* V .* data.frt));
    denom1 = 1 + (data.KmCa ./ Ca_ss).^2;
    denom2 = 1 + data.ksat .* exp((data.eta - 1) .* V .* data.frt);
    denom3 = data.KmCao .* Na_ss.^3 + data.KmNao^3 .* Ca_ss + data.KmNai^3 .* settings.Ca_o .* (1 + Ca_ss./data.KmCai);
    denom4 = data.KmCai .* settings.Na_o^3 .* (1 + (Na_ss./data.KmNai).^3) + Na_ss.^3 .* settings.Ca_o + settings.Na_o^3 .* Ca_ss;
    
    INaCass = (settings.INaCaB) * (num) ./ (denom1 .* denom2 .* (denom3 + denom4));    
end


function [cai] = Buff_Ca_i(ca_t, data, settings, fTnIP)
    kmtrpn = (1 - fTnIP) .* data.kmtrpn + fTnIP .* 1.5 .* data.kmtrpn;
    
 	bmyo = data.cmdnbar+data.trpnbar-ca_t+kmtrpn+data.kmcmdn;
	cmyo = data.kmcmdn.*kmtrpn -ca_t.*(kmtrpn+data.kmcmdn)+data.trpnbar.*data.kmcmdn+data.cmdnbar.*kmtrpn;
	dmyo = -kmtrpn*data.kmcmdn.*ca_t;
      
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