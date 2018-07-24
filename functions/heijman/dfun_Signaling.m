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
% *   function will calculate dy/dt for the adrenergic signaling pathway  *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   constants_Signaling.m                                               *
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

function dy = dfun_Signaling(t, y, flags, settings)
% State vector contains:
%  1: Gs_aGTP_CAV
%  2: Gs_aGTP_ECAV
%  3: Gs_a_GTP_CYT
%  4: Gs_bg_CAV
%  5: Gs_bg_ECAV
%  6: Gs_bg_CYT
%  7: Gs_aGDP_CAV
%  8: Gs_aGDP_ECAV
%  9: Gs_aGDP_CYT
% 10: cAMP_CAV
% 11: cAMP_ECAV
% 12: cAMP_CYT
% 13: R_pkap_tot_CAV
% 14: R_pkap_tot_ECAV
% 15: R_pkap_tot_CYT
% 16: R_grkp_tot_CAV
% 17: R_grkp_tot_ECAV
% 18: R_grkp_tot_CYT

% 19: RLC_CAV
% 20: L2RC_CAV
% 21: L2R_CAV
% 22: C_CAV
% 23: PKI_CAV
% 24: RLC_ECAV
% 25: L2RC_ECAV
% 26: L2R_ECAV
% 27: C_ECAV
% 28: PKI_ECAV
% 29: RLC_CYT
% 30: L2RC_CYT
% 31: L2R_CYT
% 32: C_CYT
% 33: PKI_CYT

% 34: PDE3_P_CAV
% 35: PDE3_P_CYT
% 36: PDE4_P_CAV
% 37: PDE4_P_ECAV
% 38: PDE4_P_CYT

% 39: Inhib1_P_CYT

% 40: fLCC_P
% 41: fIKS_P
% 42: fPLB_P
% 43: fTnI_P
% 44: fINa_P
% 45: fINaK_P
% 46: fRyR_P
% 47: fIKur_P

% 48: Rb2_pkap_tot_CAV
% 49: Rb2_grkp_tot_CAV
% 50: Gi_aGTP_CAV
% 51: Gi_bg_CAV
% 52: Gi_aGDP_CAV
% 53: Rb2_pkap_tot_ECAV
% 54: Rb2_grkp_tot_ECAV
% 55: Gi_aGTP_ECAV
% 56: Gi_bg_ECAV
% 57: Gi_aGDP_ECAV

L = settings.ISO;
data = settings.dataSignaling;

dy = zeros(size(y,1), size(y,2));

Gs_abg_CAV = (data.f_Gs_CAV * data.Gs_tot * (data.V_tot / data.V_CAV)) - y(1) - y(7);
Gi_abg_CAV = (data.f_Gi_CAV * data.Gi_tot * (data.V_tot / data.V_CAV)) - y(50) - y(52);
R_b1_comp = data.f_Rb1_CAV * data.R_b1_tot * (data.V_tot / data.V_CAV);
R_b2_comp = data.f_Rb2_CAV * data.R_b2_tot * (data.V_tot / data.V_CAV);
[Rb1Gs_CAV, LRb1Gs_CAV, Rb2Gs_CAV, LRb2Gs_CAV, Rb2Gi_CAV, LRb2Gi_CAV, dRb1_pkap_tot, dRb1_grkp_tot, dRb2_pkap_tot, dRb2_grkp_tot] = Mod_LigRecepGprotWithBeta2(R_b1_comp, y(13), y(16), R_b2_comp, y(48), y(49), L, Gs_abg_CAV, Gi_abg_CAV, y(22), data.GRK_CAV, data, settings);

dy(13) = dRb1_pkap_tot;
dy(16) = dRb1_grkp_tot;
dy(48) = dRb2_pkap_tot;
dy(49) = dRb2_grkp_tot;

Gs_abg_ECAV = (data.f_Gs_ECAV * data.Gs_tot * (data.V_tot / data.V_ECAV)) - y(2) - y(8);
Gi_abg_ECAV = ((1 - data.f_Gi_CAV) * data.Gi_tot * (data.V_tot / data.V_ECAV)) - y(55) - y(57);
R_b1_comp = data.f_Rb1_ECAV * data.R_b1_tot * (data.V_tot / data.V_ECAV);
R_b2_comp = (1 - data.f_Rb2_CAV) * data.R_b2_tot * (data.V_tot / data.V_ECAV);
[Rb1Gs_ECAV, LRb1Gs_ECAV, Rb2Gs_ECAV, LRb2Gs_ECAV, Rb2Gi_ECAV, LRb2Gi_ECAV, dRb1_pkap_tot, dRb1_grkp_tot, dRb2_pkap_tot, dRb2_grkp_tot] = Mod_LigRecepGprotWithBeta2(R_b1_comp, y(14), y(17), R_b2_comp, y(53), y(54), L, Gs_abg_ECAV, Gi_abg_ECAV, y(27), data.GRK_ECAV, data, settings);
dy(14) = dRb1_pkap_tot;
dy(17) = dRb1_grkp_tot;
dy(53) = dRb2_pkap_tot;
dy(54) = dRb2_grkp_tot;

Gs_f_CYT = ((1 - data.f_Gs_CAV - data.f_Gs_ECAV) * data.Gs_tot * (data.V_tot / data.V_CYT)) - y(3) - y(9);
R_comp = (1 - data.f_Rb1_CAV - data.f_Rb1_ECAV) * data.R_b1_tot * (data.V_tot / data.V_CYT);
[RG_CYT, LRG_CYT, dR_pkap_tot, dR_grk_tot] = Mod_LigRecepGprot(R_comp, y(15), y(18), L, Gs_f_CYT, y(32), data.GRK_CYT, data, settings);
dy(15) = dR_pkap_tot;
dy(18) = dR_grk_tot;

%% ========================================================================

dy = Mod_GprotAct(y, dy, Rb1Gs_CAV + data.k_GsAct_b2 * Rb2Gs_CAV, LRb1Gs_CAV + data.k_GsAct_b2 * LRb2Gs_CAV, Rb2Gi_CAV, LRb2Gi_CAV, Rb1Gs_ECAV + data.k_GsAct_b2 * Rb2Gs_ECAV, LRb1Gs_ECAV + data.k_GsAct_b2 * LRb2Gs_ECAV, Rb2Gi_ECAV, LRb2Gi_ECAV, RG_CYT, LRG_CYT, data);
[dcAMPAC56_CAV, dcAMPAC56_CYT, dcAMPAC47_ECAV, dcAMPAC47_CYT] = Mod_AC(y, data);

%% ========================================================================

[dy, dcAMP_PKA_CAV, dcAMP_PKA_ECAV, dcAMP_PKA_CYT] = Mod_PKA(y, dy, data);

dy = Mod_cAMP(y, dy, dcAMPAC56_CAV, dcAMPAC56_CYT, dcAMPAC47_ECAV, dcAMPAC47_CYT, dcAMP_PKA_CAV, dcAMP_PKA_ECAV, dcAMP_PKA_CYT, data);
dy = Mod_PDE_Phosphorylation(y, dy, data);
[dy, PP1_CYT] = Mod_PP1_Inhibition(y, dy, data);
dy = Mod_Channel_Phosphorylation(y, dy, PP1_CYT, data);

end

function [RG, LRG, dR_pkap_tot, dR_grkp_tot] = Mod_LigRecepGprot(R_tot, R_pkap_tot, R_grkp_tot, L, G_abg, PKAC, GRK, data, settings)        
    R_np_tot = R_tot - R_pkap_tot - R_grkp_tot;

    % NP distribution
    a = ((data.K_b1_h + L)*(data.K_b1_l + L))/data.K_b1_l;
    b = G_abg*data.K_b1_h - R_np_tot*(data.K_b1_h + L) + G_abg*L + data.K_b1_c*data.K_b1_h + (data.K_b1_c*data.K_b1_h*L)/data.K_b1_l;
    c = -R_np_tot * data.K_b1_c * data.K_b1_h;

%    a = data.K_b1_h * data.K_b1_l + data.K_b1_l * L + data.K_b1_h * L + L*L;
%    b = data.K_b1_l * data.K_b1_c * L + data.K_b1_l * data.K_b1_c * data.K_b1_h + G_abg * (L + data.K_b1_h) - data.K_b1_h * data.K_b1_l * L * R_np_tot;
%    c = -R_np_tot * data.K_b1_l * data.K_b1_c * data.K_b1_h;        
    
    Rf = (-b + sqrt(b*b - 4 * a * c)) / (2 * a);
    Gf = G_abg / (1 + Rf / data.K_b1_c + L * Rf / (data.K_b1_c * data.K_b1_h));    
    LR = (L * Rf) / data.K_b1_l;
    RG = (Rf * Gf) / data.K_b1_c;
    LRG = (L * Rf * Gf) / (data.K_b1_h * data.K_b1_c);
    
    % Group state derivatives
    dR_pkap_tot = data.rate_bds * data.k_b1_pkap * PKAC * R_np_tot - data.rate_bds * data.k_b1_dp * R_pkap_tot;
    dR_grkp_tot = data.rate_bds * data.k_b1_grkp * GRK * (LR + LRG) - data.rate_bds * data.k_b1_grkdp * R_grkp_tot;
    
%     if isfield(settings, 'f_Beta1_Block'), RG = RG * (1 - settings.f_Beta1_Block); end
    if isfield(settings, 'f_Beta1_Block')
        RG = RG + settings.f_Beta1_Block * LRG; 
        LRG = LRG * (1 - settings.f_Beta1_Block); 
    end
end

function [Rb1Gs, LRb1Gs, Rb2Gs, LRb2Gs, Rb2Gi, LRb2Gi, dRb1_pkap_tot, dRb1_grkp_tot, dRb2_pkap_tot, dRb2_grkp_tot] = Mod_LigRecepGprotWithBeta2(Rb1_tot, Rb1_pkap_tot, Rb1_grkp_tot, Rb2_tot, Rb2_pkap_tot, Rb2_grkp_tot, L, Gs_abg, Gi_abg, PKAC, GRK, data, settings)          
    Rb1_np_tot = Rb1_tot - Rb1_pkap_tot - Rb1_grkp_tot;
    Rb2_np_tot = Rb2_tot - Rb2_pkap_tot - Rb2_grkp_tot;
    
    % Beta 2 PKA phosphorylated distribution, coupling to Gi
    KF2 = data.K_b2_f; KA2 = data.K_b2_a; KN2 = data.K_b2_n;
    a = ((KF2 + L)*(KN2 + L))/KN2;
    b = Gi_abg*KF2 - Rb2_pkap_tot*(KF2 + L) + Gi_abg*L + KA2*KF2 + (KA2*KF2*L)/KN2;
    c = -Rb2_pkap_tot * KA2 * KF2;
    
    Rb2_pkap_f = (-b + sqrt(b*b - 4*a*c)) / (2*a);
    Gi_f = Gi_abg / (1 + Rb2_pkap_f / KA2 + L * Rb2_pkap_f / (KA2 * KF2));
    Rb2Gi = Rb2_pkap_f * Gi_f / KA2;
    LRb2Gi = L * Rb2Gi / KF2;
        
    % Beta 1 and Beta 2 Nonphosphorylated receptors distribution
    % a = KL1*KL2*(KH1 + L)*(KH2 + L), so b, c and d are divided by this
    KL1 = data.K_b1_l; KC1 = data.K_b1_c; KH1 = data.K_b1_h;
    KL2 = data.K_b2_l; KC2 = data.K_b2_c; KH2 = data.K_b2_h;
    
    b = Rb1_np_tot + (KC1*KH1*KL2*L^2 - Gs_abg*(KL1*KL2*L^2 + KH1*KH2*KL1*KL2 + KH1*KL1*KL2*L + KH2*KL1*KL2*L) + KC2*KH2*KL1*L^2 + KL1*KL2*L^2*Rb2_np_tot + KC1*KH1*KH2*KL1*KL2 + KC2*KH1*KH2*KL1*KL2 + KC1*KH1*KH2*KL2*L + KC2*KH1*KH2*KL1*L + KC1*KH1*KL1*KL2*L + KC2*KH2*KL1*KL2*L + KH1*KH2*KL1*KL2*Rb2_np_tot + KH1*KL1*KL2*L*Rb2_np_tot + KH2*KL1*KL2*L*Rb2_np_tot)/(KL1*KL2*L^2 + KH1*KH2*KL1*KL2 + KH1*KL1*KL2*L + KH2*KL1*KL2*L);
    c = (KC1*KH1*L^2*Rb2_np_tot - Gs_abg*KC1*KH1*L^2 - Gs_abg*KC1*KH1*KH2*KL1 - Gs_abg*KC2*KH1*KH2*KL1 - Gs_abg*KC1*KH1*KH2*L - Gs_abg*KC1*KH1*KL1*L - Gs_abg*KC2*KH2*KL1*L + KC1*KC2*KH1*KH2*KL1 + KC1*KC2*KH1*KH2*L + KC1*KH1*KH2*KL1*Rb2_np_tot + KC2*KH1*KH2*KL1*Rb1_np_tot + KC1*KH1*KH2*L*Rb2_np_tot + KC1*KH1*KL1*L*Rb2_np_tot + KC2*KH2*KL1*L*Rb1_np_tot)/(KL1*(KH1 + L)*(KH2 + L)) + (KC1*KC2*KH1*KH2*L^2 - Gs_abg*KC2*KH2*KL1*L^2 + KC2*KH2*KL1*L^2*Rb1_np_tot - Gs_abg*KC2*KH1*KH2*KL1*L + KC1*KC2*KH1*KH2*KL1*L + KC2*KH1*KH2*KL1*L*Rb1_np_tot)/(KL1*KL2*(KH1 + L)*(KH2 + L));
    d = (Gs_abg*KC1*KC2*KH1*KH2*(KL1 + L)*(KL2 + L))/(KL1*KL2*(KH1 + L)*(KH2 + L));

    V = 1/3 * (3 * c - b*b);
    W = 1/27 * (2 * b^3 - 9*b*c -27*d);
    phi= acos(sqrt((W^2/4)/(-V^3/27)));
    sol = 2*sqrt(-V/3)*cos(phi/3)-b/3;
    
    phi = (0.5 * d + (b * c / 6.0) - (b * b * b / 27) + (-(b ^ 3 * d / 27) - (b ^2 * c^2 / 108) + (b * c * d / 6) + (c ^ 3 / 27) + (d ^ 2 / 4)) ^ 0.5) ^ (1 / 3.0);
    Gs_f = abs((phi - (1 / phi) * ((c / 3) - (b * b / 9)) - (b / 3)));
    Rb1_f = Rb1_np_tot / (1 + L/KL1 + (1 / KC1 + L / (KC1 * KH1)) * Gs_f);
    Rb2_f = Rb2_np_tot / (1 + L/KL2 + (1 / KC2 + L / (KC2 * KH2)) * Gs_f);
    
    LRb1 = (L * Rb1_f) / KL1;
    Rb1Gs = (Rb1_f * Gs_f) / KC1;
    LRb1Gs = (L * Rb1_f * Gs_f) / (KH1 * KC1);    
    
    LRb2 = (L * Rb2_f) / KL2;
    Rb2Gs = (Rb2_f * Gs_f) / KC2;
    LRb2Gs = (L * Rb2_f * Gs_f) / (KH2 * KC2);
        
    % Group state derivatives
    dRb1_pkap_tot = data.rate_bds * data.k_b1_pkap * PKAC * Rb1_np_tot - data.rate_bds * data.k_b1_dp * Rb1_pkap_tot;
    dRb1_grkp_tot = data.rate_bds * data.k_b1_grkp * GRK * (LRb1 + LRb1Gs) - data.rate_bds * data.k_b1_grkdp * Rb1_grkp_tot;
    dRb2_pkap_tot = data.rate_bds * data.k_b1_pkap * PKAC * Rb2_np_tot - data.rate_bds * data.k_b1_dp * Rb2_pkap_tot;
    dRb2_grkp_tot = data.rate_bds * data.k_b1_grkp * GRK * (LRb2 + LRb2Gs) - data.rate_bds * data.k_b1_grkdp * Rb2_grkp_tot;
    
%     if isfield(settings, 'f_Beta1_Block'), Rb1Gs = Rb1Gs * (1 - settings.f_Beta1_Block); end
    if isfield(settings, 'f_Beta1_Block')
        Rb1Gs = Rb1Gs + settings.f_Beta1_Block * LRb1Gs;
        LRb1Gs = LRb1Gs * (1 - settings.f_Beta1_Block); 
    end
%     if isfield(settings, 'f_Beta2_Block'), Rb2Gs = Rb2Gs * (1 - settings.f_Beta2_Block); end
    if isfield(settings, 'f_Beta2_Block')
        Rb2Gs = Rb2Gs + settings.f_Beta2_Block * LRb2Gs;
        LRb2Gs = LRb2Gs * (1 - settings.f_Beta2_Block);
        
        Rb2Gi = Rb2Gi + settings.f_Beta2_Block * LRb2Gi;
        LRb2Gi = LRb2Gi * (1 - settings.f_Beta2_Block);
    end
%     if isfield(settings, 'f_Beta2_Block'), Rb2Gi = Rb2Gi * (1 - settings.f_Beta2_Block); end
    if isfield(settings, 'f_Beta2_Block'),  end
end

function dy = Mod_GprotAct(y, dy, RGsTot_CAV, LRGsTot_CAV, Rb2Gi_CAV, LRb2Gi_CAV, RG_ECAV, LRG_ECAV, Rb2Gi_ECAV, LRb2Gi_ECAV, RG_CYT, LRG_CYT, data)
    dy(1) = RGsTot_CAV * data.kact2_Gs + LRGsTot_CAV * data.kact1_Gs - y(1) * data.khydr_Gs;
    dy(50) = Rb2Gi_CAV * data.kact2_Gi + LRb2Gi_CAV * data.kact1_Gi - y(50) * data.khydr_Gi;
    dy(2) = RG_ECAV * data.kact2_Gs + LRG_ECAV * data.kact1_Gs - y(2) * data.khydr_Gs;
    dy(55) = Rb2Gi_ECAV * data.kact2_Gi + LRb2Gi_ECAV * data.kact1_Gi - y(55) * data.khydr_Gi;
    dy(3) = RG_CYT * data.kact2_Gs + LRG_CYT * data.kact1_Gs - y(3) * data.khydr_Gs;
    
    dy(4) = RGsTot_CAV * data.kact2_Gs + LRGsTot_CAV * data.kact1_Gs - y(4) * y(7) * data.kreas_Gs;
    dy(51) = Rb2Gi_CAV * data.kact2_Gi + LRb2Gi_CAV * data.kact1_Gi - y(51) * y(52) * data.kreas_Gi;    
    dy(5) = RG_ECAV * data.kact2_Gs + LRG_ECAV * data.kact1_Gs - y(5) * y(8) * data.kreas_Gs;
    dy(56) = Rb2Gi_ECAV * data.kact2_Gi + LRb2Gi_ECAV * data.kact1_Gi - y(56) * y(57) * data.kreas_Gi;
    dy(6) = RG_CYT * data.kact2_Gs + LRG_CYT * data.kact1_Gs - y(6) * y(9) * data.kreas_Gs;
    
    dy(7) = y(1) * data.khydr_Gs - y(4) * y(7) * data.kreas_Gs;
    dy(52) = y(50) * data.khydr_Gi - y(51) * y(52) * data.kreas_Gi;
    dy(8) = y(2) * data.khydr_Gs - y(5) * y(8) * data.kreas_Gs;
    dy(57) = y(55) * data.khydr_Gi - y(56) * y(57) * data.kreas_Gi;
    dy(9) = y(3) * data.khydr_Gs - y(6) * y(9) * data.kreas_Gs;    
end

function [dcAMPAC56_CAV, dcAMPAC56_CYT, dcAMPAC47_ECAV, dcAMPAC47_CYT] = Mod_AC(y, data)
    kAC56_CAV = (data.AC56_basal + (y(1) ^ data.AC56_hill_Gs) / (data.AC56_Km_Gs + y(1) ^ data.AC56_hill_Gs)) * (1 - (1 - data.AC56_V_GsGi * y(1) ^ data.AC56_hill_GsGi / (data.AC56_Km_GsGi + y(1) ^ data.AC56_hill_GsGi)) * y(51) / (data.AC56_Km_Gi + y(51)));
    dcAMPAC56_CAV = (kAC56_CAV * data.AC56_CAV * data.AF56 * data.ATP) / (data.KmATP + data.ATP);

    kAC56_CYT = (data.AC56_basal + (y(3) ^ data.AC56_hill_Gs) / (data.AC56_Km_Gs + y(3) ^ data.AC56_hill_Gs)) * (1 - (1 - data.AC56_V_GsGi * y(3) ^ data.AC56_hill_GsGi / (data.AC56_Km_GsGi + y(3) ^ data.AC56_hill_GsGi)) * 0.0 / (data.AC56_Km_Gi + 0.0));
    dcAMPAC56_CYT = (kAC56_CYT * data.AC56_CYT * data.AF56 * data.ATP) / (data.KmATP + data.ATP);
    
    kAC47_ECAV = (data.AC47_basal + (y(2) ^ data.AC47_hill_Gs) / (data.AC47_Km_Gs + y(2) ^ data.AC47_hill_Gs)) .* (1 + (data.AC47_Vmax_Gi .* y(56) .^ data.AC47_hill_Gi) ./ (data.AC47_Km_Gi + y(56) .^ data.AC47_hill_Gi));    
    dcAMPAC47_ECAV = (kAC47_ECAV * data.AC47_ECAV * data.AF47 * data.ATP) / (data.KmATP + data.ATP);
    
    kAC47_CYT = data.AF47*(data.AC47_basal + (y(3) ^ data.AC47_hill_Gs) / (data.AC47_Km_Gs + y(3) ^ data.AC47_hill_Gs));    
    dcAMPAC47_CYT = (kAC47_CYT * data.AC47_CYT * data.ATP) / (data.KmATP + data.ATP);
end

function dy = Mod_cAMP(y, dy, dcAMPAC56_CAV, dcAMPAC56_CYT, dcAMPAC47_ECAV, dcAMPAC47_CYT, dcAMP_PKA_CAV, dcAMP_PKA_ECAV, dcAMP_PKA_CYT, data)
    PDE2_CAV = data.f_PDE2_CAV * data.PDE2_tot * (data.V_tot / data.V_CAV);
    PDE2_ECAV = data.f_PDE2_ECAV * data.PDE2_tot * (data.V_tot / data.V_ECAV);
    PDE2_CYT = (1 - data.f_PDE2_CAV - data.f_PDE2_ECAV) * data.PDE2_tot * (data.V_tot / data.V_CYT);

    PDE3_CAV = data.f_PDE3_CAV * data.PDE3_tot * (data.V_tot / data.V_CAV);         %Note f_PDE3_ECAV = 0;
    PDE3_CYT = (1 - data.f_PDE3_CAV) * data.PDE3_tot * (data.V_tot / data.V_CYT);

    PDE4_CAV = data.f_PDE4_CAV * data.PDE4_tot * (data.V_tot / data.V_CAV);
    PDE4_ECAV = data.f_PDE4_ECAV * data.PDE4_tot * (data.V_tot / data.V_ECAV);
    PDE4_CYT = (1 - data.f_PDE4_CAV - data.f_PDE4_ECAV) * data.PDE4_tot * (data.V_tot / data.V_CYT);

    dCAMP_PDEs = [Mod_PDE_Sub(y(10), PDE2_CAV         , 0, data.kPDE2, data.KmPDE2); ...
                  Mod_PDE_Sub(y(10), PDE3_CAV - y(34) , y(34), data.kPDE3, data.KmPDE3); ...
                  Mod_PDE_Sub(y(10), PDE4_CAV - y(36) , y(36), data.kPDE4, data.KmPDE4); ...
                  Mod_PDE_Sub(y(11), PDE2_ECAV        , 0, data.kPDE2, data.KmPDE2); ...
                  Mod_PDE_Sub(y(11), PDE4_ECAV - y(37), y(37), data.kPDE4, data.KmPDE4); ...
                  Mod_PDE_Sub(y(12), PDE2_CYT         , 0, data.kPDE2, data.KmPDE2); ...
                  Mod_PDE_Sub(y(12), PDE3_CYT - y(35) , y(35), data.kPDE3, data.KmPDE3); ...
                  Mod_PDE_Sub(y(12), PDE4_CYT - y(38) , y(38), data.kPDE4, data.KmPDE4); ...
                  ];
              
    dy(10) = dcAMP_PKA_CAV + data.rate_camp * (dcAMPAC56_CAV - (dCAMP_PDEs(1) + dCAMP_PDEs(2) + dCAMP_PDEs(3))) - data.J_CAV_ECAV * ((y(10) - y(11)) / data.V_CAV) - data.J_CAV_CYT * ((y(10) - y(12)) / data.V_CAV);
    dy(11) = dcAMP_PKA_ECAV + data.rate_camp * (dcAMPAC47_ECAV - (dCAMP_PDEs(4) + dCAMP_PDEs(5))) + data.J_CAV_ECAV * ((y(10) - y(11)) / data.V_ECAV) - data.J_ECAV_CYT * ((y(11) - y(12)) / data.V_ECAV);
    dy(12) = dcAMP_PKA_CYT + data.rate_camp * (dcAMPAC47_CYT + dcAMPAC56_CYT - (dCAMP_PDEs(6) + dCAMP_PDEs(7) + dCAMP_PDEs(8))) + data.J_CAV_CYT * ((y(10) - y(12)) / data.V_CYT) + data.J_ECAV_CYT * ((y(11) - y(12)) / data.V_CYT);
end

function dcAMP_PDEx_C = Mod_PDE_Sub(cAMP, PDE, PDEp, kPDE, KmPDE)
    dcAMP_PDEx_C = (kPDE * PDE * cAMP + 3.0 * kPDE * PDEp * cAMP) / (KmPDE + cAMP);
end

function [dy, dcAMP_PKA_CAV, dcAMP_PKA_ECAV, dcAMP_PKA_CYT] = Mod_PKA(y, dy, data)        
    PKA_CAV = data.f_PKA_CAV * data.PKA_tot * (data.V_tot / data.V_CAV);
    PKA_ECAV = data.f_PKA_ECAV * data.PKA_tot * (data.V_tot / data.V_ECAV);
    PKA_CYT = (1 - data.f_PKA_CAV - data.f_PKA_ECAV) * data.PKA_tot * (data.V_tot / data.V_CYT);
    
    PKI_CAV = data.f_PKI_CAV * data.PKI_tot * (data.V_tot / data.V_CAV);
    PKI_ECAV = data.f_PKI_ECAV * data.PKI_tot * (data.V_tot / data.V_ECAV);
    PKI_CYT = (1 - data.f_PKI_CAV - data.f_PKI_ECAV) * data.PKI_tot * (data.V_tot / data.V_CYT);
    
    kPKAb1 = data.k_PKA_f1 * data.K_PKAII_1;
    kPKAb2 = data.k_PKA_f2 * data.K_PKAII_2;
    kPKAb3 = data.k_PKA_f3 * data.K_PKAII_3;
    kPKAbpki = data.k_PKA_fPKI * data.K_PKA_PKI;
    
    % CAV Compartment    
    cAMP = y(10); LRC = y(19); L2RC = y(20); L2R = y(21); C = y(22); PKIC = y(23); RC = PKA_CAV - LRC - L2RC - L2R; PKI = PKI_CAV - PKIC;
    
    dcAMP_PKA_CAV = -data.k_PKA_f1 * RC * cAMP + kPKAb1 * LRC - data.k_PKA_f2 * LRC * cAMP + kPKAb2 * L2RC;
    dy(19) = data.k_PKA_f1 * RC * cAMP - kPKAb1 * LRC + kPKAb2 * L2RC - data.k_PKA_f2 * LRC * cAMP;      
    dy(20) = data.k_PKA_f2 * LRC * cAMP - kPKAb2 * L2RC - data.k_PKA_f3 * L2RC + kPKAb3 * L2R * C;       
    dy(21) = data.k_PKA_f3 * L2RC - kPKAb3 * L2R * C;
    dy(22) = data.k_PKA_f3 * L2RC - kPKAb3 * L2R * C + kPKAbpki * PKIC - data.k_PKA_fPKI * PKI * C;
    dy(23) = -kPKAbpki * PKIC + data.k_PKA_fPKI * PKI * C;
    
    % ECAV Compartment
    cAMP = y(11); LRC = y(24); L2RC = y(25); L2R = y(26); C = y(27); PKIC = y(28); RC = PKA_ECAV - LRC - L2RC - L2R; PKI = PKI_ECAV - PKIC;
    
    dcAMP_PKA_ECAV = -data.k_PKA_f1 * RC * cAMP + kPKAb1 * LRC - data.k_PKA_f2 * LRC * cAMP + kPKAb2 * L2RC;
    dy(24) = data.k_PKA_f1 * RC * cAMP - kPKAb1 * LRC + kPKAb2 * L2RC - data.k_PKA_f2 * LRC * cAMP;      
    dy(25) = data.k_PKA_f2 * LRC * cAMP - kPKAb2 * L2RC - data.k_PKA_f3 * L2RC + kPKAb3 * L2R * C;       
    dy(26) = data.k_PKA_f3 * L2RC - kPKAb3 * L2R * C;
    dy(27) = data.k_PKA_f3 * L2RC - kPKAb3 * L2R * C + kPKAbpki * PKIC - data.k_PKA_fPKI * PKI * C;
    dy(28) = -kPKAbpki * PKIC + data.k_PKA_fPKI * PKI * C;
    
    % CYT Compartment
    kPKAb1 = data.k_PKA_f1 * data.K_PKAI_1;
    kPKAb2 = data.k_PKA_f2 * data.K_PKAI_2;
    kPKAb3 = data.k_PKA_f3 * data.K_PKAI_3;
    kPKAbpki = data.k_PKA_fPKI * data.K_PKA_PKI;
    
    cAMP = y(12); LRC = y(29); L2RC = y(30); L2R = y(31); C = y(32); PKIC = y(33); RC = PKA_CYT - LRC - L2RC - L2R; PKI = PKI_CYT - PKIC;
    
    dcAMP_PKA_CYT = -data.k_PKA_f1 * RC * cAMP + kPKAb1 * LRC - data.k_PKA_f2 * LRC * cAMP + kPKAb2 * L2RC;
    dy(29) = data.k_PKA_f1 * RC * cAMP - kPKAb1 * LRC + kPKAb2 * L2RC - data.k_PKA_f2 * LRC * cAMP;      
    dy(30) = data.k_PKA_f2 * LRC * cAMP - kPKAb2 * L2RC - data.k_PKA_f3 * L2RC + kPKAb3 * L2R * C;       
    dy(31) = data.k_PKA_f3 * L2RC - kPKAb3 * L2R * C;
    dy(32) = data.k_PKA_f3 * L2RC - kPKAb3 * L2R * C + kPKAbpki * PKIC - data.k_PKA_fPKI * PKI * C;
    dy(33) = -kPKAbpki * PKIC + data.k_PKA_fPKI * PKI * C;
end

function dy = Mod_PDE_Phosphorylation(y, dy, data)
    PDE3_CAV = data.f_PDE3_CAV * data.PDE3_tot * (data.V_tot / data.V_CAV);         %Note f_PDE3_ECAV = 0;
    PDE3_CYT = (1 - data.f_PDE3_CAV) * data.PDE3_tot * (data.V_tot / data.V_CYT);
    
    PDE4_CAV = data.f_PDE4_CAV * data.PDE4_tot * (data.V_tot / data.V_CAV);
    PDE4_ECAV = data.f_PDE4_ECAV * data.PDE4_tot * (data.V_tot / data.V_ECAV);
    PDE4_CYT = (1 - data.f_PDE4_CAV - data.f_PDE4_ECAV) * data.PDE4_tot * (data.V_tot / data.V_CYT);
    
    data.kbPDEp = data.KPDEp * data.kfPDEp;
    dy(34) = data.kfPDEp * y(22) * (PDE3_CAV - y(34)) - data.kbPDEp * y(34);
    dy(35) = data.kfPDEp * y(32) * (PDE3_CYT - y(35)) - data.kbPDEp * y(35);
    dy(36) = data.kfPDEp * y(22) * (PDE4_CAV - y(36)) - data.kbPDEp * y(36);
    dy(37) = data.kfPDEp * y(27) * (PDE4_ECAV - y(37)) - data.kbPDEp * y(37);
    dy(38) = data.kfPDEp * y(32) * (PDE4_CYT - y(38)) - data.kbPDEp * y(38);
end

function [dy, PP1_CYT] = Mod_PP1_Inhibition(y, dy, data)
    I1_tot_CYT = data.f_PP1_Inh1 / (1 - data.f_PP1_Inh1) * data.KI1 + data.f_PP1_Inh1 * data.PP1_tot_CYT;
    dy(39) = (data.kfI1p * y(32) * (I1_tot_CYT - y(39))) / (data.KmI1p + (I1_tot_CYT - y(39))) - (data.kbI1p * data.PP2A_CYT * y(39)) / (data.KmI1dp + y(39));
    PP1_CYT = 0.5 * sqrt( (data.KI1 + y(39)- data.PP1_tot_CYT)^2 + 4 * data.PP1_tot_CYT * data.KI1) + 0.5 * data.PP1_tot_CYT - 0.5 * data.KI1 - 0.5 * y(39);
end

function [dy] = Mod_Channel_Phosphorylation(y, dy, PP1_CYT, data)       
    %Substrates without AKAP
    dy(42) = ( (data.k_PKAPLB * y(32) * (1 - y(42))) / (data.Km_PKAPLB + (1 - y(42))) - (data.k_PP1PLB * PP1_CYT * y(42)) / (data.Km_PP1PLB + y(42)) );
    dy(43) = (data.k_PKATnI * y(32) * (1 - y(43))) / (data.Km_PKATnI + (1 - y(43))) - (data.k_PPTnI * data.PP2A_CYT * y(43)) / (data.Km_PPTnI + y(43));
    dy(44) = ( (data.k_PKAINa * y(22) * (1 - y(44))) / (data.Km_PKAINa + (1 - y(44))) - (data.k_PP1INa * data.PP_CAV * y(44)) / (data.Km_PP1INa + y(44)));   
    dy(45) = (data.k_PKAINaK * y(22) * (1 - y(45))) / (data.Km_PKAINaK + (1 - y(45))) - (data.k_PP1INaK * data.PP_CAV * y(45)) / (data.Km_PP1INaK + y(45));
    dy(47) = data.temp_rate_IKur * ( (data.k_PKAIKur * y(27) * (1 - y(47))) / (data.Km_PKAIKur + (1 - y(47))) - (data.k_PPIKur * data.PP1_ECAV * y(47)) / (data.Km_PPIKur + y(47)) );
    
    %Substrates with AKAP
    dy(41) = ((data.k_PKAIKS * y(27) * (data.IKs_AKAP_PKA_PP1 - y(41))) / (data.Km_PKAIKS / data.AKAP_Factor + (data.IKs_AKAP_PKA_PP1 - y(41)))) ...
        - ((data.k_PP1IKS * data.PP1_ECAV * y(41)) / (data.Km_PP1IKS / data.AKAP_Factor + y(41)));
    
    
    dy(40) = data.temp_rate_LCC * ( ((data.k_PKALCC * y(22) * (data.ICaL_AKAP_PKA_PP - y(40))) / (data.Km_PKALCC / data.AKAP_Factor + (data.ICaL_AKAP_PKA_PP - y(40)))) ...
        - ((data.k_PPLCC * data.PP_CAV * y(40)) / (data.Km_PPLCC / data.AKAP_Factor + y(40))) );
    
    dy(46) = data.temp_rate_RyR * ( ((data.k_PKARyR * y(22) * (data.RyR_AKAP_PKA_PP - y(46))) / (data.Km_PKARyR / data.AKAP_Factor + (data.RyR_AKAP_PKA_PP - y(46)))) ...
        - ((data.k_PPRyR * data.PP_CAV * y(46)) / (data.Km_PPRyR / data.AKAP_Factor + y(46))) );
    
%     if y(40) < 0, dy(40) = 0.001; end  
%     if y(40) > ICaL_AKAP_PKA_PP, dy(40) = -0.001; end
%     if y(46) < 0, dy(46) = 0.001; end 
% %     if y(42) > 1, dy(42) = -0.001; end
end