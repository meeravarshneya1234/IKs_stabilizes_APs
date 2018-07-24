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
% *   Returns a structure with all paramaters used in calculating         *
% *   BAR signaling state derrivatives.                                   *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   -                                                                   *
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

function data = constants_Signaling()      
    l = 0.01;       % Length of the cell (cm)
    a = 0.0011;     % Radius of the cell (cm)
    vcell = 1000*pi*a*a*l;     %   3.801e-5 uL   % Cell volume (uL)
    data.V_CAV = 0.02 * vcell;
    data.V_ECAV = 0.04 * vcell;
    data.V_CYT = 0.678 * vcell;
    data.V_tot = vcell; 
        
    data.khydr_Gs = 0.8;
    data.kreas_Gs = 1.21E3; 
    data.khydr_Gi = data.khydr_Gs;
    data.kreas_Gi = data.kreas_Gs;
    data.ATP = 5E3;
    data.KmATP = 315;
    
    data.AC47_hill_Gs = 1.0043;
    data.AC47_basal = 0.03135;
    data.AC47_Km_Gs = 0.031544;
    data.AC47_Vmax_Gi = 0;
    data.AC47_hill_Gi = 0.8921;
    data.AC47_Km_Gi = 0.053733;
    
    data.kPDE2 = 20;
    data.KmPDE2 = 50;
    data.kPDE3 = 2.5;
    data.KmPDE3 = 0.8;
    data.kPDE4 = 4.0;
    data.KmPDE4 = 1.4;

    data.k_PKA_f1 = 100; data.K_PKAII_1 = 2.4984; data.K_PKAI_1 = 0.1088;
    data.k_PKA_f2 = 100; data.K_PKAII_2 = 11.359; data.K_PKAI_2 = 0.4612;
    data.k_PKA_f3 = 100; data.K_PKAII_3 = 0.3755; data.K_PKAI_3 = 0.3755;
    data.k_PKA_fPKI = 50; data.K_PKA_PKI = 0.01 / 50;
                
    data.KI1 = 1E-3;
    data.PP1_tot_CYT = 0.89;
    data.I1_tot_CYT = 0.3;
    data.kfI1p = 1.0145e-002;
    data.KmI1p = 1.4690e-003;
    data.kbI1p = 3.5731e-003;
    data.PP2A_CYT = 1;
    data.KmI1dp = 1.9526e-005;

    data.PP_CAV = 0.25;
    data.PP1_ECAV = 0.1;   
    
    data.K_PP1_Yotiao = 0.01; % 10 nM affinity similar to PKA based onZakhary et al.
    data.K_PKA_Yotiao = 0.01;
    data.K_IKs_Yotiao = 1E-4; %1E-4; % From Saucerman et al. 2004 
    data.IKs_tot = 0.025; data.Yotiao_tot = data.IKs_tot;
    
    data.K_PP1_AKAP_ICaL = 0.01; % 10 nM affinity similar to PKA based onZakhary et al.
    data.K_PKA_AKAP_ICaL = 0.01;
    data.K_ICaL_AKAP = 1E-4; % From Saucerman et al. 2004 
    data.ICaL_tot = 0.025; data.AKAP_ICaL_tot = data.ICaL_tot;
    
    data.K_PP1_AKAP_RyR = 0.01; % 10 nM affinity similar to PKA based onZakhary et al.
    data.K_PKA_AKAP_RyR = 0.01;
    data.K_RyR_AKAP = 1E-4; % From Saucerman et al. 2004 
    data.RyR_tot = 5*0.025; data.AKAP_RyR_tot = data.RyR_tot;           
        
    data.J_CAV_ECAV = 5E-15 * 1E6;
    data.J_CAV_CYT = 7.5E-14 * 1E6;
    data.J_ECAV_CYT = 0.9E-14 * 1E6;
    data.AKAP_Factor = 1; 
    
    data.rate_bds = 0.35;
    data.rate_camp = 1.0;
    data.k_GsAct_b2 = 1.0;
        
    data.Gi_tot = 0.5;
    data.R_b2_tot = 0.15 * 0.025;
    data.R_b1_tot = 0.85 * 0.025;
    data.Gs_tot = 224 * data.R_b1_tot;
    data.f_Rb1_CAV = 8.1161e-002;
    data.f_Rb1_ECAV = 4.8744e-001;
    data.f_Rb2_CAV = 0.85;
    data.f_Gi_CAV = 0.85;

    data.K_b1_h = 0.062; 
    data.K_b1_l = 0.567; 
    data.K_b1_c = 2.4490e+000;

    data.K_b2_h = 0.012;
    data.K_b2_l = 1.053;
    data.K_b2_c = 1.8463e+000;
    data.K_b2_a = 1.6655e+000;
    data.K_b2_f = 0.1; 
    data.K_b2_n = 1.053;

    data.f_Gs_CAV = 1.1071e-003;
    data.f_Gs_ECAV = 5.6640e-001;
    data.kact1_Gs = 4.9054e+000;
    data.kact2_Gs = 2.5945e-001;
    data.kact1_Gi = 4;
    data.kact2_Gi = 0.05;
        
    data.AC_tot = 3 * data.R_b1_tot;
    data.f_AC56_AC47 = 1 / (1 + 0.35);
    data.f_AC56_CAV = 8.7459e-002; 
    data.f_AC47_ECAV = 1.6479e-001;

    data.AC56_basal = 3.7696e-002; 
    data.AC56_hill_Gs = 1.3574; 
    data.AC56_Km_Gs = 0.0852;
    data.AC56_V_GsGi = 0.8569;
    data.AC56_hill_GsGi = 0.6623;
    data.AC56_Km_GsGi = 0.4824;
    data.AC56_Km_Gi = 0.0465;

    data.AF56 = 4.1320e+001; 
    data.AF47 = 3.3757e+000; 

    data.f_CYT_PDE3 = 0.35; % 35% of cytosolic PDE is PDE Type 3
    data.f_PDE4_PART = 0.125; % 7 times more cytosolic PDE4 than particulate PDE4
    data.r_PDE3_PDE4 = 3.71; % ratio of PDE3 to PDe4 in particulate fraction is 78 : 21
    data.f_PART_PDE = 0.2;                

    data.PDE2_tot = 2.9268e-002;
    data.f_PDE2_CAV = 1.6957e-001;
    data.f_PDE2_ECAV = 2.1257e-004;
    data.f_PDE4_CAV = 1.2481e-001; 
    data.f_PDE3_ECAV = 0.0;

    data.k_b1_pkap = 0.0065;
    data.k_b1_grkp = 0.00133;
    data.K_b1_PKA = 1.5629e-001;
    data.k_b1_grkdp = 0.0009833;

    data.kfPDEp = 0.0196;
    data.KPDEp = 5.2218e-001; 
    data.GRK_CAV = 1.0; 
    data.GRK_ECAV = 1.0;
    data.GRK_CYT = 1.0; 

    data.PKA_tot = 0.5;
    data.f_PKA_CAV = 0.0388; 
    data.f_PKA_ECAV = 0.1;
    data.PKI_tot = 0.2 * data.PKA_tot;
    data.f_PKI_CAV = data.f_PKA_CAV;
    data.f_PKI_ECAV = data.f_PKA_ECAV;          
                    
    data.temp_rate_LCC = 1.0;               
    data.k_PKALCC = 5.1009e-004; 
    data.Km_PKALCC = 1.2702e-006;
    data.k_PPLCC = 6.9030e-004; 
    data.Km_PPLCC = 6.3064e-003;

    data.temp_rate_RyR = 1.0;         
    data.k_PKARyR = 2.5548e-003; 
    data.Km_PKARyR = 6.6298e-005;
    data.k_PPRyR = 3.8257e-003; 
    data.Km_PPRyR = 4.3003e-02; 

    data.k_PKAIKS = 1.6305e-001;
    data.Km_PKAIKS = 9.9794e-005;
    data.k_PP1IKS = 1.0542; 
    data.Km_PP1IKS = 1.1147e-004; 

    data.k_PKAPLB = 1.1348e-001;
    data.Km_PKAPLB = 9.8854e-004;
    data.k_PP1PLB = 4.8302e-001;
    data.Km_PP1PLB = 8.0737e-001;

    data.k_PKATnI = 1.0408e-001; 
    data.Km_PKATnI = 2.7143e-005; 
    data.k_PPTnI = 5.2633e-002; 
    data.Km_PPTnI = 2.6714e-001; 

    data.k_PKAINa = 1.3680e-002; 
    data.Km_PKAINa = 1.0988e-001; 
    data.k_PP1INa = 5.2811e-02; 
    data.Km_PP1INa = 7.8605e+00; 

    data.k_PKAINaK = 1.5265e-002; 
    data.Km_PKAINaK = 1.1001e-03; 
    data.k_PP1INaK = 9.2455e-02; 
    data.Km_PP1INaK = 5.7392e+00; 

    data.temp_rate_IKur = 1.0;
    data.k_PKAIKur =  6.9537e-002; 
    data.Km_PKAIKur = 2.7623e-001; 
    data.k_PPIKur = 3.1700e-001; 
    data.Km_PPIKur = 2.3310e-003; 
        
    data.PP1_tot_CYT = 0.2;
    data.f_PP1_Inh1 = 0.3;
        
    data.IBMX = 0;
    
    % ======================================
    % Calculate compartmental values based on fractions and whole cell
    % levels.
    % ======================================
    data.k_b1_dp = data.K_b1_PKA * data.k_b1_pkap;
    r_PDE3_CYT = (data.f_CYT_PDE3 / (1 - data.f_CYT_PDE3));
    f_PDE2_PART = data.f_PDE2_CAV + data.f_PDE2_ECAV;
    alpha_PDE3 = r_PDE3_CYT * (data.f_PDE4_PART * (1 + data.r_PDE3_PDE4 - data.r_PDE3_PDE4 * f_PDE2_PART - data.f_PART_PDE) + f_PDE2_PART * (data.f_PART_PDE - 1)) + data.r_PDE3_PDE4 * data.f_PDE4_PART * (data.f_PART_PDE - f_PDE2_PART);
    beta_PDE3 = data.f_PDE4_PART * (1 + data.r_PDE3_PDE4 + data.f_PART_PDE * (r_PDE3_CYT - data.r_PDE3_PDE4)) - data.f_PART_PDE * (1 + r_PDE3_CYT);
    data.PDE3_tot = (alpha_PDE3 / beta_PDE3) * data.PDE2_tot;
    data.PDE4_tot = ((data.f_PART_PDE - f_PDE2_PART) * data.PDE2_tot + data.f_PART_PDE * data.PDE3_tot) / ((1 + data.r_PDE3_PDE4) * data.f_PDE4_PART - data.f_PART_PDE);
    data.f_PDE3_CAV = (data.r_PDE3_PDE4 * data.f_PDE4_PART * data.PDE4_tot) / data.PDE3_tot;
    data.f_PDE4_ECAV = data.f_PDE4_PART - data.f_PDE4_CAV;

    data.PDE2_tot = data.PDE2_tot * (1 - (data.IBMX ^ 1.167 / (21.58 + data.IBMX ^ 1.167)));
    data.PDE3_tot = data.PDE3_tot * (1 - (data.IBMX ^ 0.7629 / (2.642 + data.IBMX ^ 0.7629)));
    data.PDE4_tot = data.PDE4_tot * (1 - (data.IBMX ^ 0.9024 / (11.89 + data.IBMX ^ 0.9024)));

    data.AC56_CAV = data.f_AC56_CAV * data.f_AC56_AC47 * data.AC_tot * (data.V_tot / data.V_CAV);
    data.AC56_CYT = (1- data.f_AC56_CAV) * data.f_AC56_AC47 * data.AC_tot * (data.V_tot / data.V_CYT);
    data.AC47_ECAV = data.f_AC47_ECAV * (1 - data.f_AC56_AC47) * data.AC_tot * (data.V_tot / data.V_ECAV);
    data.AC47_CYT = (1- data.f_AC47_ECAV) * (1 - data.f_AC56_AC47) * data.AC_tot * (data.V_tot / data.V_CYT);

    PKA_CAV = data.f_PKA_CAV * data.PKA_tot * (data.V_tot / data.V_CAV);
    PKA_ECAV = data.f_PKA_ECAV * data.PKA_tot * (data.V_tot / data.V_ECAV);
        
    %IKs: compartment with 1 AKAP
    PP1_ECAV_free = 0.5 * data.K_PP1_Yotiao * (-(1 + (1 / data.K_PP1_Yotiao) * (data.Yotiao_tot - data.PP1_ECAV)) + sqrt((1 + (1 / data.K_PP1_Yotiao) * (data.Yotiao_tot - data.PP1_ECAV))^2 + 4 * data.PP1_ECAV * (1 / data.K_PP1_Yotiao)));
    IKs_free = 0.5 * data.K_IKs_Yotiao * (-(1 + (1 / data.K_IKs_Yotiao) * (data.Yotiao_tot - data.IKs_tot)) + sqrt((1 + (1 / data.K_IKs_Yotiao) * (data.Yotiao_tot - data.IKs_tot))^2 + 4 * data.IKs_tot * (1 / data.K_IKs_Yotiao)));
    R_free = 0.5 * data.K_PKA_Yotiao * (-(1 + (1 / data.K_PKA_Yotiao) * (data.Yotiao_tot - PKA_ECAV)) + sqrt((1 + (1 / data.K_PKA_Yotiao) * (data.Yotiao_tot - PKA_ECAV))^2 + 4 * PKA_ECAV * (1 / data.K_PKA_Yotiao)));
    Yotiao_free = (data.Yotiao_tot  - data.IKs_tot + IKs_free) / ((1 / data.K_PP1_Yotiao) * PP1_ECAV_free + (1 / data.K_PKA_Yotiao) * ((1 / data.K_PP1_Yotiao) * PP1_ECAV_free + 1) * R_free + 1);
    
    data.IKs_AKAP_PKA = (1 / (data.K_IKs_Yotiao * data.K_PKA_Yotiao)) * IKs_free * Yotiao_free * R_free;
    data.IKs_AKAP_PKA_PP1 = data.IKs_AKAP_PKA * (1 / data.K_PP1_Yotiao) * PP1_ECAV_free;
    
    %ICaL and RyR: compartment with 2 AKAPs
    b_pp1 = data.AKAP_ICaL_tot + data.AKAP_RyR_tot + data.K_PP1_AKAP_ICaL + data.K_PP1_AKAP_RyR - data.PP_CAV;
    c_pp1 = data.AKAP_ICaL_tot * data.K_PP1_AKAP_RyR + data.AKAP_RyR_tot * data.K_PP1_AKAP_ICaL + data.K_PP1_AKAP_ICaL * data.K_PP1_AKAP_RyR - data.PP_CAV * (data.K_PP1_AKAP_ICaL + data.K_PP1_AKAP_RyR);
    d_pp1 = data.PP_CAV * data.K_PP1_AKAP_ICaL * data.K_PP1_AKAP_RyR;
    phi_pp1 = (0.5 * d_pp1 + (b_pp1 * c_pp1 / 6.0) - (b_pp1 * b_pp1 * b_pp1 / 27) + (-(b_pp1 ^ 3 * d_pp1 / 27) - (b_pp1 ^2 * c_pp1^2 / 108) + (b_pp1 * c_pp1 * d_pp1 / 6) + (c_pp1 ^ 3 / 27) + (d_pp1 ^ 2 / 4)) ^ 0.5) ^ (1 / 3.0);
    PP1_CAV_free = (phi_pp1 - (1 / phi_pp1) * ((c_pp1 / 3) - (b_pp1 * b_pp1 / 9)) - (b_pp1 / 3));
    PP1_CAV_free = real(PP1_CAV_free);
    
    b_pka = data.AKAP_ICaL_tot + data.AKAP_RyR_tot + data.K_PKA_AKAP_ICaL + data.K_PKA_AKAP_RyR - PKA_CAV;
    c_pka = data.AKAP_ICaL_tot * data.K_PKA_AKAP_RyR + data.AKAP_RyR_tot * data.K_PKA_AKAP_ICaL + data.K_PKA_AKAP_ICaL * data.K_PKA_AKAP_RyR - PKA_CAV * (data.K_PKA_AKAP_ICaL + data.K_PKA_AKAP_RyR);
    d_pka = PKA_CAV * data.K_PKA_AKAP_ICaL * data.K_PKA_AKAP_RyR;
    phi_pka = (0.5 * d_pka + (b_pka * c_pka / 6.0) - (b_pka * b_pka * b_pka / 27) + (-(b_pka ^ 3 * d_pka / 27) - (b_pka ^2 * c_pka^2 / 108) + (b_pka * c_pka * d_pka / 6) + (c_pka ^ 3 / 27) + (d_pka ^ 2 / 4)) ^ 0.5) ^ (1 / 3.0);
    R_CAV_free = (phi_pka - (1 / phi_pka) * ((c_pka / 3) - (b_pka * b_pka / 9)) - (b_pka / 3));
    R_CAV_free = real(R_CAV_free);
    
    ICaL_free = 0.5 * data.K_ICaL_AKAP * (-(1 + (1 / data.K_ICaL_AKAP) * (data.AKAP_ICaL_tot - data.ICaL_tot)) + sqrt((1 + (1 / data.K_ICaL_AKAP) * (data.AKAP_ICaL_tot - data.ICaL_tot))^2 + 4 * data.ICaL_tot * (1 / data.K_ICaL_AKAP)));
    RyR_free = 0.5 * data.K_RyR_AKAP * (-(1 + (1 / data.K_RyR_AKAP) * (data.AKAP_RyR_tot - data.RyR_tot)) + sqrt((1 + (1 / data.K_RyR_AKAP) * (data.AKAP_RyR_tot - data.RyR_tot))^2 + 4 * data.RyR_tot * (1 / data.K_RyR_AKAP)));

    AKAP_ICaL_free = (data.AKAP_ICaL_tot  - data.ICaL_tot + ICaL_free) / ((1 / data.K_PP1_AKAP_ICaL) * PP1_CAV_free + (1 / data.K_PKA_AKAP_ICaL) * ((1 / data.K_PP1_AKAP_ICaL) * PP1_CAV_free + 1) * R_CAV_free + 1);
    AKAP_RyR_free = (data.AKAP_RyR_tot  - data.RyR_tot + RyR_free) / ((1 / data.K_PP1_AKAP_RyR) * PP1_CAV_free + (1 / data.K_PKA_AKAP_RyR) * ((1 / data.K_PP1_AKAP_RyR) * PP1_CAV_free + 1) * R_CAV_free + 1);
    
	data.ICaL_AKAP_PKA = (1 / (data.K_ICaL_AKAP * data.K_PKA_AKAP_ICaL)) * ICaL_free * AKAP_ICaL_free * R_CAV_free;
    data.ICaL_AKAP_PKA_PP = data.ICaL_AKAP_PKA * (1 / data.K_PP1_AKAP_ICaL) * PP1_CAV_free;
        
    data.RyR_AKAP_PKA = (1 / (data.K_RyR_AKAP * data.K_PKA_AKAP_RyR)) * RyR_free * AKAP_RyR_free * R_CAV_free;
    data.RyR_AKAP_PKA_PP = data.RyR_AKAP_PKA * (1 / data.K_PP1_AKAP_RyR) * PP1_CAV_free;
end