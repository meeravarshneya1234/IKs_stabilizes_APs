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
% *   function will calculate dy/dt which can be used in an ODE solver.   *
% *   Calculations are performed in two subroutines; one for the          *
% *   adrenergic signaling pathway (state variables 89 and onwards), and  *
% *   one for the electrophysiology (state variables 1-88), connected     *
% *   using the fraction of phosphorylated channels.                      *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   constants_Signaling.m, dfun_Signaling.m, dfun_Electrophysiol.m      *
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

function dy = funHRdBA(t, y, flags, settings,block)
    if isempty(t) && isempty(y)
        if strcmp(flags, 'jpattern') == 1
            dy = settings.sparseJacobian;
            return;
        end        
    end
        
    dySignaling = zeros(57,1);
    if settings.runSignalingPathway == 1
        dySignaling = 0.001 * dfun_Signaling(t, y(89:end), flags, settings); % conversion from 1 / sec to 1 / ms

        sigdata = settings.dataSignaling;
        settings.fICaLP = block.ICaL*(min(max(((y(128) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1));
        settings.fIKsP = block.IKs*min(max(((y(129) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);    
        settings.fPLBP = block.PLB*min(max((y(130) - 6.591000e-001) / (9.945000e-001 - 6.591000e-001), 0), 1);    
        settings.fTnIP = block.TnI*min(max((y(131) - 6.735188e-001) / (9.991797e-001 - 6.735188e-001), 0), 1);    
        settings.fINaP = block.INa*min(max((y(132) - 2.394795e-001) / (9.501431e-001 - 2.394795e-001), 0), 1);    
        settings.fINaKP = block.INaK*min(max((y(133) - 1.263453e-001) / (9.980137e-001 - 1.263453e-001), 0), 1);    
        settings.fRyRP = block.RyR*min(max(((y(134) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);            
        settings.fIKurP = block.IKur*min(max((y(135) - 5.893798e-002) / (3.937470e-001 - 5.893798e-002), 0), 1);
    end
    
    dyElectrophysiol = zeros(88,1);
    if settings.runElectrophysiol == 1
        settings.Whole_cell_PP1 = 0.1371;
        if settings.runSignalingPathway == 1                 
            PP1_CYT = 0.5 * sqrt( (sigdata.KI1 + y(127) - sigdata.PP1_tot_CYT)^2 + 4 * sigdata.PP1_tot_CYT * sigdata.KI1) + 0.5 * sigdata.PP1_tot_CYT - 0.5 * sigdata.KI1 - 0.5 * y(127);
            settings.Whole_cell_PP1 = sigdata.V_CAV / sigdata.V_tot * sigdata.PP_CAV + sigdata.V_ECAV / sigdata.V_tot * sigdata.PP1_ECAV + sigdata.V_CYT / sigdata.V_tot * PP1_CYT;
        end
        dyElectrophysiol = dfun_Electrophysiol(t, y(1:88), flags, settings);
    end
    
    dy = [dyElectrophysiol; dySignaling];
end