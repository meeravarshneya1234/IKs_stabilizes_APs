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
% *   When called with a partial (potentially empty) settings structure,  *
% *   this function will supply the default values for all other          *
% *   parameters.                                                         *
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

function settings = setDefaultSettings(settings)
    if ~isfield(settings, 'applyVoltageClamp'), settings.applyVoltageClamp = 0; end
    if ~isfield(settings, 'vcRate'), settings.vcRate = -100; end
    if ~isfield(settings, 'runSignalingPathway'), settings.runSignalingPathway = 1; end
    if ~isfield(settings, 'runElectrophysiol'), settings.runElectrophysiol = 1; end
    if ~isfield(settings, 'APDRepLevel'), settings.APDRepLevel = 0.95; end        
    if ~isfield(settings, 'showProgress'), settings.showProgress = 1; end
    if ~isfield(settings, 'bcl'), settings.bcl = 1000; end
    if ~isfield(settings, 'freq'), settings.freq = 10; end
    if ~isfield(settings, 'LastBCL'), settings.LastBCL = settings.freq; end
    if ~isfield(settings, 'storeLast'), settings.storeLast = NaN; end    
    if ~isfield(settings, 'K_o'), settings.K_o = 5.4; end
    if ~isfield(settings, 'Ca_o'), settings.Ca_o = 1.8; end
    if ~isfield(settings, 'Na_o'), settings.Na_o = 140; end
    if ~isfield(settings, 'Cl_o'), settings.Cl_o = 100; end
    if ~isfield(settings, 'Istim'), settings.Istim = -80; end
    if ~isfield(settings, 'stimdur'), settings.stimdur = 0.5; end    
  
    if ~isfield(settings, 'IKsNPParams'), settings.IKsNPParams = [7.3990e-003,3.1196e-002,8.0019e-001,5.6992e-003,4.1520e-002,1.3489e+000,3.8839e-001,-1.5019e-001,6.0693e-001,9.0654e-002,-1.1157e-001,2.8330e-002,3.1124e-003,-5.1660e-002,1.5522e+000,2.7304e-003,4.4198e-004,-1.2022e+000,4.0173e-004,2.0873e-004,1.9561e-001]; end
    if ~isfield(settings, 'IKsPParams'), settings.IKsPParams =  [9.9415e-003,4.4809e-002,5.8172e-001,3.3201e-003,9.4217e-002,9.5364e-001,5.6356e-001,-1.7986e-001,5.8381e-001,6.5700e-002,-1.1899e-001,1.2406e-002,3.8525e-004,-6.4118e-002,7.7992e-001,4.6171e-003,2.3730e-004,-1.9742e+000,2.2652e-004,2.4689e-004,1.9561e-001]; end                                                                  
    if ~isfield(settings, 'LCCNPParams'), settings.LCCNPParams = [0.59,0.8,0.052,13,0.132,13.56,9.45,70,49.1,10.349,26.553,0.213,10.807,17.5,3,0.2474,13.825,6.3836,0.001,14.9186,1E-6,1.552E-4,1.1E-3,3.3696,1.2E-2]; end
    if ~isfield(settings, 'LCCPParams'), settings.LCCPParams = [0.59,0.8,0.052,13,0.132,-4.7980,7.5699,70,49.1,10.349,38.494,0.213,10.807,29.979,3.1775,0.1,32.5,18.0,0.0001,6.0,1E-6,2.579E-04,0.002,10,0.01]; end
    if ~isfield(settings, 'INaNPParams'), settings.INaNPParams = [2.15*8.25*1.1, 0.13, 10.66+16.7434, -11.1, 0.3, -2.535E-7, -0.1, 32, 0.135, 80+7, -6.8, 3.56, 0.079, 3.1E5, 0.35, -1.2714E5, 0.2444, -6.948E-5, -0.04391, 37.78, 0.311, 79.23, 0.1212, -0.01052, -0.1378, 40.14, -47.1-11.3729, -13.7299, -7]; end   
    if ~isfield(settings, 'RyRP_Amp'), settings.RyRP_Amp = 1.9925; end
    if ~isfield(settings, 'RyRP_Tau'), settings.RyRP_Tau = 0.5357; end
    if ~isfield(settings, 'TauTr'), settings.TauTr = 75; end
    if ~isfield(settings, 'observeICaLSS'), settings.observeICaLSS = 0; end
    
    if ~isfield(settings, 'ICaLB'), settings.ICaLB = 1; end
    if ~isfield(settings, 'IKsB'), settings.IKsB = 1; end    
    if ~isfield(settings, 'IKrB'), settings.IKrB = 1; end
    if ~isfield(settings, 'INaKB'), settings.INaKB =1; end
    if ~isfield(settings, 'INaCaB'), settings.INaCaB = 1; end
    if ~isfield(settings, 'IKpB'), settings.IKpB = 1; end
    if ~isfield(settings, 'IK1B'), settings.IK1B = 1; end
    if ~isfield(settings, 'INabB'), settings.INabB = 1; end
    if ~isfield(settings, 'ITo1B'), settings.ITo1B = 1; end
    if ~isfield(settings, 'ITo2B'), settings.ITo2B = 1; end
    if ~isfield(settings, 'INaB'), settings.INaB = 1; end
    if ~isfield(settings, 'INaLB'), settings.INaLB = 1; end
    if ~isfield(settings, 'IClB'), settings.IClB = 1; end
    if ~isfield(settings, 'IpCaB'), settings.IpCaB = 1; end
    if ~isfield(settings, 'CTKClB'), settings.CTKClB = 0.0; end
    if ~isfield(settings, 'CTNaClB'), settings.CTNaClB = 0.0; end
    if ~isfield(settings, 'ICabB'), settings.ICabB = 1; end
    if ~isfield(settings, 'IrelB'), settings.IrelB = 1; end
    if ~isfield(settings, 'IupB'), settings.IupB = 1; end
    if ~isfield(settings, 'ItrB'), settings.ItrB = 0.0; end
    if ~isfield(settings, 'IleakB'), settings.IleakB = 1; end
    if ~isfield(settings, 'IdiffB'), settings.IdiffB = 0.0; end    
    if ~isfield(settings, 'CAMKIIB'), settings.CAMKIIB = 0; end        
end