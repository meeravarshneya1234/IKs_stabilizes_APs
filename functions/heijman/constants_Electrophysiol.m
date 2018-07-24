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
% *   electrophysiological state derrivatives.                            *
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

function data = constants_Electrophysiol()   
    % General parameters
    % ================================================
    data.F=96487;   % Faraday constant
    data.R=8314;    % Gas constant
    data.T=310;     % Temperature

    l = 0.01;       % Length of the cell (cm)
    a = 0.0011;     % Radius of the cell (cm)

    vcell = 1000*3.14*a*a*l;     %   3.801e-5 uL   Cell volume (uL)
    ageo = 2*3.14*a*a+2*pi*a*l;  %   7.671e-5 cm^2 Geometric membrane area (cm^2)
    Acap = ageo*2;               %   1.534e-4 cm^2 Capacitive membrane area (cm^2)

    data.vmyo = vcell*0.678;     % Myoplasm volume (uL)
    data.vsr = vcell*0.06;       % SR volume (uL)
    data.vnsr = vcell*0.0552;    % NSR volume (uL)  
    data.vjsr=  vcell*0.0048;    % JSR volume (uL)
    data.vss= vcell*0.02;        % Subspace volume (uL)
    data.vssCaL = vcell * 0.002; % L-type subspace volume (uL)
    
    data.AF=Acap/data.F;
    data.frt=data.F/data.T/data.R;
    data.RTF=1./data.frt;
    
    % Sodium
    % ================================================
    data.GNa=9.075;  
    data.GNaL=65e-4;

    data.KmCa=1.25e-4;     %% Na-Ca Exchanger Current
    data.NCXmax=4.5;
    data.ksat = 0.32;
    data.eta = 0.27;
    data.KmNai=12.3; data.KmNao=87.5;
    data.KmCai=0.0036; data.KmCao=1.3;

    % Calcium
    % ================================================
    data.ibarpca =0.0575; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
    data.kmpca = 5e-4;    % Half-saturation concentration of sarcolemmal Ca pump (mM)
    data.gcab=1.995084e-7;% ICab
    
    data.k_JSR_Leak = 1.75E-04;
    data.k_JSR_Leak_P = 5.0E-04;
    data.Km_JSR_Leak = 20;
    data.Km_JSR_Leak_P = 1.1;
    
    data.cmdnbar = 0.050;   % Buffer constants
    data.trpnbar = 0.070;   
    data.kmcmdn = 0.00238;  
    data.kmtrpn = 0.0005;
    data.kmcsqn=0.8;
    data.csqnbar=10;
    data.BSRmax=0.047; data.KmBSR=0.00087; 
    data.BSLmax=1.124; data.KmBSL=0.0087;

    data.KmCaMK=0.25;       % Iup
	data.iupmax=0.004375; 
	data.Kmup=0.00092; 
    data.nsrmax=15.0;
    
    data.dKmPLBCAMK = 0.00017;
    data.dKmPLBPKA = 0.46 * data.Kmup;
    data.dKmPLBBoth = data.dKmPLBPKA;    
    data.dJupmax=0.75;

    data.CaMK0=0.05;        % CAMKII signaling
    data.Km=1.5e-3;
    data.betaCamK=6.8e-4;
    data.alphaCamK=0.05;

    data.tau_PLB_CaMKII = 100000;
    data.tau_RyR_CaMKII = 10000;
    
    % Potassium
    % ================================================
    data.GK1max = 0.5;
    data.GKpmax = 3.84E-03; %2.76E-3; 
    data.GKrmax = 0.0138542;        
    data.prnak=0.01833;     % IKs sodium / potassium permeability ratio
       
    data.kmnai_NP = 2.6;    % Half-saturation concentration of NaK pump (mM)
    data.kmnai_P = 1.846;   % Half-saturation concentration of NaK pump (mM)
    data.kmko = 1.5;        % Half-saturation concentration of NaK pump (mM)
    data.ibarnak = 1.4;     % Max. current through Na-K pump (uA/uF)
    
    data.gitodv = 0.4975;   % ITo1

    % Chloride
    % ================================================
    data.Kmto2=0.1502;      % ITo2
    data.PCl=9e-7;
    data.GClb=2.25e-4;      % IClb

    data.CTKClmax=1.77E-5;
    data.CTNaClmax=2.46108E-5;   
end