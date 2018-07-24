% -- parameters for each of the models used in the study
function [p,c] = parameters(model_name,celltype)

%% FOX
if strcmp(model_name,'Fox')
    
    % Physical constants
    p.F = 96.5 ;                   % Faraday constant, coulombs/mmol
    p.R = 8.314 ;                  % gas constant, J/K
    p.T = 273+25 ;                 % absolute temperature, K
    p.RTF = p.R*p.T/p.F ;
    
    % Cell geometry constants
    p.Acap = 1.534e-4 ;            % cm^2
    p.Cm = 1 ;                     % uF/cm^2
    p.Vmyo = 25.84e-6 ;            % uL
    p.VSR = 2e-6 ;                 % uL
    
    % Fixed ionic concentrations
    % Initial conditions of others listed below
    p.Ko = 4000 ;                  % uM
    p.Nao = 138000 ;               % uM
    p.Cao = 2000 ;                 % uM
    
    p.Ki = 148400 ;                % uM
    p.Nai = 10000 ;                % uM
    
    c.GNa_ = 12.8 ;                  % ms/uF
    c.GK1_ = 2.8 ;                   % ms/uF
    c.GKr_ = 0.0136 ;                % ms/uF
    c.GKs_ = 0.0245 ;                % ms/uF
    c.GKp_ = 0.002216 ;               % ms/uF
    c.Gto_ = 0.23815 ;               % ms/uF
    c.GNab_ = 0.0031 ;               % ms/uF
    c.GCab_ = 3.842e-4 ;             % ms/uF
    
    c.PCa_ = 2.26e-5 ;                % cm/ms
    c.PCaK = 5.79e-7 ;               % cm/ms
    c.Prel = 6 ;                     % ms-1
    c.Pleak = 1e-6 ;                 % ms-1
    
    c.INaK_ = 0.693 ;                % uA/uF
    p.ICahalf = -0.265 ;             % uA/uF
    c.IpCa_ = 0.05 ;                 % uA/uF
    c.kNaCa = 1500 ;                 % uA/uF
    c.Vup = 0.1 ;                    % uM/ms
    
    % % values needed for calculation of NCX current
    p.eta = 0.35 ;
    p.ksat = 0.2 ;
    
    % Several saturation constants, nomenclature gets tricky
    % NCX
    p.KmNa = 87500 ;                 % uM
    p.KmCa = 1380 ;                  % uM
    % Ca-inactivation of L-type current
    p.KmfCa = 0.18 ;                 % uM
    % IK1
    p.KmK1 = 13000 ;                 % uM
    % Na-K pump
    p.KmNai = 10000 ;                % uM
    p.KmKo = 1500 ;                  % uM
    % sarcolemmal Ca pump
    p.KmpCa = 0.05 ;                 % uM
    % SERCA
    p.Kmup = 0.32 ;                  % uM
    
    % Buffering parameters
    p.CMDNtot = 10 ;                 % uM
    p.KmCMDN = 2 ;                   % uM
    p.CSQNtot = 10000 ;              % uM
    p.KmCSQN = 600 ;                 % uM
    %% HUND
elseif strcmp(model_name,'Hund')
    % Physical constants
    p.F = 96485 ;                  % Faraday's constant, C/mol
    p.R = 8314 ;                   % gas constant, mJ/K
    p.T = 273+25 ;                 % absolute temperature, K
    p.RTF = p.R*p.T/p.F ;
    
    % Cell geometry constants
    p.Acap = 1.534e-4 ;            % cm^2
    p.Vmyo = 25.84e-6 ;            % uL
    p.Vss = 0.76e-6 ;              % uL
    p.VJSR = 0.182e-6 ;            % uL
    p.VNSR = 2.098e-6 ;            % uL
    
    % Fixed ionic concentrations
    % Initial conditions of others listed below
    p.Ko = 5.4 ;                  % uM
    p.Nao = 140 ;               % uM
    p.Clo = 100 ;
    p.Cao = 1.8 ;                 % uM
    
    % % New global variables that will be defined in 'hund_ode.m
    c.GNa_= 8.25 ;                          % mS/cm^2
    c.GNaL_= 6.5e-3 ;
    c.PCa_ = 2.43e-4 ;                       % cm/s
    c.PCab = 1.995084e-7 ;                  % cm/s
    p.gamma_Cao = 0.341 ;                   % dimensionless
    p.gamma_Cai = 1 ;                       % dimensionless
    p.dtaufCa_max = 10 ;                    % ms
    c.GK1_ = 0.5 ;
    c.GKr_ = 0.0138542 ;
    c.GKs_ = 0.0248975 ;
    p.pKNa = 0.01833 ;                  % relative permeability of IKs, Na to K
    c.GKp_ = 2.76e-3 ;
    c.Gto_ = 0.19 ;
    p.Kmto2 = 0.1502 ;
    c.GClb_ = 2.25e-4 ;
    p.CTKCl_ = 7.0756e-6 ;
    p.CTNaCl_ = 9.8443e-6 ;
    c.PCl = 4e-7 ;
    p.KmNa_NaK = 10 ;             % Half-saturation concentration of NaK pump (mM)
    p.KmK_NaK = 1.5 ;             % Half-saturation concentration of NaK pump (mM)
    c.INaK_ = 0.61875 ;           % Max. current through Na-K pump (uA/uF)
    c.kNaCa = 4.5 ;
    p.ksat = 0.27 ;
    p.eta = 0.35 ;
    p.KmNai = 12.3 ;
    p.KmNao = 87.5 ;
    p.KmCai = 0.0036 ;
    p.KmCao = 1.3 ;
    p.KmCa_allo = 1.25e-4 ;
    c.Krel_ = 3000 ;
    % This Km is same for CaMK effects on release, ICa gating, and uptake
    p.KmCaMK = 0.15 ;
    p.dtau_rel_max = 10 ;
    c.IpCa_ = 0.0575 ; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
    p.KmpCa = 0.5e-3 ; % Half-saturation concentration of sarcolemmal Ca pump (mM)
    c.Vup = 0.004375 ;
    p.Kmup = 0.00092 ;
    p.CaNSR_max = 15.0 ;
    p.tau_diff = 0.2 ;
    p.tau_transfer = 120 ;
    p.dKmPLBmax = 0.00017 ;
    p.dJupmax = 0.75 ;
    
    % % Buffers in cytosol
    p.TRPNtot = 70e-3 ;
    p.KmTRPN = 0.5e-3 ;
    p.CMDNtot = 50e-3 ;
    p.KmCMDN = 2.38e-3 ;
    
    % % Buffers in JSR
    p.CSQNtot = 10 ;
    p.KmCSQN = 0.8 ;
    
    % % Buffers in subspace
    p.BSRtot = 0.047 ;
    p.KmBSR = 0.00087 ;
    p.BSLtot = 1.124 ;
    p.KmBSL = 0.0087 ;
    p.CaMK0 = 0.05 ;
    p.KmCaM = 0.0015 ;
    %% LIVSHITZ
elseif strcmp(model_name,'Livshitz')
    % Physical constants
    p.F = 96485 ;                  % Faraday's constant, C/mol
    p.R = 8314 ;                   % gas constant, mJ/K
    p.T = 273+37 ;                 % absolute temperature, K
    p.RTF = p.R*p.T/p.F ;
    p.FRT = 1/p.RTF ;
    
    % Cell Geometry
    p.length_cell = 0.01;       % Length of the cell (cm)
    p.radius = 0.0011;     % Radius of the cell (cm)
    p.Vcell = 1000*pi*p.radius*p.radius*p.length_cell ;     %   3.801e-5 uL   % Cell volume (uL)
    p.Ageo = 2*pi*p.radius*p.radius+2*pi*p.radius*p.length_cell ;  %   7.671e-5 cm^2    % Geometric membrane area (cm^2)
    p.Acap = 2*p.Ageo ;             %   1.534e-4 cm^2    % Capacitive membrane area (cm^2)
    p.Vmyo = p.Vcell*0.68 ;    % Myoplasm volume (uL)
    p.Vmito = p.Vcell*0.24 ;  % Mitochondria volume (uL)
    % data.vsr = Vcell*0.06;    % SR volume (uL)
    p.VNSR = p.Vcell*0.0552 ;   % NSR volume (uL)
    p.VJSR =  p.Vcell*0.0048 ;   % JSR volume (uL)
    p.Vss = p.Vcell*0.02 ;
    
    % Fixed ionic concentrations
    % Initial conditions of others listed below
    p.Ko = 4.5 ;                  % uM
    p.Nao = 140 ;               % uM
    p.Cao = 1.8 ;                 % uM
    
    % % Na current
    c.GNa_= 16 ;                          % mS/cm^2
    c.GNab = 0.004 ;
    % GNaL_= 6.5e-3 ;
    
    % % Ca current
    c.PCa_ = 5.4e-4 ;                        % cm/s
    c.PCa_Na = 6.75e-7 ;                    % cm/s
    c.PCa_K = 1.93e-7 ;                     % cm/s
    c.PCab = 1.995084e-7 ;                  % cm/s
    p.gamma_Cao = 0.341 ;                   % dimensionless
    p.gamma_Cai = 1 ;                       % dimensionless
    p.gamma_Nao = 0.75 ;                    % dimensionless
    p.gamma_Nai = 0.75 ;                    % dimensionless
    p.gamma_Ko = 0.75 ;                     % dimensionless
    p.gamma_Ki = 0.75 ;                     % dimensionless
    % hLca = 1 ;                            % dimensionless, Hill coefficient
    % % Need to make sure this variable not used elsewhere
    p.KmCa = 6e-4 ;                         % Half saturation constant, mM
    % % T-type & background currents
    c.GCaT = 0.05 ;
    c.GCab = 0.003016 ;
    
    c.GK1_ = 0.75;
    c.GKr_ = 0.02614 ;
    c.GKs_ = 0.433 ;
    p.pKNa = 0.01833 ;                  % relative permeability of IKs, Na to K
    c.GKp_ = 5.52e-3 ;
    
    c.INaK_ = 2.25 ;             % Max. current through Na-K pump (uA/uF)
    p.KmNa_NaK = 10 ;             % Half-saturation concentration of NaK pump (mM)
    p.KmK_NaK = 1.5 ;             % Half-saturation concentration of NaK pump (mM)
    
    c.kNCX = 0.00025 ;
    p.ksat = 0.0001 ;
    p.eta = 0.15 ;
    
    c.alpha_rel = 0.125 ;
    p.Krel_inf = 1 ;
    p.hrel = 9 ;
    p.beta_tau = 4.75 ;
    p.Krel_tau = 0.0123 ;
    
    c.IpCa_ = 1.15 ; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
    p.KmpCa = 5e-4 ; % Half-saturation concentration of sarcolemmal Ca pump (mM)
    
    c.Vserca = 8.75e-3 ;              % mM/ms
    p.Kmserca = 9.0e-4 ;              % mM
    p.CaNSR_max = 15.0 ;
    p.tau_transfer = 120 ;
    
    % % Buffers in cytosol
    p.TRPNtot = 70e-3 ;
    p.KmTRPN = 0.5e-3 ;
    
    % % % % troponin buffering might be more complicated now, need to check this
    % trpnf = 40;   % forward  buffered in TRPN (mM)
    % trpnb = 0.02;   % backward  TRPN (mM)
    % cmdnf = 100;   % forward  buffered in TRPN (mM)
    % cmdnb = 0.238;   % backward  TRPN (mM)
    
    p.CMDNtot = 50e-3 ;
    p.KmCMDN = 2.38e-3 ;
    
    % % Buffers in JSR
    p.CSQNtot = 10 ;
    p.KmCSQN = 0.8 ;
    
    %% DEVENYI
elseif strcmp(model_name,'Devenyi')
    % Physical constants
    p.F = 96485 ;                  % Faraday's constant, C/mol
    p.R = 8314 ;                   % gas constant, mJ/K
    p.T = 273+37 ;                 % absolute temperature, K
    p.RTF = p.R*p.T/p.F ;
    p.FRT = 1/p.RTF ;
    
    p.length_cell = 0.01;       % Length of the cell (cm)
    p.radius = 0.0011;     % Radius of the cell (cm)
    p.Vcell = 1000*pi*p.radius*p.radius*p.length_cell ;     %   3.801e-5 uL   % Cell volume (uL)
    p.Ageo = 2*pi*p.radius*p.radius+2*pi*p.radius*p.length_cell ;  %   7.671e-5 cm^2    % Geometric membrane area (cm^2)
    p.Acap = 2*p.Ageo ;             %   1.534e-4 cm^2    % Capacitive membrane area (cm^2)
    p.Vmyo = p.Vcell*0.68 ;    % Myoplasm volume (uL)
    % Vmito = Vcell*0.24 ;  % Mitochondria volume (uL)
    % data.vsr = Vcell*0.06;    % SR volume (uL)
    p.VNSR = p.Vcell*0.0552 ;   % NSR volume (uL)
    p.VJSR =  p.Vcell*0.0048 ;   % JSR volume (uL)
    p.Vss = p.Vcell*0.02 ;
    
    % Fixed ionic concentrations
    % Initial conditions of others listed below
    %%%Concentrations from Francis's experiments
    % Ko = 5.4;
    % Nao = 137;
    % Cao = 2.0;
    p.Ko = 4.5 ;                  % uM
    p.Nao = 140 ;               % uM
    p.Cao = 1.8 ;
    
    % % Na current
    c.GNa_= 16 ;                          % mS/cm^2
    c.GNab_ = 0.004 ;
    % GNaL_= 6.5e-3 ;
    
    % % Ca current
    c.PCa_ = 5.4e-4 ;                        % cm/s
    c.PCa_Na = 6.75e-7 ;                    % cm/s
    c.PCa_K = 1.93e-7 ;                     % cm/s
    % PCab = 1.995084e-7 ;                  % cm/s
    p.KCaL = 1; %Scale factor applied to all of PCa, PCa_Na, and PCa_K
    p.gamma_Cao = 0.341 ;                   % dimensionless
    p.gamma_Cai = 1 ;                       % dimensionless
    p.gamma_Nao = 0.75 ;                    % dimensionless
    p.gamma_Nai = 0.75 ;                    % dimensionless
    p.gamma_Ko = 0.75 ;                     % dimensionless
    p.gamma_Ki = 0.75 ;                     % dimensionless
    p.KmCa = 6e-4 ;                         % Half saturation constant, mM
    
    % % T-type & background currents
    c.CaT = 0.05 ;
    c.Cab = 0.003016 ;
    
    c.GK1_ = 0.75;
    c.GKr_ = 0.02614 ;
    c.GKs_ = 0.433 ;
    p.pKNa = 0.01833 ;                  % relative permeability of IKs, Na to K
    c.GKp_ = 5.52e-3 ;
    
    c.INaK_max = 2.25 ;             % Max. current through Na-K pump (uA/uF)
    p.KmNa_NaK = 10 ;             % Half-saturation concentration of NaK pump (mM)
    p.KmK_NaK = 1.5 ;             % Half-saturation concentration of NaK pump (mM)
    
    c.kNCX = 0.00025 ;
    p.ksat = 0.0001 ;
    p.eta = 0.15 ;
    
    c.alpha_rel = 0.125 ;
    p.Krel_inf = 1 ;
    p.hrel = 9 ;
    p.beta_tau = 4.75 ;
    p.Krel_tau = 0.0123 ;
    
    c.IpCa_max = 1.15 ; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
    p.KmpCa = 5e-4 ; % Half-saturation concentration of sarcolemmal Ca pump (mM)
    
    c.Vserca = 8.75e-3 ;              % mM/ms
    p.Kmserca = 9.0e-4 ;              % mM
    p.CaNSR_max = 15.0 ;
    p.tau_transfer = 120 ;
    
    p.TRPNtot = 70e-3 ;
    p.KmTRPN = 0.5e-3 ;
    p.CMDNtot = 50e-3 ;
    p.KmCMDN = 2.38e-3 ;
    
    % % Buffers in JSR
    p.CSQNtot = 10 ;
    p.KmCSQN = 0.8 ;
    
    % New Parameterization based on the Livshitz Model
    scaling = [...
        0.9667548;... G.Na_
        1.9912335;... G.Nab
        0.5249614;... G.KCaL
        0.3216265;... G.CaT
        1.5588248;... G.Cab
        1.0728723;... G.K1_
        1.7822444;... G.Kr_
        0.0423702;... G.Ks_
        0.0589097;... G.Kp_
        1.6446156;... G.INaK_max
        1.0957701;... G.kNCX
        0.2702221;... G.IpCa_max
        2.1396191]; %G.Vserca
    
    c.GNa_ = c.GNa_ * scaling(1);
    c.GNab_ = c.GNab_ * scaling(2);
    p.KCaL = p.KCaL * scaling(3);
    c.CaT = c.CaT * scaling(4);
    c.Cab = c.Cab * scaling(5);
    c.GK1_ = c.GK1_ * scaling(6);
    c.GKr_ = c.GKr_ * scaling(7);
    c.GKs_ = c.GKs_ * scaling(8);
    c.GKp_ = c.GKp_ * scaling(9);
    c.INaK_max = c.INaK_max * scaling(10);
    c.kNCX = c.kNCX * scaling(11);
    c.IpCa_max = c.IpCa_max * scaling(12);
    c.Vserca = c.Vserca*scaling(13);
    %% SHANNON
elseif strcmp(model_name,'Shannon')
    % % Physical parameters
    p.F = 96486.7;   % coulomb_per_mole (in model_parameters)
    p.R2 = 8314.3;   % joule_per_kilomole_kelvin (R in model_parameters)
    p.T = 310.0;   % kelvin (in model_parameters)
    
    % % Cell geometry
    p.Cm_per_area = 2.0e-6;   % farad_per_cm2 (in model_parameters)
    p.cell_length = 100.0;   % micrometre (in model_parameters)
    p.cell_radius = 10.25;   % micrometre (in model_parameters)
    p.Fx_Ks_SL = 0.89;   % dimensionless (in IKs)
    p.Fx_Ks_jct = 0.11;   % dimensionless (in IKs)
    p.pKNa = 0.01833;   % dimensionless (in IKs)
    
    % % Geometric and diffusional parameters
    p.A_SL_cytosol = 1.3e-4;   % cm2 (in ion_diffusion)
    p.A_jct_SL = 3.01e-6;   % cm2 (in ion_diffusion)
    p.D_Ca_SL_cytosol = 1.22e-6;   % dm2_per_second (in ion_diffusion)
    p.D_Ca_jct_SL = 1.64e-6;   % dm2_per_second (in ion_diffusion)
    p.D_Na_SL_cytosol = 1.79e-5;   % dm2_per_second (in ion_diffusion)
    p.D_Na_jct_SL = 1.09e-5;   % dm2_per_second (in ion_diffusion)
    p.x_SL_cytosol = 0.45;   % micrometre (in ion_diffusion)
    p.x_jct_SL = 0.5;   % micrometre (in ion_diffusion)
    
    % % Ionic concentrations
    p.Cao = 1.8;   % millimolar (in model_parameters)
    p.Cli = 15.0;   % millimolar (in model_parameters)
    p.Clo = 150.0;   % millimolar (in model_parameters)
    p.Ki = 135.0;   % millimolar (in model_parameters)
    p.Ko = 5.4;   % millimolar (in model_parameters)
    p.Mgi = 1.0;   % millimolar (in model_parameters)
    p.Nao = 140.0;   % millimolar (in model_parameters)
    
    p.Fx_NaBk_SL = 0.89;   % dimensionless (in INab)
    p.Fx_NaBk_jct = 0.11;   % dimensionless (in INab)
    c.G_NaBk = 0.297e-3;   % milliS_per_microF (in INab)
    p.Fx_Na_SL = 0.89;   % dimensionless (in INa)
    p.Fx_Na_jct = 0.11;   % dimensionless (in INa)
    c.G_INa = 16.0;   % milliS_per_microF (in INa)
    c.G_tof = 0.02;   % milliS_per_microF (in Itof)
    c.G_tos = 0.06;   % milliS_per_microF (in Itos)
    
    p.Fx_Cl_SL = 0.89;   % dimensionless (in ICl_Ca)
    p.Fx_Cl_jct = 0.11;   % dimensionless (in ICl_Ca)
    c.G_Cl = 0.109625;   % milliS_per_microF (in ICl_Ca)
    p.Kd_ClCa = 0.1;   % millimolar (in ICl_Ca)
    c.G_ClBk = 0.009;   % milliS_per_microF (in IClb)
    
    % % ICaL parameters
    p.Fx_ICaL_SL = 0.1;   % dimensionless (in ICaL)
    p.Fx_ICaL_jct = 0.9;   % dimensionless (in ICaL)
    c.PCa_ = 5.4e-4;   % litre_per_farad_millisecond (in ICaL)
    p.PK = 2.7e-7;   % litre_per_farad_millisecond (in ICaL)
    p.PNa = 1.5e-8;   % litre_per_farad_millisecond (in ICaL)
    p.Q10_CaL = 1.8;   % dimensionless (in ICaL)
    p.gamma_Cai = 0.341;   % dimensionless (in ICaL)
    p.gamma_Cao = 0.341;   % dimensionless (in ICaL)
    p.gamma_Ki = 0.75;   % dimensionless (in ICaL)
    p.gamma_Ko = 0.75;   % dimensionless (in ICaL)
    p.gamma_Nai = 0.75;   % dimensionless (in ICaL)
    p.gamma_Nao = 0.75;   % dimensionless (in ICaL)
    
    c.GKr_ = 0.03 ;
    c.GKs_ = 0.07 ;
    c.GK1_ = 0.9 ;
    
    % % SR Ca release parameters
    c.KSRleak = 5.348e-6;   % per_millisecond (in Jleak_SR)
    p.H2 = 1.787;   % dimensionless (H in Jpump_SR)
    p.Kmf = 0.000246;   % millimolar (in Jpump_SR)
    p.Kmr = 1.7;   % millimolar (in Jpump_SR)
    p.Q10_SRCaP = 2.6;   % dimensionless (in Jpump_SR)
    c.V_max_3 = 286.0e-6;   % millimolar_per_millisecond (V_max in Jpump_SR)
    p.EC50_SR = 0.45;   % millimolar (in Jrel_SR)
    p.HSR = 2.5;   % dimensionless (in Jrel_SR)
    p.Max_SR = 15.0;   % dimensionless (in Jrel_SR)
    p.Min_SR = 1.0;   % dimensionless (in Jrel_SR)
    p.kiCa = 0.5;   % per_millimolar_per_millisecond (in Jrel_SR)
    p.kim = 0.005;   % per_millisecond (in Jrel_SR)
    p.koCa = 10.0;   % per_millimolar2_per_millisecond (in Jrel_SR)
    p.kom = 0.06;   % per_millisecond (in Jrel_SR)
    c.ks = 25.0;   % per_millisecond (in Jrel_SR)
    
    % % Na-K pump paramters
    p.Fx_NaK_SL = 0.89;   % dimensionless (in INaK)
    p.Fx_NaK_jct = 0.11;   % dimensionless (in INaK)
    p.H_NaK = 4.0;   % dimensionless (in INaK)
    c.I_NaK_max = 1.91;   % microA_per_microF (in INaK)
    p.Km_Ko = 1.5;   % millimolar (in INaK)
    p.Km_Nai = 11.0;   % millimolar (in INaK)
    p.Q10_Km_Nai = 1.49;   % dimensionless (in INaK)
    p.Q10_NaK = 1.63;   % dimensionless (in INaK)
    
    % % NCX parameters
    p.Fx_NCX_SL = 0.89;   % dimensionless (in INaCa)
    p.Fx_NCX_jct = 0.11;   % dimensionless (in INaCa)
    p.HNa = 3.0;   % dimensionless (in INaCa)
    p.K_mCai = 0.00359;   % millimolar (in INaCa)
    p.K_mCao = 1.3;   % millimolar (in INaCa)
    p.K_mNai = 12.29;   % millimolar (in INaCa)
    p.K_mNao = 87.5;   % millimolar (in INaCa)
    p.Kd_act = 0.000256;   % millimolar (in INaCa)
    p.Q10_NCX = 1.57;   % dimensionless (in INaCa)
    c.V_max_2 = 9.0;   % microA_per_microF (V_max in INaCa)
    p.eta = 0.35;   % dimensionless (in INaCa)
    c.ksat = 0.27;   % dimensionless (in INaCa)
    
    % % SL Ca pump and background currents
    p.Fx_CaBk_SL = 0.89;   % dimensionless (in ICab)
    p.Fx_CaBk_jct = 0.11;   % dimensionless (in ICab)
    c.G_CaBk = 0.0002513;   % milliS_per_microF (in ICab)
    p.Fx_SLCaP_SL = 0.89;   % dimensionless (in ICap)
    p.Fx_SLCaP_jct = 0.11;   % dimensionless (in ICap)
    p.H1 = 1.6;   % dimensionless (H in ICap)
    p.Km = 0.0005;   % millimolar (in ICap)
    p.Q10_SLCaP = 2.35;   % dimensionless (in ICap)
    c.V_maxAF = 0.0673;   % microA_per_microF (in ICap)
    
    % % Buffering parameters
    p.Bmax_Calsequestrin = 0.14;   % millimolar (in Ca_buffer)
    p.Bmax_SLB_SL = 0.0374;   % millimolar (in Ca_buffer)
    p.Bmax_SLB_jct = 0.0046;   % millimolar (in Ca_buffer)
    p.Bmax_SLHigh_SL = 0.0134;   % millimolar (in Ca_buffer)
    p.Bmax_SLHigh_jct = 0.00165;   % millimolar (in Ca_buffer)
    p.koff_Calsequestrin = 65.0;   % per_millisecond (in Ca_buffer)
    p.koff_SLB = 1.3;   % per_millisecond (in Ca_buffer)
    p.koff_SLHigh = 30.0e-3;   % per_millisecond (in Ca_buffer)
    p.kon_Calsequestrin = 100.0;   % per_millimolar_per_millisecond (in Ca_buffer)
    p.kon_SL = 100.0;   % per_millimolar_per_millisecond (in Ca_buffer)
    
    p.Bmax_SL = 1.65;   % millimolar (in Na_buffer)
    p.Bmax_jct = 7.561;   % millimolar (in Na_buffer)
    p.koff = 1.0e-3;   % per_millisecond (in Na_buffer)
    p.kon = 0.0001;   % per_millimolar_per_millisecond (in Na_buffer)
    p.Bmax_Calmodulin = 0.024;   % millimolar (in cytosolic_Ca_buffer)
    p.Bmax_Myosin_Ca = 0.14;   % millimolar (in cytosolic_Ca_buffer)
    p.Bmax_Myosin_Mg = 0.14;   % millimolar (in cytosolic_Ca_buffer)
    p.Bmax_SRB = 0.0171;   % millimolar (in cytosolic_Ca_buffer)
    p.Bmax_TroponinC = 0.07;   % millimolar (in cytosolic_Ca_buffer)
    p.Bmax_TroponinC_Ca_Mg_Ca = 0.14;   % millimolar (in cytosolic_Ca_buffer)
    p.Bmax_TroponinC_Ca_Mg_Mg = 0.14;   % millimolar (in cytosolic_Ca_buffer)
    p.koff_Calmodulin = 238.0e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.koff_Myosin_Ca = 0.46e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.koff_Myosin_Mg = 0.057e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.koff_SRB = 60.0e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.koff_TroponinC = 19.6e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.koff_TroponinC_Ca_Mg_Ca = 0.032e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.koff_TroponinC_Ca_Mg_Mg = 3.33e-3;   % per_millisecond (in cytosolic_Ca_buffer)
    p.kon_Calmodulin = 34.0;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
    p.kon_Myosin_Ca = 13.8;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
    p.kon_Myosin_Mg = 15.7e-3;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
    p.kon_SRB = 100.0;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
    p.kon_TroponinC = 32.7;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
    p.kon_TroponinC_Ca_Mg_Ca = 2.37;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
    p.kon_TroponinC_Ca_Mg_Mg = 3.0e-3;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)

    %% TT04
elseif strcmp(model_name,'TT04')
    p.celltype = celltype;
    % Physical constants
    p.F = 96.4867 ;                % Faraday constant, coulombs/mmol
    p.R = 8.314 ;                  % gas constant, J/(K mol)
    p.T = 273+37 ;                 % absolute temperature, K
    p.RTF = p.R*p.T/p.F ;
    
    p.Cm = 1 ;                      % uF/cm2
    p.S = 0.2 ;                     % surface area to volume ratio
    
    p.Vmyo = 16.404 ;               % pL
    p.VSR = 1.094 ;                 % pL
    
    p.Ko = 5400 ;                   % uM
    p.Nao = 140000 ;                % uM
    p.Cao = 2000 ;                  % uM
    
    c.GNa_ = 14.838 ;               % nS/pF
    c.PCa_ = 1.75e-4 ;              % cm3 uF-1 s-1
    c.GK1_ = 5.405 ;                % nS/pF
    
    c.GpK_ = 0.0146 ;               % nS/pF
    c.GpCa_ = 0.025 ;               % nS/pF
    c.GNab_ = 2.9e-4 ;              % nS/pF
    c.GCab_ = 5.92e-4 ;             % nS/pF
    c.GKr_ = 0.096 ;                % nS/pF
    p.pKNa = 0.03 ;                  % relative permeability, Na to K
    
    % Maximum rates of intracellular transport mechanisms
    c.INaK_ = 1.362 ;                % pA/pF
    c.kNaCa = 1000 ;                 % pA/pF
    c.Iup_ = 0.425 ;                 % uM/ms
    c.Vleak = 8e-5 ;                 % ms-1
    
    % % values needed for calculation of NCX current
    % Several saturation constants, nomenclature gets tricky
    % NCX
    p.KmNa = 87500 ;                 % uM
    p.KmCa = 1380 ;                  % uM
    p.ksat = 0.1 ;                   % unitless
    p.alpha_ncx = 2.5 ;              % unitless
    p.eta = 0.35 ;                   % unitless, actually gamma in paper
    % Na-K pump
    p.KmNai = 40000 ;                % uM
    p.KmKo = 1000 ;                  % uM
    % Sarcolemmal Ca pump
    p.KpCa = 0.5 ;                    % uM
    % SERCA
    p.Kmup = 0.25 ;                  % uM
    
    c.arel = 16.464 ;                 % uM/ms
    p.brel = 250 ;                   % uM
    c.crel = 8.232 ;                  % uM/ms
    
    % % % Buffering parameters
    % % % Considered generic cytotolic and SR buffers
    p.Bufc = 150 ;                % uM
    p.Kbufc = 1 ;                   % uM
    p.BufSR = 10000 ;              % uM
    p.KbufSR = 300 ;                 % uM
    
    if strcmp(celltype,'epi')==1
        % epicardial cell
        c.Gto_ = 0.294 ;                % nS/pF
        c.GKs_ = 0.245 ;                % nS/pF
    elseif strcmp(celltype,'mid')==1
        % M cell
        c.Gto_ = 0.294 ;                % nS/pF
        c.GKs_ = 0.062 ;                 % nS/pF
    elseif strcmp(celltype,'endo')==1
        % endocardial cell
        c.Gto_ = 0.073 ;                % nS/pF
        c.GKs_ = 0.245 ;                % nS/pF
    else
        fprintf('Invalid cell type entered. Please re-enter cell type and try again.')
    end
    %% TT06
elseif strcmp(model_name,'TT06') 
    p.celltype = celltype;
    % Physical constants
    p.F = 96.4853415 ;              % Faraday's constant, coulombs/mmol
    p.R = 8.314472 ;                % gas constant, J/(K mol)
    p.T = 310.0;                    % absolute temperature, K
    p.RTF = p.R*p.T/p.F ;
    
    p.Cm = 2 ;                      % uF/cm2
    p.Acap = 5.6297*3.280e-5 ;             % cm2
    % % The factor of 5.6297 inserted as empirical correction factor
    % % ten Tusscher-Panfilov paper is vague about volumes, capacitances, etc.
    % % C code (from website) includes for instance, the variable
    % % CAPACITANCE=0.185;
    % % Units?  how is this used exactly?
    % % To account for this, I computed change in [Ca] with unitary current
    % % in both models.  Surface area needed to be increased by factor of 5.6297 to
    % % make the two derivatives equal
    
    p.Vmyo = 16.404 ;               % pL
    p.VSR = 1.094 ;                 % pL
    p.Vss = 0.05468 ;               % pL
    % Paper says these in units of (um)^3, but this is obviously wrong
    % 1 (um)^3 = 1 fL; myoplasmic volume more like 16 pL, or 16000 fL
    
    p.Ko = 5400 ;                   % uM
    p.Nao = 140000 ;                % uM
    p.Cao = 2000 ;                  % uM
    
    c.GNa_ = 14.838 ;               % nS/pF
    c.PCa_ = 3.980e-5 ;             % cm uF-1 ms-1
    c.GK1_ = 5.405 ;                % nS/pF
    c.GKr_ = 0.153 ;                % nS/pF
    
    p.GpK_ = 0.0146 ;               % nS/pF
    c.GpCa_ = 0.1238 ;              % nS/pF
    c.GNab_ = 2.9e-4 ;              % nS/pF
    c.GCab_ = 5.92e-4 ;             % nS/pF
    
    if strcmp(celltype,'endo')==1
        c.Gto_ = 0.073 ;                % nS/pF
        c.GKs_ = 0.392 ;                % nS/pF
    elseif strcmp(celltype,'mid')==1
        c.Gto_ = 0.294 ;                % nS/pF
        c.GKs_ = 0.098 ;                % nS/pF
    elseif strcmp(celltype,'epi')==1
        c.Gto_ = 0.294 ;                % nS/pF
        c.GKs_ = 0.392 ;                % nS/pF
    else
        fprintf('Invalid cell type entered. Please re-enter cell type and try again.')
    end
    
    p.pKNa = 0.03 ;                  % relative permeability, Na to K
    
    % Maximum rates of intracellular transport mechanisms
    c.INaK_ = 2.724 ;                % pA/pF
    c.kNaCa = 1000 ;                 % pA/pF
    c.Iup_ = 6.375 ;                 % uM/ms
    
    % % values needed for calculation of NCX current
    % Several saturation constants, nomenclature gets tricky
    % NCX
    p.KmNa = 87500 ;                 % uM
    p.KmCa = 1380 ;                  % uM
    p.ksat = 0.1 ;                   % unitless
    p.alpha_ncx = 2.5 ;              % unitless
    p.eta = 0.35 ;                   % unitless, actually gamma in paper
    % Na-K pump
    p.KmNai = 40000 ;                % uM
    p.KmKo = 1000 ;                  % uM
    % Sarcolemmal Ca pump
    p.KpCa = 0.5 ;                    % uM
    % SERCA
    p.Kmup = 0.25 ;                  % uM
    
    % Vrel = 40.8 ;                 % ms^-1
    c.Vrel = 0.102 ;                % ms^-1 as per KWTT source code
    p.Vxfer = 0.0038 ;              % ms^-1
    c.Vleak = 3.6e-4 ;              % ms^-1
    % Note, units given in paper must be incorrect
    % Paper says that these are in mM/ms
    % But equations multiply rate by concentration gradient
    % Thus no need to adjust rates for different concentration units
    
    p.k1_ryr_prime = 0.15e-6 ;        % uM-2 ms-1
    p.k2_ryr_prime = 0.045e-3 ;       % uM-1 ms-1
    p.k3_ryr = 0.06 ;                 % ms-1
    % k4_ryr = 1.5e-5 ;               % ms-1
    p.k4_ryr = 0.005 ;                % ms-1 as per KWTT source code
    p.maxsr = 2.5 ;                   % dimensionless
    p.minsr = 1 ;                     % dimensionless
    p.EC_ryr = 1500 ;                 % uM
    
    % % % Buffering parameters
    % % % Considered generic cytotolic and SR buffers
    p.Bufc = 200 ;                  % uM
    p.Kbufc = 1 ;                   % uM
    p.Bufss = 400 ;                 % uM
    p.Kbufss = 0.25 ;               % uM
    p.BufSR = 10000 ;               % uM
    p.KbufSR = 300 ;                % uM
    
    %% GRANDI
elseif strcmp(model_name,'Grandi')
    p.celltype = celltype;
    % Constants
    p.R = 8314;       % [J/kmol*K]
    p.Frdy = 96485;   % [C/mol]
    p.Temp = 310;     % [K]
    p.FoRT = p.Frdy/p.R/p.Temp;
    p.Cmem = 1.3810e-10;   % [F] membrane capacitance
    p.Qpow = (p.Temp-310)/10;
    
    % Cell geometry
    p.cellLength = 100;     % cell length [um]
    p.cellRadius = 10.25;   % cell radius [um]
    p.junctionLength = 160e-3;  % junc length [um]
    p.junctionRadius = 15e-3;   % junc radius [um]
    p.distSLcyto = 0.45;    % dist. SL to cytosol [um]
    p.distJuncSL = 0.5;  % dist. junc to SL [um]
    p.DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
    p.DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
    p.DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
    p.DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec]
    p.Vcell = pi*p.cellRadius^2*p.cellLength*1e-15;    % [L]
    p.Vmyo = 0.65*p.Vcell; p.Vsr = 0.035*p.Vcell; p.Vsl = 0.02*p.Vcell; p.Vjunc = 0.0539*.01*p.Vcell;
    p.SAjunc = 20150*pi*2*p.junctionLength*p.junctionRadius;  % [um^2]
    p.SAsl = pi*2*p.cellRadius*p.cellLength;          % [um^2]
    %J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
    %J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
    %J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
    %J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
    %J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 1.1074e-13
    %J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.5714e-12
    %J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 7.36e-13
    %J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.3056e-11
    % tau's from c-code, not used here
    p.J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
    p.J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
    p.J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
    p.J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11
    
    % Fractional currents in compartments
    p.Fjunc = 0.11;   p.Fsl = 1-p.Fjunc;
    p.Fjunc_CaL = 0.9; p.Fsl_CaL = 1-p.Fjunc_CaL;
    
    % Fixed ion concentrations
    p.Cli = 15;   % Intracellular Cl  [mM]
    p.Clo = 150;  % Extracellular Cl  [mM]
    p.Ko = 5.4;   % Extracellular K   [mM]
    p.Nao = 140;  % Extracellular Na  [mM]
    p.Cao = 1.8;  % Extracellular Ca  [mM]1.8
    p.Mgi = 1;    % Intracellular Mg  [mM]
    
    % Na transport parameters
    p.KmNaip = 11;         % [mM]11
    p.KmKo =1.5;         % [mM]1.5
    p.Q10NaK = 1.63;
    p.Q10KmNai = 1.39;
    
    % Cl current parameters
    p.KdClCa = 100e-3;    % [mM]
    
    % I_Ca parameters
    %%%% EAS: eventually modify this so that pNa/PCa_=constant
    %%%%%     for now do not worry about it
    p.pNa = 0.50*1.5e-8;       % [cm/sec]
    p.pK = 0.50*2.7e-7;        % [cm/sec]
    p.Q10CaL = 1.8;
    
    % %% Ca transport parameters
    p.KmCai = 3.59e-3;    % [mM]
    p.KmCao = 1.3;        % [mM]
    p.KmNai = 12.29;      % [mM]
    p.KmNao = 87.5;       % [mM]
    p.ksat = 0.32;        % [none]
    p.nu = 0.27;          % [none]
    p.Kdact = 0.150e-3;   % [mM]
    p.Q10NCX = 1.57;      % [none]
    p.KmPCa = 0.5e-3;     % [mM]
    p.Q10SLCaP = 2.35;    % [none]
    
    % SR flux parameters
    p.Q10SRCaP = 2.6;          % [none]
    p.Kmf = 0.246e-3;          % [mM] default
    %Kmf = 0.175e-3;          % [mM]
    p.Kmr = 1.7;               % [mM]L cytosol
    p.hillSRCaP = 1.787;       % [mM]
    p.koCa = 10;               % [mM^-2 1/ms]   %default 10   modified 20
    p.kom = 0.06;              % [1/ms]
    p.kiCa = 0.5;              % [1/mM/ms]
    p.kim = 0.005;             % [1/ms]
    p.ec50SR = 0.45;           % [mM]
    
    % Buffering parameters
    % koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
    p.Bmax_Naj = 7.561;       % [mM] % Bmax_Naj = 3.7; (c-code difference?)  % Na buffering
    p.Bmax_Nasl = 1.65;       % [mM]
    p.koff_na = 1e-3;         % [1/ms]
    p.kon_na = 0.1e-3;        % [1/mM/ms]
    p.Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
    p.koff_tncl = 19.6e-3;    % [1/ms]
    p.kon_tncl = 32.7;        % [1/mM/ms]
    p.Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity
    p.koff_tnchca = 0.032e-3; % [1/ms]
    p.kon_tnchca = 2.37;      % [1/mM/ms]
    p.koff_tnchmg = 3.33e-3;  % [1/ms]
    p.kon_tnchmg = 3e-3;      % [1/mM/ms]
    p.Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
    p.koff_cam = 238e-3;      % [1/ms]
    p.kon_cam = 34;           % [1/mM/ms]
    p.Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
    p.koff_myoca = 0.46e-3;   % [1/ms]
    p.kon_myoca = 13.8;       % [1/mM/ms]
    p.koff_myomg = 0.057e-3;  % [1/ms]
    p.kon_myomg = 0.0157;     % [1/mM/ms]
    p.Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
    p.koff_sr = 60e-3;        % [1/ms]
    p.kon_sr = 100;           % [1/mM/ms]
    p.Bmax_SLlowsl = 37.4e-3*p.Vmyo/p.Vsl;        % [mM]    % SL buffering
    p.Bmax_SLlowj = 4.6e-3*p.Vmyo/p.Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
    p.koff_sll = 1300e-3;     % [1/ms]
    p.kon_sll = 100;          % [1/mM/ms]
    p.Bmax_SLhighsl = 13.4e-3*p.Vmyo/p.Vsl;       % [mM]
    p.Bmax_SLhighj = 1.65e-3*p.Vmyo/p.Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
    p.koff_slh = 30e-3;       % [1/ms]
    p.kon_slh = 100;          % [1/mM/ms]
    p.Bmax_Csqn = 140e-3*p.Vmyo/p.Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
    p.koff_csqn = 65;         % [1/ms]
    p.kon_csqn = 100;         % [1/mM/ms]
    
    c.GNa=23;
    c.GNaB = 0.597e-3;    % [mS/uF] 0.897e-3
    c.PCa_ = 0.50*5.4e-4;       % [cm/sec]
    c.GCaB = 5.513e-4;    % [uA/uF] 3
    
    c.GKr_ = 0.035 ;
    c.GKs_ = 0.0035 ;
    c.GK1_ = 0.35 ;
    c.gkp = 2*0.001;
    c.GClCa =0.5* 0.109625;   % [mS/uF]
    c.GClB = 1*9e-3;        % [mS/uF]
    c.IbarNaK = 1.0*1.8;%1.90719;     % [uA/uF]
    c.IbarNCX = 1.0*4.5;      % [uA/uF]5.5 before - 9 in rabbit
    c.ks = 25;                 % [1/ms]
    c.Kleak_ = 5.348e-6 ;
    c.Vmax_SRCaP = 1.0*5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
    c.IbarSLCaP = 0.0673; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
    
    if strcmp(celltype,'epi')==1
        c.GtoSlow=1.0*0.13*0.12; %epi
        c.GtoFast=1.0*0.13*0.88; %epi0.88
    elseif strcmp(celltype,'endo')==1
        c.GtoSlow=0.13*0.3*0.964; %endo
        c.GtoFast=0.13*0.3*0.036; %endo
    else
        fprintf('Invalid cell type entered. Please re-enter cell type and try again.')
    end
    
    %% O'HARA
elseif strcmp(model_name,'Ohara') 
    
    p.celltype = celltype;
    p.Nao=140.0;
    p.Cao=1.8;
    p.Ko=5.4;
    
    %physical constants
    p.R=8314.0;
    p.T=310.0;
    p.F=96485.0;
    p.Cm=1.0;                     %uF
    
    %cell geometry
    p.L=0.01;
    p.rad=0.0011;
    p.vcell=1000*3.14*p.rad*p.rad*p.L;
    p.Ageo=2*3.14*p.rad*p.rad+2*3.14*p.rad*p.L;
    p.Acap=2*p.Ageo;
    p.vmyo=0.68*p.vcell;
    p.vnsr=0.0552*p.vcell;
    p.vjsr=0.0048*p.vcell;
    p.vss=0.02*p.vcell;
    
    %jsr constants
    p.bt=4.75;
    p.a_rel=0.5*p.bt;
    
    % computed quantities that do not change during simulation
    c.GNa=75;
    c.GNaL=0.0075;
    c.Gto=0.02;
    c.GKr_=0.046;
    c.GKs_=0.0034;
    c.GK1=0.1908;
    c.Gncx=0.0008;
    c.GKb=0.003;
    c.GpCa=0.0005;
    c.PCa_=0.0001;
    c.Pnak=30;
    
    if  strcmp(celltype,'epi')==1
        c.GNaL=c.GNaL*0.6;
        c.Gto=c.Gto*4.0;
        c.GKr_=c.GKr_*1.3;
        c.GKs_=c.GKs_*1.4;
        c.GK1=c.GK1*1.2;
        c.Gncx=c.Gncx*1.1;
        c.GKb=c.GKb*0.6;
        c.PCa_=c.PCa_*1.2;
        c.Pnak=c.Pnak*0.9;
        
    elseif  strcmp(celltype,'mid')==1
        c.Gto=c.Gto*4.0;
        c.GKr_=c.GKr_*0.8;
        c.GK1=c.GK1*1.3;
        c.Gncx=c.Gncx*1.4;
        c.PCa_=c.PCa_*2.5;
        c.Pnak=c.Pnak*0.7;
        
    end
    c.PNab=3.75e-10;
    c.PCab=2.5e-8;
    
    c.SERCA_total = 1 ;
    c.RyR_total = 1 ;
      
else
    fprintf('Please enter a valid model name.')
end
