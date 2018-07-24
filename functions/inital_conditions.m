% -- parameters for each of the models used in the study
function ic = inital_conditions(model_name)

%% FOX
if strcmp(model_name,'Fox')
    ic.V = -94.7 ;
    ic.Cai = 0.0472 ;
    ic.CaSR = 320 ;
    ic.f = 0.983 ;
    ic.d = 1e-4 ;
    ic.m = 2.46e-4 ;
    ic.h = 0.99869 ;
    ic.j = 0.99887 ;
    ic.fCa = 0.942 ;
    ic.xKr = 0.229 ;
    ic.xKs = 1e-4 ;
    ic.xto = 3.742e-5 ;
    ic.yto = 1 ;
    
    %% HUND
elseif strcmp(model_name,'Hund')
    ic.V = -94.7 ;
    ic.Cai = 0.0822e-3 ;
    ic.Cass = 0.0822e-3 ;
    ic.CaJSR = 1.25 ;
    ic.CaNSR = 1.25 ;
    ic.Nai = 9.71 ;
    ic.Ki = 142.82 ;
    ic.Cli = 19.53 ;
    ic.m = 2.46e-4 ;
    ic.h = 0.99869 ;
    ic.j = 0.99887 ;
    ic.hL = 0.99869 ;
    ic.d = 1e-4 ;
    ic.dp = 1 ;
    ic.f = 0.983 ;
    ic.f2 = 0.983 ;
    ic.fCa = 0.942 ;
    ic.fCa2 = 0.942 ;
    ic.xKr = 0.229 ;
    ic.xs1 = 1e-4 ;
    ic.xs2 = 1e-4 ;
    ic.ato = 3.742e-5 ;
    ic.ito1 = 1 ;
    ic.ito2 = 1 ;
    ic.aa = 0 ;
    ic.ro = 0 ;
    ic.ri = 0 ;
    ic.CAMK_trap = 0.001 ;
    
    %% LIVSHITZ
elseif strcmp(model_name,'Livshitz')
    ic.V = -84.7 ;
    ic.Cai = 0.0822e-3 ;
    ic.CaNSR = 1.25 ;
    ic.CaJSR = 1.25 ;
    ic.Nai = 9.71 ;
    ic.Ki = 142.82 ;
    ic.m = 2.46e-4 ;
    ic.h = 0.99869 ;
    ic.j = 0.99887 ;
    ic.d = 1e-4 ;
    ic.f = 0.983 ;
    ic.b = 1e-4 ;
    ic.g = 0.983 ;
    ic.xKr = 0.229 ;
    ic.xs1 = 1e-4 ;
    ic.xs2 = 1e-4 ;
    ic.Jrel = 1e-4 ;
    
    %% DEVENYI
elseif strcmp(model_name,'Devenyi')
    ic.V = -84.7 ;
    ic.Cai = 0.0822e-3 ;
    ic.CaNSR = 1.25 ;
    ic.CaJSR = 1.25 ;
    ic.Nai = 9.71 ;
    ic.Ki = 142.82 ;
    ic.m = 2.46e-4 ;
    ic.h = 0.99869 ;
    ic.j = 0.99887 ;
    ic.d = 1e-4 ;
    ic.f = 0.983 ;
    ic.b = 1e-4 ;
    ic.g = 0.983 ;
    ic.xKr = 0.229 ;
    ic.xs1 = 1e-4 ;
    ic.xs2 = 1e-4 ;
    ic.Jrel = 1e-4 ;
    
    %% SHANNON
elseif strcmp(model_name,'Shannon')
    ic.Ca_Calsequestrin = 1.1865 ;
    ic.Ca_SL = 0.0001064 ;
    ic.Ca_SLB_SL = 0.0098686 ;
    ic.Ca_SLB_jct = 0.0077808 ;
    ic.Ca_SLHigh_SL = 0.11444 ;
    ic.Ca_SLHigh_jct = 0.077504 ;
    ic.Ca_SR = 0.54561 ;
    ic.Ca_jct = 0.00017484 ;
    ic.Cai = 8.735e-005 ;
    ic.d = 6.9975e-006 ;
    ic.fCaB_SL = 0.015353 ;
    ic.fCaB_jct = 0.024609 ;
    ic.f = 1.0007 ;
    ic.Xr = 0.0084716 ;
    ic.Xs = 0.006874 ;
    ic.h = 0.98714 ;
    ic.j = 0.99182 ;
    ic.m = 0.0013707 ;
    ic.X_tof = 0.0040112 ;
    ic.Y_tof = 0.99463 ;
    ic.R_tos = 0.38343 ;
    ic.X_tos = 0.0040113 ;
    ic.Y_tos = 0.29352 ;
    ic.I = 9.272e-008 ;
    ic.O = 7.1126e-007 ;
    ic.R1 = 0.88467 ;
    ic.Na_SL = 8.8741 ;
    ic.Na_SL_buf = 0.77612 ;
    ic.Na_jct = 8.8728 ;
    ic.Na_jct_buf = 3.5571 ;
    ic.Nai = 8.8745 ;
    ic.V = -85.7197 ;
    ic.Ca_Calmodulin = 0.00029596 ;
    ic.Ca_Myosin = 0.0019847 ;
    ic.Ca_SRB = 0.0021771 ;
    ic.Ca_TroponinC = 0.0089637 ;
    ic.Ca_TroponinC_Ca_Mg = 0.118 ;
    ic.Mg_Myosin = 0.1375 ;
    ic.Mg_TroponinC_Ca_Mg = 0.010338 ;
    
    %% TT04
elseif strcmp(model_name,'TT04')
    ic.V = -85 ;
    ic.m = 0 ;
    ic.h = 1 ;
    ic.j = 1 ;
    ic.d = 0 ;
    ic.f = 1 ;
    ic.fCa = 1 ;
    ic.r = 0 ;
    ic.s = 1 ;
    ic.xs = 0 ;
    ic.xr1 = 1 ;
    ic.xr2 = 1 ;
    ic.g = 0 ;
    ic.Cai = 0.1 ;
    ic.CaSR = 1000 ;
    ic.Nai = 10000 ;
    ic.Ki = 130000 ;
    
    %% TT06
elseif strcmp(model_name,'TT06') || strcmp(model_name,'TT06_opt')
    ic.V = -83.5092 ;
    ic.m = 0.0025 ;
    ic.h = 0.6945 ;
    ic.j = 0.6924 ;
    ic.d = 4.2418e-005 ;
    ic.f = 0.9697 ;
    ic.f2 = 0.9784 ;
    ic.fCass = 0.9999 ;
    ic.r = 3.2195e-008 ;
    ic.s = 1.0000 ;
    ic.xs = 0.0038 ;
    ic.xr1 = 3.1298e-004 ;
    ic.xr2 = 0.4534 ;
    ic.Rbar_ryr = 0.9816 ;
    ic.Cai = 0.1061 ;
    ic.Cass = 0.2381 ;
    ic.CaSR = 3.6426e+003 ;
    ic.Nai = 3.8067e+003 ;
    ic.Ki = 1.2369e+005 ;
    
    %% GRANDI
elseif strcmp(model_name,'Grandi') || strcmp(model_name,'Grandi_opt')
    ic.mo=1.405627e-3;
    ic.ho= 9.867005e-1;
    ic.jo=9.915620e-1;
    ic.do=7.175662e-6;
    ic.fo=1.000681;
    ic.fcaBjo=2.421991e-2;
    ic.fcaBslo=1.452605e-2;
    ic.xtoso=4.051574e-3;
    ic.ytoso=9.945511e-1;
    ic.xtofo=4.051574e-3;
    ic.ytofo= 9.945511e-1;
    ic.xkro=8.641386e-3;
    ic.xkso= 5.412034e-3;
    ic.RyRro=8.884332e-1;
    ic.RyRoo=8.156628e-7;
    ic.RyRio=1.024274e-7;
    ic.NaBjo=3.539892;
    ic.NaBslo=7.720854e-1;
    ic.TnCLo=8.773191e-3;
    ic.TnCHco=1.078283e-1;
    ic.TnCHmo=1.524002e-2;
    ic.CaMo=2.911916e-4;
    ic.Myoco=1.298754e-3;
    ic.Myomo=1.381982e-1;
    ic.SRBo=2.143165e-3;
    ic.SLLjo=9.566355e-3;
    ic.SLLslo=1.110363e-1;
    ic.SLHjo=7.347888e-3;
    ic.SLHslo=7.297378e-2;
    ic.Csqnbo= 1.242988;
    ic.Ca_sro=0.1e-1; %5.545201e-1;
    ic.Najo=9.06;%8.80329;
    ic.Naslo=9.06;%8.80733;
    ic.Naio=9.06;%8.80853;
    ic.Kio=120;
    ic.Cajo=1.737475e-4;
    ic.Caslo= 1.031812e-4;
    ic.Caio=8.597401e-5;
    ic.V=-8.09763e+1;
    
    %% O'HARA
elseif strcmp(model_name,'Ohara') || strcmp(model_name,'Ohara_opt')
    ic.V=-87;
    ic.Nai=7;
    ic.Nass=ic.Nai;
    ic.Ki=145;
    ic.Kss=ic.Ki;
    ic.Cai=1.0e-4;
    ic.Cass=ic.Cai;
    ic.Cansr=1.2;
    ic.Cajsr=ic.Cansr;
    ic.m=0;
    ic.hf=1;
    ic.hs=1;
    ic.j=1;
    ic.hsp=1;
    ic.jp=1;
    ic.mL=0;
    ic.hL=1;
    ic.hLp=1;
    ic.a=0;
    ic.iF=1;
    ic.iS=1;
    ic.ap=0;
    ic.iFp=1;
    ic.iSp=1;
    ic.d=0;
    ic.ff=1;
    ic.fs=1;
    ic.fcaf=1;
    ic.fcas=1;
    ic.jca=1;
    ic.nca=0;
    ic.ffp=1;
    ic.fcafp=1;
    ic.xrf=0;
    ic.xrs=0;
    ic.xs1=0;
    ic.xs2=0;
    ic.xk1=1;
    ic.Jrelnp=0;
    ic.Jrelp=0;
    ic.CaMKt=0;
    
    %%
else
    fprintf('Please enter a valid model name.')
end
