function [ICaL,IKr,IKs]=currents(model_name,state_vars,p,c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- "Slow delayed rectifier current protects ventricular myocytes from
% arrhythmic dynamics across multiple species: a computational study" ---%

% By: Varshneya,Devenyi,Sobie 
% For questions, please contact Dr.Eric A Sobie -> eric.sobie@mssm.edu 
% or put in a pull request or open an issue on the github repository:
% https://github.com/meeravarshneya1234/IKs_stabilizes_APs.git. 

%--- Note:
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
                            %% -- currents.m -- %%
% Description: Calculates the ICaL, IKs, and IKr currents for a single AP. 

% Inputs:
% --> model_name - [string] name of model for which the parameters are needed.
% --> state_vars - [cell array] state variables calculated in dydt, used to
% determine resulting currents
% --> p - [struct array] model parameters 
% --> c - [struct array] model parameters that were varied to create a
% population

% Outputs:
% --> ICaL - [array] L-type calcium current 
% --> IKs  - [array] slow delayed rectifier potassium current  
% --> IKr -  [array] rapid delayed rectifier potassium current 
%--------------------------------------------------------------------------

%% FOX
if strcmp(model_name,'Fox')
    [V,Cai,CaSR,f,d,m,h,j,fCa,xKr,xKs,xto,yto] = deal(state_vars{:});
    
    % % % ICaL
    ICa_ = c.PCa_/p.Cm*4*p.F/p.RTF*V.*(Cai.*exp(2*V/p.RTF) - 0.341*p.Cao)./ ...
        (exp(2*V/p.RTF) - 1) ;
    ICaL = ICa_.*d.*f.*fCa ;
    
    % % % IKr
    EK = p.RTF*log(p.Ko/p.Ki) ;
    R_V = 1./(2.5*exp(0.1*(V+28)) + 1) ;
    IKr = c.GKr_*R_V.*xKr*sqrt(p.Ko/4000).*(V - EK) ;
    
    % % % IKs
    EKs = p.RTF*log((p.Ko+0.01833*p.Nao)/(p.Ki+0.01833*p.Nai)) ;
    IKs = c.GKs_*xKs.^2.*(V - EKs) ;
    
    %% HUND
elseif strcmp(model_name,'Hund')
    [V,Cai,Cass,CaNSR,CaJSR,Nai,Ki,Cli,m,h,j,hL, ...
        d,dp,f,f2,fCa,fCa2,xKr,xs1,xs2,ato,ito1,ito2,aa,ro,ri,CAMK_trap] = deal(state_vars{:});
    
    % % % ICaL
    ICa_ = c.PCa_*4*p.F/p.RTF*(V-15).* ...
        (p.gamma_Cai*Cass.*exp(2*(V-15)/p.RTF) - p.gamma_Cao*p.Cao)./ ...
        (exp(2*(V-15)/p.RTF) - 1) ;
    ICaL = ICa_.*d.^(dp).*f.*f2.*fCa.*fCa2 ;
    
    % % % IKr
    EK = p.RTF*log(p.Ko./Ki) ;
    RKr = 1./(exp((V+10)/15.4) + 1) ;
    IKr = c.GKr_*sqrt(p.Ko/5.4)*xKr.*RKr.*(V - EK) ;
    
    % % % IKs
    EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)./(Ki + p.pKNa*Nai)) ;
    IKs = c.GKs_*(1+0.6./((3.8e-5./Cai).^(1.4)+1)).*xs1.*xs2.*(V - EKs) ;
    
    %% LIVSHITZ
elseif strcmp(model_name,'Livshitz')
    [V,Cai,CaNSR,CaJSR,Nai,Ki,m,h,j,d,f,b,g,xKr,xs1,xs2,Jrel] = deal(state_vars{:});
    
    % % % ICaL
    ICa_ = c.PCa_*4*p.F*p.FRT.*V.* ...
        (p.gamma_Cai.*Cai.*exp(2*V*p.FRT) - p.gamma_Cao*p.Cao)./ ...
        (exp(2*V*p.FRT) - 1) ;
    ICaK_ = c.PCa_K*p.F*p.FRT.*V.* ...
        (p.gamma_Ki.*Ki.*exp(V*p.FRT) - p.gamma_Ko*p.Ko)./ ...
        (exp(V*p.FRT) - 1) ;
    ICaNa_ = c.PCa_Na*p.F*p.FRT.*V.* ...
        (p.gamma_Nai.*Nai.*exp(V*p.FRT) - p.gamma_Nao*p.Nao)./ ...
        (exp(V*p.FRT) - 1) ;
    fCa = 1./(Cai/p.KmCa + 1) ;
    ICa = ICa_.*d.*f.*fCa ;
    ICaL_K = ICaK_.*d.*f.*fCa ;
    ICaL_Na = ICaNa_.*d.*f.*fCa ;
    ICaL = ICa+ICaL_K+ICaL_Na;
    
    % % % IKr
    EK = p.RTF*log(p.Ko./Ki) ;
    RKr = 1./(exp((V+9)/22.4) + 1) ;
    IKr = c.GKr_*sqrt(p.Ko./5.4).*xKr.*RKr.*(V - EK) ;
    
    % % % IKs
    EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)./(Ki + p.pKNa*Nai)) ;
    IKs = c.GKs_*(1+0.6./((3.8e-5./Cai).^(1.4)+1)).*xs1.*xs2.*(V - EKs) ;
    
    
    %% DEVENYI
elseif strcmp(model_name,'Devenyi')
    
    [V,Cai,CaNSR,CaJSR,Nai,Ki,m,h,j,d,f,b,g,xKr,xs1,xs2,Jrel]  = deal(state_vars{:});
    
    % % % ICaL
    ICa_ = p.KCaL*c.PCa_*4*p.F*p.FRT.*V.* ...
        (p.gamma_Cai.*Cai.*exp(2.*V*p.FRT) - p.gamma_Cao*p.Cao)./ ...
        (exp(2.*V*p.FRT) - 1) ;
    ICaK_ = p.KCaL*c.PCa_K*p.F*p.FRT.*V.* ...
        (p.gamma_Ki.*Ki.*exp(V*p.FRT) - p.gamma_Ko*p.Ko)./ ...
        (exp(V.*p.FRT) - 1) ;
    ICaNa_ = p.KCaL*c.PCa_Na*p.F*p.FRT.*V.* ...
        (p.gamma_Nai.*Nai.*exp(V*p.FRT) - p.gamma_Nao*p.Nao)./ ...
        (exp(V.*p.FRT) - 1) ;
    fCa = 1./(Cai./p.KmCa + 1) ;
    ICa = ICa_.*d.*f.*fCa ;
    ICaL_K = ICaK_.*d.*f.*fCa ;
    ICaL_Na = ICaNa_.*d.*f.*fCa ;
    ICaL = ICa+ICaL_K+ICaL_Na;
    
    % % % IKr
    EK = p.RTF.*log(p.Ko./Ki) ;
    RKr = 1./(exp((V+9)/22.4) + 1) ;
    IKr =c.GKr_*sqrt(p.Ko/5.4).*xKr.*RKr.*(V - EK) ;
    
    % % % IKs
    EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)./(Ki + p.pKNa*Nai)) ;
    IKs = c.GKs_*(1+0.6./((3.8e-5./Cai).^(1.4)+1)).*xs1.*xs2.*(V - EKs) ;
    
    %% SHANNON
elseif strcmp(model_name,'Shannon')
    
    [Ca_Calsequestrin,Ca_SL,Ca_SLB_SL,Ca_SLB_jct,Ca_SLHigh_SL, ...
        Ca_SLHigh_jct,Ca_SR,Ca_jct,Cai,d,fCaB_SL,fCaB_jct,f,Xr,Xs,h,j,m,X_tof, ...
        Y_tof,R_tos,X_tos,Y_tos,I,O,R1,Na_SL,Na_SL_buf,Na_jct,Na_jct_buf,Nai,V, ...
        Ca_Calmodulin,Ca_Myosin,Ca_SRB,Ca_TroponinC,Ca_TroponinC_Ca_Mg,Mg_Myosin, ...
        Mg_TroponinC_Ca_Mg]  = deal(state_vars{:});
    
    % % % ICaL
    Q_CaL = p.Q10_CaL^((p.T-310.0)/10.0);
    temp = 0.45.*d.*f.*Q_CaL.*V*p.F^2.0/(p.R2*p.T);
    fCa_jct = 1.0-fCaB_jct;
    i_CaL_Ca_jct = temp.*fCa_jct*p.Fx_ICaL_jct*c.PCa_*4.0.*(p.gamma_Cai*Ca_jct.*exp(2.0*V*p.F/(p.R2*p.T))-p.gamma_Cao*p.Cao)./(exp(2.0*V*p.F/(p.R2*p.T))-1.0);
    fCa_SL = 1.0-fCaB_SL;
    i_CaL_Ca_SL = temp.*fCa_SL*p.Fx_ICaL_SL*c.PCa_*4.0.*(p.gamma_Cai.*Ca_SL.*exp(2.0*V*p.F/(p.R2*p.T))-p.gamma_Cao*p.Cao)./(exp(2.0*V*p.F/(p.R2*p.T))-1.0);
    i_CaL_Na_SL = temp.*fCa_SL*p.Fx_ICaL_SL*p.PNa.*(p.gamma_Nai.*Na_SL.*exp(V*p.F/(p.R2*p.T))-p.gamma_Nao*p.Nao)./(exp(V*p.F/(p.R2*p.T))-1.0);
    i_CaL_Na_jct = temp.*fCa_jct*p.Fx_ICaL_jct*p.PNa.*(p.gamma_Nai.*Na_jct.*exp(V*p.F/(p.R2*p.T))-p.gamma_Nao*p.Nao)./(exp(V*p.F/(p.R2*p.T))-1.0);
    i_CaL_K = temp.*(fCa_SL*p.Fx_ICaL_SL+fCa_jct*p.Fx_ICaL_jct).*p.PK.*(p.gamma_Ki*p.Ki.*exp(V*p.F/(p.R2*p.T))-p.gamma_Ko*p.Ko)./(exp(V*p.F/(p.R2*p.T))-1.0);
    ICaL = i_CaL_Ca_SL+i_CaL_Ca_jct+i_CaL_Na_SL+i_CaL_Na_jct+i_CaL_K;
    
    % % % IKr
    E_K = p.R2*p.T/p.F*log(p.Ko/p.Ki);
    G_IKr = c.GKr_*sqrt(p.Ko/5.4);
    Rr = 1.0./(1.0+exp((33.0+V)/22.4));
    IKr = G_IKr*Xr.*Rr.*(V-E_K);
    
    % % % IKs
    pCa_jct = -log10(Ca_jct/1.0)+3.0;
    pCa_SL = -log10(Ca_SL/1.0)+3.0;
    G_Ks_jct = c.GKs_.*(0.057+0.19./(1.0+exp((-7.2+pCa_jct)./0.6)));
    G_Ks_SL = c.GKs_.*(0.057+0.19./(1.0+exp((-7.2+pCa_SL)./0.6)));
    E_Ks_jct = p.R2*p.T./p.F*log((p.Ko+p.pKNa*p.Nao)./(p.Ki+p.pKNa*Na_jct));
    E_Ks_SL = p.R2*p.T./p.F*log((p.Ko+p.pKNa*p.Nao)./(p.Ki+p.pKNa*Na_SL_buf));
    E_Ks = p.R2*p.T./p.F*log((p.Ko+p.pKNa*p.Nao)./(p.Ki+p.pKNa*Nai));
    i_Ks_jct = p.Fx_Ks_jct.*G_Ks_jct.*Xs.^2.*(V-E_Ks);
    i_Ks_SL = p.Fx_Ks_SL.*G_Ks_SL.*Xs.^2.0.*(V-E_Ks);
    IKs = i_Ks_jct+i_Ks_SL;
    
    %% TT04
elseif strcmp(model_name,'TT04')
    
    [V,m,h,j,d,f,fCa,r,s,xs,xr1,xr2,g,Cai,CaSR,Nai,Ki]= deal(state_vars{:});
    
    % % % ICaL
    ICaL= c.PCa_*d.*f.*fCa*4*p.F/p.RTF.*V.*(Cai.*exp(2*V/p.RTF) - 0.341*p.Cao)./ ...
        (exp(2*V/p.RTF) - 1) ;
    
    % % % IKr
    EK = p.RTF*log(p.Ko./Ki) ;
    IKr = c.GKr_*sqrt(p.Ko/5400)*xr1.*xr2.*(V-EK) ;
    
    % % % IKs
    EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)./(Ki + p.pKNa*Nai)) ;
    IKs = c.GKs_*xs.^2.*(V-EKs) ;
    
    
    %% TT06
elseif strcmp(model_name,'TT06') || strcmp(model_name,'TT06_opt')
    
    [V,m,h,j,d,f,f2,fCass,r,s,xs,xr1,xr2, Rbar_ryr,Cai,Cass,CaSR,Nai,Ki]...
        = deal(state_vars{:});
    
    % % % ICaL
    ICaL = c.PCa_*d.*f.*f2.*fCass*4*p.F/p.RTF.*(V-15).* ...
        (0.25*Cai.*exp(2*(V-15)/p.RTF) - p.Cao)./ ...
        (exp(2*(V-15)/p.RTF) - 1) ;
    
    % % % IKr
    EK = p.RTF*log(p.Ko./Ki) ;
    IKr = c.GKr_*sqrt(p.Ko/5400)*xr1.*xr2.*(V-EK) ;
    
    % % % IKs
    EKs = p.RTF*log((p.Ko + p.pKNa*p.Nao)./(Ki + p.pKNa*Nai)) ;
    IKs = c.GKs_*xs.^2.*(V-EKs) ;
    
    %% GRANDI
elseif strcmp(model_name,'Grandi') || strcmp(model_name,'Grandi_opt')
    
    [m, h, j, d, f, fcaBj, fcaBsl, xtos, ytos, ...
        xtof, ytof, xkr, xks,RyRr, RyRo, RyRi, NaBj, NaBsl, ...
        TnCL, TnCHc, TnCHm, CaM, Myoc, Myom,SRB, SLLj, SLLsl,...
        SLHj, SLHsl, Csqnb, Ca_sr, Naj, Nasl, Nai, Ki, Caj,...
        Casl, Cai, V] = deal(state_vars{:});
    
    % % % ICaL
    fcaCaMSL=0;
    fcaCaj= 0;
    ibarca_j = c.PCa_*4*(V* p.Frdy* p.FoRT).*(0.341*Caj.*exp(2*V*p.FoRT)-0.341*p.Cao)./(exp(2*V*p.FoRT)-1);
    ibarca_sl = c.PCa_*4*(V*p.Frdy*p.FoRT).*(0.341*Casl.*exp(2*V*p.FoRT)-0.341*p.Cao)./(exp(2*V*p.FoRT)-1);
    ibark = p.pK*(V*p.Frdy*p.FoRT).*(0.75*Ki.*exp(V*p.FoRT)-0.75*p.Ko)./(exp(V*p.FoRT)-1);
    ibarna_j = p.pNa*(V*p.Frdy*p.FoRT).*(0.75*Naj.*exp(V*p.FoRT)-0.75*p.Nao)./(exp(V*p.FoRT)-1);
    ibarna_sl = p.pNa*(V*p.Frdy*p.FoRT).*(0.75*Nasl.*exp(V*p.FoRT)-0.75*p.Nao)./(exp(V*p.FoRT)-1);
    I_Ca_junc = (p.Fjunc_CaL*ibarca_j.*d.*f.*((1-fcaBj)+fcaCaj).*p.Q10CaL^p.Qpow)*0.45*1;
    I_Ca_sl = (p.Fsl_CaL*ibarca_sl.*d.*f.*((1-fcaBsl)+fcaCaMSL).*p.Q10CaL^p.Qpow)*0.45*1;
    I_Ca = I_Ca_junc+I_Ca_sl;
    I_CaK = (ibark.*d.*f.*(p.Fjunc_CaL*(fcaCaj+(1-fcaBj))+p.Fsl_CaL*(fcaCaMSL+(1-fcaBsl)))*p.Q10CaL^p.Qpow)*0.45*1;
    I_CaNa_junc = (p.Fjunc_CaL*ibarna_j.*d.*f.*((1-fcaBj)+fcaCaj)*p.Q10CaL^p.Qpow)*0.45*1;
    I_CaNa_sl = (p.Fsl_CaL*ibarna_sl.*d.*f.*((1-fcaBsl)+fcaCaMSL)*p.Q10CaL^p.Qpow)*.45*1;
    I_CaNa = I_CaNa_junc+I_CaNa_sl;
    ICaL = I_Ca+I_CaK+I_CaNa;
    
    % % % IKr
    pNaK = 0.01833;
    ek = (1/p.FoRT)*log(p.Ko./Ki);	        % [mV]
    gkr = c.GKr_*sqrt(p.Ko/5.4);
    rkr = 1./(1+exp((V+74)/24));
    IKr = gkr.*xkr.*rkr.*(V-ek);
    
    % % % IKs
    eks = (1/p.FoRT)*log((p.Ko+pNaK*p.Nao)./(Ki+pNaK.*Nai));
    gks_junc = c.GKs_ ;   %[mS/uF]
    gks_sl = c.GKs_ ;     %[mS/uF]
    Fjunc = 0.11;
    Fsl = 1-Fjunc;
    I_ks_junc = Fjunc*gks_junc.*xks.^2.*(V-eks);
    I_ks_sl = Fsl.*gks_sl.*xks.^2.*(V-eks);
    IKs = I_ks_junc+I_ks_sl;
    
    %% O'HARA
elseif strcmp(model_name,'Ohara') || strcmp(model_name,'Ohara_opt')
    
    [V, Nai, Nass, Ki, Kss, Cai, Cass, Cansr, Cajsr, m, hf, hs, j, hsp, jp, mL, hL,...
        hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp,...
        xrf, xrs, xs1, xs2, xk1, Jrelnp, Jrelp, CaMKt] = deal(state_vars{:});
    
    % % % ICaL
    KmCaMK=0.15;
    CaMKo=0.05;
    KmCaM=0.0015;
    CaMKb=CaMKo*(1.0-CaMKt)./(1.0+KmCaM./Cass);
    CaMKa=CaMKb+CaMKt;
    vffrt=V*p.F*p.F/(p.R*p.T);
    vfrt=V*p.F/(p.R*p.T);
    Aff=0.6;
    Afs=1.0-Aff;
    f=Aff*ff+Afs*fs;
    Afcaf=0.3+0.6./(1.0+exp((V-10.0)/10.0));
    Afcas=1.0-Afcaf;
    fca=Afcaf.*fcaf+Afcas.*fcas;
    fp=Aff.*ffp+Afs.*fs;
    fcap=Afcaf.*fcafp+Afcas.*fcas;
    PhiCaL=4.0*vffrt.*(Cass.*exp(2.0*vfrt)-0.341*p.Cao)./(exp(2.0*vfrt)-1.0);
    PhiCaNa=1.0*vffrt.*(0.75*Nass.*exp(1.0*vfrt)-0.75*p.Nao)./(exp(1.0*vfrt)-1.0);
    PhiCaK=1.0*vffrt.*(0.75*Kss.*exp(1.0*vfrt)-0.75*p.Ko)./(exp(1.0*vfrt)-1.0);
    PCap=1.1*c.PCa_;
    PCaNa=0.00125*c.PCa_;
    PCaK=3.574e-4*c.PCa_;
    PCaNap=0.00125*PCap;
    PCaKp=3.574e-4*PCap;
    fICaLp=(1.0./(1.0+KmCaMK./CaMKa));
    ICa_L=(1.0-fICaLp).*c.PCa_.*PhiCaL.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCap.*PhiCaL.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
    ICaNa=(1.0-fICaLp).*PCaNa.*PhiCaNa.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCaNap.*PhiCaNa.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
    ICaK=(1.0-fICaLp).*PCaK.*PhiCaK.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCaKp.*PhiCaK.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
    ICaL = ICa_L + ICaNa + ICaK;
    
    % % % IKr
    EK=(p.R*p.T/p.F)*log(p.Ko./Ki);
    Axrf=1.0./(1.0+exp((V+54.81)/38.21));
    Axrs=1.0-Axrf;
    xr=Axrf.*xrf+Axrs.*xrs;
    rkr=1.0./(1.0+exp((V+55.0)/75.0))*1.0./(1.0+exp((V-10.0)/30.0));
    IKr=c.GKr_*sqrt(p.Ko/5.4).*xr.*rkr.*(V-EK);
    
    % % % IKs
    PKNa=0.01833;
    EKs=(p.R*p.T/p.F)*log((p.Ko+PKNa*p.Nao)./(Ki+PKNa.*Nai));
    KsCa=1.0+0.6./(1.0+(3.8e-5./Cai).^1.4);
    IKs=c.GKs_.*KsCa.*xs1.*xs2.*(V-EKs);
    
    %%
else
    fprintf('Please enter a valid model name.')
end