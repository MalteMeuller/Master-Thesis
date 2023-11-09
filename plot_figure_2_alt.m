% this code plots figure 2 based on the IRF used in irf_VAR_MCMC_STRUCT_ALT

clear all
close all

flag_plot=0; %change here depending on which plot


%%
load('code_fullsample_linear.mat')

[irf_DimpLIN,irf_CImeanLIN,CIupLIN,CIlowLIN,IRFseLIN]=irf_VAR_MCMC_struct_alt(beta0E,Omega0E,Omega0mat,beta0matF,s_desc);

Omega0Ex1=zeros(size(Omega0E));
for ii=1:length(Amat(s_MCMC.dropfirstT:end-1,:))
    [Omega0Ex]=vec2matScov_linear(Amat(s_MCMC.dropfirstT+ii-1,:)',s_desc.Omega_length_cov);
    Omega0Ex1=Omega0Ex1+Omega0Ex;
end

Omega0E=Omega0Ex1/length(Amat(s_MCMC.dropfirstT:end-1,:));




load('code_fullsample.mat')

[irf_DimpREC,irf_CImeanREC,CIupREC,CIlowREC,IRFseREC]=irf_VAR_MCMC_struct_alt(beta1E,Omega1E,Omega1mat,beta1matF,s_desc);
[irf_DimpEXP,irf_CImeanEXP,CIupEXP,CIlowEXP,IRFseEXP]=irf_VAR_MCMC_struct_alt(beta0E,Omega0E,Omega0mat,beta0matF,s_desc);

Omega0Ex1=zeros(size(Omega0E));
Omega1Ex1=zeros(size(Omega0E));
for ii=1:length(Amat(s_MCMC.dropfirstT:end-1,:))
    [Omega0Ex,Omega1Ex,thetaEx]=vec2matScov(Amat(s_MCMC.dropfirstT+ii-1,:)',s_desc.Omega_length_cov);
    Omega0Ex1=Omega0Ex1+Omega0Ex;
    Omega1Ex1=Omega1Ex1+Omega1Ex;
end

Omega1E=Omega1Ex1/length(Amat(s_MCMC.dropfirstT:end-1,:));
Omega0E=Omega0Ex1/length(Amat(s_MCMC.dropfirstT:end-1,:));






%%

%=====================================================================
%             Impulse responses in the linear model (CI lin)
%=====================================================================
if flag_plot==0

    figure('Position',[100 100 600 325])
    figure(1)
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupLIN(s_desc.KPIF,:),fliplr(CIlowLIN(s_desc.KPIF,:))]; %which shocks to use CI on
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 20], [0 0], 'k', 'LineWidth', 1)
    hold off
    title('KIX shock => KPIF response')
    xlim([1 20])
    legend('95% CI LIN','Linear','High','Low')

    %CI high
    figure('Position',[100 100 600 325])
    figure(2)
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupEXP(s_desc.KPIF,:),fliplr(CIlowEXP(s_desc.KPIF,:))];
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 20], [0 0], 'k', 'LineWidth', 1)
    hold off
    title('KIX shock => KPIF response')
    xlim([1 20])
    legend('95% CI HIGH','Linear','High','Low')

    % low inflation  
    figure('Position',[100 100 600 325])
    figure(3)
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupREC(s_desc.KPIF,:),fliplr(CIlowREC(s_desc.KPIF,:))];
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 20], [0 0], 'k', 'LineWidth', 1)
    hold off
    title('KIX shock => KPIF response')
    xlim([1 20])
    legend('95% CI LOW','Linear','High','Low')

  


end

%% calculate the exchange rate pass through 
%Linear model
T_lin_s = sum(irf_CImeanLIN(s_desc.KPIF, 1:2)) ;
N_lin_s = sum(irf_CImeanLIN(s_desc.KIX, 1:2)); %kanske ändra för att shocken inte kommer in direkt i t

T_lin_m = sum(irf_CImeanLIN(s_desc.KPIF,1:4)) ;
N_lin_m = sum(irf_CImeanLIN(s_desc.KIX, 1:4));

T_lin_l = sum(irf_CImeanLIN(s_desc.KPIF,:)) ;
N_lin_l = sum(irf_CImeanLIN(s_desc.KIX,:));

ERPT_lin_s = T_lin_s/N_lin_s
ERPT_lin_m = T_lin_m/N_lin_m
ERPT_lin_l = T_lin_l/N_lin_l

%High inflation
T_high_s = sum(irf_CImeanEXP(s_desc.KPIF, 1:2)) ;
N_high_s = sum(irf_CImeanEXP(s_desc.KIX, 1:2));

T_high_m = sum(irf_CImeanEXP(s_desc.KPIF,1:4)) ;
N_high_m = sum(irf_CImeanEXP(s_desc.KIX,1:4));

T_high_l = sum(irf_CImeanEXP(s_desc.KPIF,:)) ;
N_high_l = sum(irf_CImeanEXP(s_desc.KIX,:));

ERPT_high_s = T_high_s/N_high_s
ERPT_high_m = T_high_m/N_high_m
ERPT_high_l = T_high_l/N_high_l



%Low inflation
T_low_s = sum(irf_CImeanREC(s_desc.KPIF, 1:2)) ;
N_low_s = sum(irf_CImeanREC(s_desc.KIX, 1:2));

T_low_m = sum(irf_CImeanREC(s_desc.KPIF,1:4)) ;
N_low_m = sum(irf_CImeanREC(s_desc.KIX,1:4));

T_low_l = sum(irf_CImeanREC(s_desc.KPIF,:)) ;
N_low_l = sum(irf_CImeanREC(s_desc.KIX,:));

ERPT_low_s = T_low_s/N_low_s
ERPT_low_m = T_low_m/N_low_m
ERPT_low_l = T_low_l/N_low_l
