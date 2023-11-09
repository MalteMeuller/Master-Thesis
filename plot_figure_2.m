% this code plots figure 2 based on the IRF used in irf_VAR_MCMC_STRUCT
%NOTATIONS

%wierd things going on with the notations
clear all
close all
test=0;  %change here
notuse=0;
normal=1;
dyn=0;

%%
if normal ==1
    load('C:\Users\malte\OneDrive - Handelshögskolan i Stockholm\Thesis\MATLAB\LST-SVAR-QUARTER\KPIF\code_fullsample_KPIF.mat')
    disp("MAR")
    mean(acceptrate)
    percentile=percentile/100 ;
    
     
    %discard burn ins and sums up values for all draws
    Omega0Ex1=zeros(size(Omega0E));
    Omega1Ex1=zeros(size(Omega0E));
    for ii=1:length(Amat(s_MCMC.dropfirstT:end-1,:))
        [Omega0Ex,Omega1Ex,thetaEx]=vec2matScov(Amat(s_MCMC.dropfirstT+ii-1,:)',s_desc.Omega_length_cov);
        Omega0Ex1=Omega0Ex1+Omega0Ex;
        Omega1Ex1=Omega1Ex1+Omega1Ex;
    end
    
    %takes the average of all 80 000 draws, parameter values. 
    Omega1E=Omega1Ex1/length(Amat(s_MCMC.dropfirstT:end-1,:));
    Omega0E=Omega0Ex1/length(Amat(s_MCMC.dropfirstT:end-1,:));
    
    %high and low get from different distributions
    [irf_DimpREC,irf_CImeanREC,CIupREC,CIlowREC,IRFseREC]=irf_VAR_MCMC_struct(beta1E,Omega1E,Omega1mat,beta1matF,s_desc);
    [irf_DimpEXP,irf_CImeanEXP,CIupEXP,CIlowEXP,IRFseEXP]=irf_VAR_MCMC_struct(beta0E,Omega0E,Omega0mat,beta0matF,s_desc);
    
    %linear model, köra threshold 50/50
    percentile = 0.5 ;
    beta0E = beta1E * percentile + beta0E*(1-percentile);
    Omega0E =Omega1E * percentile + Omega0E*(1-percentile);
    Omega0mat =Omega1mat * percentile + Omega0mat*(1-percentile);
    beta0matF = beta1matF * percentile + beta0matF*(1-percentile);
    [irf_DimpLIN,irf_CImeanLIN,CIupLIN,CIlowLIN,IRFseLIN]=irf_VAR_MCMC_struct(beta0E,Omega0E,Omega0mat,beta0matF,s_desc);

    %=====================================================================
    %             Impulse responses in the linear model (CI lin)
    %=====================================================================
    figure('Position',[100 100 600 325])
    figure(1)
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupLIN(s_desc.KPIF,:),fliplr(CIlowLIN(s_desc.KPIF,:))]; %which shocks to use CI on
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 s_desc.irf_hor], [0 0], 'k', 'LineWidth', 1)
    hold off
    title('NEER shock => CPIF response')
    set(gca, 'FontSize', 12);
    xlim([1 s_desc.irf_hor])
    %ylim([-0.1 0.5])
    % Increase the font size of the legend
    legend('90% CI LIN', 'Linear', 'High', 'Low', 'FontSize', 12)

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
    title('NEER shock => CPIF response')
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.1 0.5])
    % Increase the font size of the legend
    legend('90% CI HIGH', 'Linear', 'High', 'Low', 'FontSize', 12)


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
    title('NEER shock => CPIF response')
    % Increase the font size for the Y-axis labels
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.1 0.5])
    % Increase the font size of the legend
    legend('90% CI LOW', 'Linear', 'High', 'Low', 'FontSize', 12)

   
    %ERPT
    PRICE.lin_short = sum(irf_CImeanLIN(s_desc.KPIF, 1:2)) ;
    NEER.lin_short = sum(irf_CImeanLIN(s_desc.KIX, 1:2)) ;
    ERPT.lin_short = PRICE.lin_short/NEER.lin_short
    
    PRICE.lin_medium = sum(irf_CImeanLIN(s_desc.KPIF, 1:4)) ;
    NEER.lin_medium = sum(irf_CImeanLIN(s_desc.KIX, 1:4)) ;
    ERPT.lin_medium = PRICE.lin_medium/NEER.lin_medium
    
    PRICE.lin_long = sum(irf_CImeanLIN(s_desc.KPIF, 1:20)) ;
    NEER.lin_long = sum(irf_CImeanLIN(s_desc.KIX, 1:20)) ;
    ERPT.lin_long = PRICE.lin_long/NEER.lin_long ;
    
    %High inflation
    PRICE.high_short = sum(irf_CImeanEXP(s_desc.KPIF, 1:2)) ;
    NEER.high_short = sum(irf_CImeanEXP(s_desc.KIX, 1:2)) ;
    ERPT.high_short = PRICE.high_short/NEER.high_short ;
    
    PRICE.high_medium = sum(irf_CImeanEXP(s_desc.KPIF, 1:4)) ;
    NEER.high_medium = sum(irf_CImeanEXP(s_desc.KIX, 1:4)) ;
    ERPT.high_medium = PRICE.high_medium/NEER.high_medium ;
    
    PRICE.high_long = sum(irf_CImeanEXP(s_desc.KPIF, 1:20)) ;
    NEER.high_long = sum(irf_CImeanEXP(s_desc.KIX, 1:20)) ;
    ERPT.high_long = PRICE.high_long/NEER.high_long ;
    
    %Low inflation
    PRICE.low_short = sum(irf_CImeanREC(s_desc.KPIF, 1:2)) ;
    NEER.low_short = sum(irf_CImeanREC(s_desc.KIX, 1:2)) ;
    ERPT.low_short = PRICE.lin_short/NEER.low_short ;
    
    PRICE.low_medium = sum(irf_CImeanREC(s_desc.KPIF, 1:4)) ;
    NEER.low_medium = sum(irf_CImeanREC(s_desc.KIX, 1:4)) ;
    ERPT.low_medium = PRICE.low_medium/NEER.low_medium;
    
    PRICE.low_long = sum(irf_CImeanREC(s_desc.KPIF, 1:20)) ;
    NEER.low_long = sum(irf_CImeanREC(s_desc.KIX, 1:20)) ;
    ERPT.low_long = PRICE.low_long/NEER.low_long;
    
    N = struct2table(NEER);
    N = rows2vars(N);
    N.Properties.VariableNames = ["Model","NEER"];
    
    P = struct2table(PRICE);
    P = rows2vars(P);
    P.Properties.VariableNames = ["Model","PRICE"];
    
    PT = struct2table(ERPT);
    PT = rows2vars(PT);
    PT.Properties.VariableNames = ["Model","ERPT"];
    
    Final = join(join(N, P, 'Keys', 'Model'), PT, 'Keys', 'Model')

%strA=['code_fullsample_KPIF.mat'];
%save(strA)

end



%% Dynamic feedback

%=====================================================================
%             Construct the dynamic multiplier with endo
%             exit from a regime
%=====================================================================


if dyn==1
    load('C:\Users\malte\Dropbox\Thesis\MATLAB\LST-SVAR-QUARTER\KPIF\code_fullsample_KPIF.mat')
    mean(acceptrate)
    percentile=percentile/100 ;

    clear CIlowREC CIupREC CIlowEXP CIupEXP CIlowLIN CIupLIN irf_CImeanREC irf_CImeanEXP irf_CImeanLIN
    
    %setting up three different starting values
    s_desc.dyn =  0.8
    s_desc.meanrev= mean(Z0)
    s_desc.const = mean(Z0)*(1-s_desc.dyn)
    %dynamics for how fast it returns to 0
    high=1
    low=-1
    lin=0
    


     %give differernt starting values for different regimes
     s_desc.startZirf = high; %high
    [irf_Dimp_dyn0_max,irf_CImeanEXP,CIupEXP,CIlowEXP,IRFse_dyn0_max,Zser000mat0_max,F_Z_tx000mat0_max]=...
        irf_VAR_MCMC_struct_alt_dynamic_feedback3(beta0E,Omega0E,beta1E,Omega1E,...
        Omega0mat,beta0matF,Omega1mat,beta1matF,s_desc,s_prior);
    
    s_desc.startZirf = low; %low
    [irf_Dimp_dyn0_max,irf_CImeanREC,CIupREC,CIlowREC,IRFse_dyn0_max,Zser000mat0_max,F_Z_tx000mat0_max]=...
        irf_VAR_MCMC_struct_alt_dynamic_feedback3(beta0E,Omega0E,beta1E,Omega1E,...
        Omega0mat,beta0matF,Omega1mat,beta1matF,s_desc,s_prior);


    s_desc.startZirf = lin; %linear
    [irf_Dimp_dyn0_max,irf_CImeanLIN,CIupLIN,CIlowLIN,IRFse_dyn0_max,Zser000mat0_max,F_Z_tx000mat0_max]=...
    irf_VAR_MCMC_struct_alt_dynamic_feedback3(beta0E,Omega0E,beta1E,Omega1E,...
    Omega0mat,beta0matF,Omega1mat,beta1matF,s_desc,s_prior);

    

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
    title('NEER shock => CPIF response')
   % Increase the font size for the Y-axis labels
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.06 0.08])
    % Increase the font size of the legend
    legend('90% CI LIN', 'Linear', 'High', 'Low', 'FontSize', 12)
   
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
    title('NEER shock => CPIF response')
    % Increase the font size for the Y-axis labels
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.06 0.08])
    % Increase the font size of the legend
    legend('90% CI HIGH', 'Linear', 'High', 'Low', 'FontSize', 12)

    % low inflation  
    figure('Position', [100 100 600 325])
    figure(3)
    xv = [1:s_desc.irf_hor, fliplr(1:s_desc.irf_hor)]; % length of the shocks
    yv = [CIupREC(s_desc.KPIF, :), fliplr(CIlowREC(s_desc.KPIF, :))];
    hReg = fill(xv, yv, [0.75 0.75 0.75], 'EdgeColor', 'none'); % draw region
    hold on
    plot(1:s_desc.irf_hor, irf_CImeanLIN(s_desc.KPIF, :), 'black-o', 'MarkerFaceColor', [1 1 1], 'Linewidth', 2)
    plot(1:s_desc.irf_hor, irf_CImeanEXP(s_desc.KPIF, :), 'red--', 'Linewidth', 3)
    plot(1:s_desc.irf_hor, irf_CImeanREC(s_desc.KPIF, :), 'blue:', 'Linewidth', 3)
    plot([1 20], [0 0], 'k', 'LineWidth', 1)
    hold off
    title('NEER shock => CPIF response')
    % Increase the font size for the Y-axis labels
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.06 0.08])
    % Increase the font size of the legend
    legend('90% CI LOW', 'Linear', 'High', 'Low', 'FontSize', 12)


    %ERPT
    PRICE.lin_short = sum(irf_CImeanLIN(s_desc.KPIF, 1:2)) ;
    NEER.lin_short = sum(irf_CImeanLIN(s_desc.KIX, 1:2)) ;
    ERPT.lin_short = PRICE.lin_short/NEER.lin_short
    
    PRICE.lin_medium = sum(irf_CImeanLIN(s_desc.KPIF, 1:4)) ;
    NEER.lin_medium = sum(irf_CImeanLIN(s_desc.KIX, 1:4)) ;
    ERPT.lin_medium = PRICE.lin_medium/NEER.lin_medium
    
    PRICE.lin_long = sum(irf_CImeanLIN(s_desc.KPIF, 1:20)) ;
    NEER.lin_long = sum(irf_CImeanLIN(s_desc.KIX, 1:20)) ;
    ERPT.lin_long = PRICE.lin_long/NEER.lin_long ;
    
    %High inflation
    PRICE.high_short = sum(irf_CImeanEXP(s_desc.KPIF, 1:2)) ;
    NEER.high_short = sum(irf_CImeanEXP(s_desc.KIX, 1:2)) ;
    ERPT.high_short = PRICE.high_short/NEER.high_short ;
    
    PRICE.high_medium = sum(irf_CImeanEXP(s_desc.KPIF, 1:4)) ;
    NEER.high_medium = sum(irf_CImeanEXP(s_desc.KIX, 1:4)) ;
    ERPT.high_medium = PRICE.high_medium/NEER.high_medium ;
    
    PRICE.high_long = sum(irf_CImeanEXP(s_desc.KPIF, 1:20)) ;
    NEER.high_long = sum(irf_CImeanEXP(s_desc.KIX, 1:20)) ;
    ERPT.high_long = PRICE.high_long/NEER.high_long ;
    
    %Low inflation
    PRICE.low_short = sum(irf_CImeanREC(s_desc.KPIF, 1:2)) ;
    NEER.low_short = sum(irf_CImeanREC(s_desc.KIX, 1:2)) ;
    ERPT.low_short = PRICE.lin_short/NEER.low_short ;
    
    PRICE.low_medium = sum(irf_CImeanREC(s_desc.KPIF, 1:4)) ;
    NEER.low_medium = sum(irf_CImeanREC(s_desc.KIX, 1:4)) ;
    ERPT.low_medium = PRICE.low_medium/NEER.low_medium;
    
    PRICE.low_long = sum(irf_CImeanREC(s_desc.KPIF, 1:20)) ;
    NEER.low_long = sum(irf_CImeanREC(s_desc.KIX, 1:20)) ;
    ERPT.low_long = PRICE.low_long/NEER.low_long;
    
    N = struct2table(NEER);
    N = rows2vars(N);
    N.Properties.VariableNames = ["Model","NEER"];
    
    P = struct2table(PRICE);
    P = rows2vars(P);
    P.Properties.VariableNames = ["Model","PRICE"];
    
    PT = struct2table(ERPT);
    PT = rows2vars(PT);
    PT.Properties.VariableNames = ["Model","ERPT"];
    
    Final = join(join(N, P, 'Keys', 'Model'), PT, 'Keys', 'Model')

end

 