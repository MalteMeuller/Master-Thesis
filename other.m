clear all ;
close all ;


%%  better and mroe detailed table
for i=["KPIF", "PPI", "IMPI"]

    name=sprintf(i)
    filename = "C:\Users\malte\OneDrive - Handelshögskolan i Stockholm\Thesis\MATLAB\robust_order_1\code_fullsample_"+name+".mat";
    load (filename)

    if strcmpi(name,"KPIF")
        plot_num=1; 
    elseif strcmpi(name,"PPI")
        plot_num=5;
    else %impi
        plot_num=9;
    end 

    disp("MAR")
    mean(acceptrate)
    percentile=percentile/100 ;
    
    %high and low get from different distributions, here we say=1=rec
    %ändrade här till final final omegamatrix
    [irf_DimpREC,irf_CImeanREC,CIupREC,CIlowREC,IRFseREC]=irf_VAR_MCMC_struct(beta1E,Omega1E,Omega1mat,beta1matF,s_desc);
    [irf_DimpEXP,irf_CImeanEXP,CIupEXP,CIlowEXP,IRFseEXP]=irf_VAR_MCMC_struct(beta0E,Omega0E,Omega0mat,beta0matF,s_desc);
    
    %linear model, köra threshold (60/40) eller 50/50?
    percentile = 0.5 ;
    beta0E = beta1E * percentile + beta0E*(1-percentile);
    Omega0E =Omega1E * percentile + Omega0E*(1-percentile);
    Omega0mat =Omega1mat * percentile + Omega0mat*(1-percentile);
    beta0matF = beta1matF * percentile + beta0matF*(1-percentile);
    [irf_DimpLIN,irf_CImeanLIN,CIupLIN,CIlowLIN,IRFseLIN]=irf_VAR_MCMC_struct(beta0E,Omega0E,Omega0mat,beta0matF,s_desc);

    %linear
    PRICE.lin_short = sum(irf_CImeanLIN(s_desc.KPIF, 1:2)) ;
    NEER.lin_short = sum(irf_CImeanLIN(s_desc.KIX, 1:2)) ;
    PERR.lin_short = PRICE.lin_short/NEER.lin_short
    
    PRICE.lin_medium = sum(irf_CImeanLIN(s_desc.KPIF, 1:4)) ;
    NEER.lin_medium = sum(irf_CImeanLIN(s_desc.KIX, 1:4)) ;
    PERR.lin_medium = PRICE.lin_medium/NEER.lin_medium
    
    PRICE.lin_long = sum(irf_CImeanLIN(s_desc.KPIF, 1:20)) ;
    NEER.lin_long = sum(irf_CImeanLIN(s_desc.KIX, 1:20)) ;
    PERR.lin_long = PRICE.lin_long/NEER.lin_long ;
    
    %High inflation
    PRICE.high_short = sum(irf_CImeanEXP(s_desc.KPIF, 1:2)) ;
    NEER.high_short = sum(irf_CImeanEXP(s_desc.KIX, 1:2)) ;
    PERR.high_short = PRICE.high_short/NEER.high_short ;
    
    PRICE.high_medium = sum(irf_CImeanEXP(s_desc.KPIF, 1:4)) ;
    NEER.high_medium = sum(irf_CImeanEXP(s_desc.KIX, 1:4)) ;
    PERR.high_medium = PRICE.high_medium/NEER.high_medium ;
    
    PRICE.high_long = sum(irf_CImeanEXP(s_desc.KPIF, 1:20)) ;
    NEER.high_long = sum(irf_CImeanEXP(s_desc.KIX, 1:20)) ;
    PERR.high_long = PRICE.high_long/NEER.high_long ;
    
    %Low inflation
    PRICE.low_short = sum(irf_CImeanREC(s_desc.KPIF, 1:2)) ;
    NEER.low_short = sum(irf_CImeanREC(s_desc.KIX, 1:2)) ;
    PERR.low_short = PRICE.lin_short/NEER.low_short ;
    
    PRICE.low_medium = sum(irf_CImeanREC(s_desc.KPIF, 1:4)) ;
    NEER.low_medium = sum(irf_CImeanREC(s_desc.KIX, 1:4)) ;
    PERR.low_medium = PRICE.low_medium/NEER.low_medium;
    
    PRICE.low_long = sum(irf_CImeanREC(s_desc.KPIF, 1:20)) ;
    NEER.low_long = sum(irf_CImeanREC(s_desc.KIX, 1:20)) ;
    PERR.low_long = PRICE.low_long/NEER.low_long;
    
    N = struct2table(NEER);
    N = rows2vars(N);
    N.Properties.VariableNames = ["Model","NEER"];
    
    P = struct2table(PRICE);
    P = rows2vars(P);
    P.Properties.VariableNames = ["Model","ERPT"];
    
    PT = struct2table(PERR);
    PT = rows2vars(PT);
    PT.Properties.VariableNames = ["Model","PERR"];
    
    Final = join(join(N, P, 'Keys', 'Model'), PT, 'Keys', 'Model');
    
    outputname = "C:\Users\malte\OneDrive - Handelshögskolan i Stockholm\Thesis\MATLAB\robust_noinf\table_"+name+".xlsx";
    writetable(Final, outputname);



    cum_lin=[];
    cum_high = [];
    cum_low = [];
    
    
    for i=[1:20]
        lin  = sum(irf_CImeanLIN(s_desc.KPIF, 1:i)) ;
        high = sum(irf_CImeanEXP(s_desc.KPIF, 1:i)) ;
        low  = sum(irf_CImeanREC(s_desc.KPIF, 1:i)) ;
        cum_lin(i)=lin ;
        cum_high(i)=high ;
        cum_low(i)=low ;
    
    
    end


    if strcmpi(name,"KPIF")
        name="CPIF";
    end 


    plotname = "Cumulative ERPT "+name+""
    % Plot the cumulative values
    x = 1:20;
    figure(plot_num);
    plot(x, cum_lin, 'black', x, cum_high, 'red', x, cum_low, 'blue', 'LineWidth', 2);
    set(gca, 'FontSize', 12);
    xlabel('Time');
    ylabel('');
    title(plotname);
    legend('Linear', 'High', 'Low', 'FontSize',12, 'Location', 'southeast');
    %ylim([0,1])
    grid on;

    
    


%% plots
    titlename = "NEER shock => "+name+" response"

    figure('Position',[100 100 600 325])
    figure(plot_num+1);
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupLIN(s_desc.KPIF,:),fliplr(CIlowLIN(s_desc.KPIF,:))]; %which shocks to use CI on
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 s_desc.irf_hor], [0 0], 'k', 'LineWidth', 1)
    hold off
    title(titlename)
    set(gca, 'FontSize', 12);
    xlim([1 s_desc.irf_hor])
    %ylim([-0.1 0.5])
    % Increase the font size of the legend
    legend('90% CI LIN', 'Linear', 'High', 'Low', 'FontSize', 12)

    %CI high
    figure('Position',[100 100 600 325])
    figure(plot_num+2);
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupEXP(s_desc.KPIF,:),fliplr(CIlowEXP(s_desc.KPIF,:))];
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 20], [0 0], 'k', 'LineWidth', 1)
    hold off
    title(titlename)
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.1 0.5])
    % Increase the font size of the legend
    legend('90% CI HIGH', 'Linear', 'High', 'Low', 'FontSize', 12)


    % low inflation  
    figure('Position',[100 100 600 325])
    figure(plot_num+3);
    xv = [1:s_desc.irf_hor,fliplr(1:s_desc.irf_hor)]; %lenght of the shocks
    yv = [CIupREC(s_desc.KPIF,:),fliplr(CIlowREC(s_desc.KPIF,:))];
    hReg = fill(xv,yv,[0.75 0.75 0.75],'EdgeColor','none'); % draw region
    hold on
    plot(1:s_desc.irf_hor,irf_CImeanLIN(s_desc.KPIF,:),'black-o','MarkerFaceColor',[1 1 1],'Linewidth',2)
    plot(1:s_desc.irf_hor,irf_CImeanEXP(s_desc.KPIF,:),'red--','Linewidth',3)
    plot(1:s_desc.irf_hor,irf_CImeanREC(s_desc.KPIF,:),'blue:','Linewidth',3)
    plot([1 20], [0 0], 'k', 'LineWidth', 1)
    hold off
    title(titlename)
    % Increase the font size for the Y-axis labels
    set(gca, 'FontSize', 12);
    xlim([1 20])
    %ylim([-0.1 0.5])
    % Increase the font size of the legend
    legend('90% CI LOW', 'Linear', 'High', 'Low', 'FontSize', 12)

  
end   