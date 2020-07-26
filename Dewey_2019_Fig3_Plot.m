%% Plot Figure 3 from Dewey et al. 2019 J Neurosci paper
% loads 'SUPPRESSION_DATA_FIG3.mat', which contains data from a representative mouse.

% The 'data' structure contains the following fields:
% data.cf; % CF (Hz)
% data.supCrits; % suppression criterion (dB) used for analysis

% and the following substructures:

% data.BM_CF_L60; % suppression data for BM responses to 60 dB SPL, ~CF tones
% data.RL_CF_L60; % suppression data for RL responses to 60 dB SPL, ~CF tones
% data.RL_4kHz_L60; % suppression data for RL responses to 60 dB SPL, ~4 kHz tones

% Each substructure includes the following fields:
%         f1: 9600 % probe frequency (Hz)
%         L1: 60 % probe level (dB SPL)
%        f2s: [1×20 double] % suppressor frequencies (Hz)
%        L2s: [30 35 40 45 50 55 60 65 70 75 80 85 90] % suppressor levels (dB SPL)

% and the following substructures:
%       vib1: [1×1 struct] % displacement data for probe alone 
%       vib2: [1×1 struct] % displacement data for probe + suppressor
%       vib3: [1×1 struct] % displacement data for suppressor alone

% Each vibration substructure includes displacement data for 'f1' (probe) and 'f2' (suppressor), 
% w/ each substructure (e.g., 'data.BM_CF_L60.vib1.f1') containing:
% 
%      mag: [20×13 double] % displacement magnitude (nm)
%       nf: [20×13 double] % mean noise floor (nm) 
%     nfsd: [20×13 double] % noise floor standard deviation (nm)
%      phi: [20×13 double] % displacement phase (radians)

%  with data organized by suppressor frequency x suppressor level (20x13)

% Suppression analysis is also found here: e.g., 'data.BM_CF_L60.supNal', which contains:
%
%     supThreshL: [20×8 double] % suppression thresholds (dB SPL)
%        supDisp: [20×8 double] % suppressor-evoked displacement at the suppression threshold (nm)
%          supNF: [20×8 double] % noise floor (mean + 3 SDs) for the suppressor-evoked displacement at suppression threshold
%     probePhase: [20×8 double] % suppressor-evoked change in phase of the response to the probe at suppression threshold

% Data are organized by suppressor frequency x suppression criterion (20x8)
% Suppression criteria are found in data.supCrits (1.5 - 12 dB in 1.5 dB steps).

clear;clc;close all;
load('SUPPRESSION_DATA_FIG3.mat');
cf = data.cf; % CF (Hz)

%% Plot Fig3A-C: Probe displacement vs. suppressor level
figNames = [{'Fig3A'}; {'Fig3B'}; {'Fig3C'}];
figN = length(figNames);

f2s = [3500 9500 14500 20500]; % suppressor frequencies to plot (Hz) 
f2N = length(f2s); % # suppressor frequencies
L2s = 30:5:90; % suppressor levels (dB SPL)
L2N = length(L2s); % # suppressor levels

%% Plotting parameters
pltH = 0.36; % plot height
lStys = [{'-.'};{'--'};{':'};{'x:'}]'; % line styles for each suppressor frequency
axlw = 1.8; % axis line width
fntSz = 16; % font size 
prbClr = [0.4 0.4 0.4]; % color for dashed line at probe level

%% Plot each figure panel
for fig_i = 1:figN
    figName = figNames{fig_i};
    panel = figName(end);
    
    switch panel
        case 'A'
            cond = 'BM_CF_L60';
            clr = 'k';
            plotLegend = 1; % plot legend? 0 = no; 1 = yes
        case 'B'
            cond = 'RL_CF_L60';
            clr = 'r';
            plotLegend = 0; % plot legend? 0 = no; 1 = yes
        case 'C'
            cond = 'RL_4kHz_L60';
            clr = [0 170 200]/255;
            plotLegend = 0; % plot legend? 0 = no; 1 = yes
    end
    
    d = data.(genvarname(cond)); % get data
    
    %% Calculate average unsuppressed probe magnitude
    f1mags = nan(1,f2N);
    for f2_i = 1:f2N
        f2 = f2s(f2_i);
        [~,f2_ix] = ismember(f2, d.f2s); % f2 index in data set     
        f1mags(f2_i)  = mean(d.vib1.f1.mag(f2_ix,:)); % mean unsuppressed probe response (nm)
    end
    f1mag_ave = mean(f1mags); % average unsuppressed probe response (nm)
    
    %% SET UP FIGURE
    h = figure('units','normalized','position',[0 0 .28 1]);
    ax1 = axes('Position',[.1 .66 .8 .27],'box','off','FontSize',fntSz,'LineWidth',axlw); hold all;
    xlim([27.5 92.5]); set(gca,'Xtick',30:10:90); ax1.XRuler.TickLength = [0.025 0.025];
    ylim([.1 20]); set(gca,'Yscale','log','Ytick',[0.1 1 10],'YtickLabels',[0.1 1 10]); ax1.YRuler.TickLength = [0.025 0.025]; ax1.YAxis.Exponent = 0; hold on;    
    plot([L2s(1) L2s(end)], [f1mag_ave f1mag_ave], '-', 'LineWidth', 1.1, 'Color', prbClr); hold on; % plot average unsuppressed probe response (nm)

    % Plot suppression by each specific suppressor frequency
    for f2_i = 1:f2N
        f2 = f2s(f2_i); % suppressor frequency
        
        [~,f2_ix] = ismember(f2, d.f2s); % f2 index in data set     
        lSty = lStys{f2_i};
        
        if f2_i < 4 % set line width
            lw = 2.5;
            mrkSz = 5;
        else
            lw = 1.5;
            mrkSz = 8;
        end
        
        f1mag1 = d.vib1.f1.mag(f2_ix,:); % unsuppressed probe response (nm)
        f1mag2 = d.vib2.f1.mag(f2_ix,:); % suppressed probe response (nm)
        f1mag2 = (f1mag2 ./ f1mag1) * f1mag_ave; % re-scale suppressed response relative to mean unsuppressed response
    
        %% PLOT
        axes(ax1);
        plot(L2s, f1mag2, lSty, 'MarkerSize', mrkSz, 'LineWidth', lw, 'Color', clr);
    end
    
    if plotLegend % Legend
        lgd = legend(' ',' ',' ',' ',' ');
        lgd.Location = 'southwest';
        legend('boxoff');
        lgd.Position = [0.18 0.73 0.2 0.1];
    end
    
    yyaxis right; % suppression (dB) axis
    ylim(20*log10([.1 20]/f1mag_ave)); ax1.YColor = 'k'; ax1.YRuler.TickLength = [0.025 0.025];
end


%% Plot Fig3D-F: Suppression thresholds and corresponding suppressor-evoked displacements
figNames = [{'Fig3D'}; {'Fig3E'}; {'Fig3F'}];
figN = length(figNames);

%% Plotting parameters
pltH = 0.32; % plot height
axlw = 1.5; % 
fntSz = 16; % font size
prbClr = [0.4 0.4 0.4];
xlims = [1.4 20];

for fig_i = 1:figN
    figName = figNames{fig_i};
    panel = figName(end);
    
    switch panel
        case 'D'
            cond = 'BM_CF_L60';
            clr = 'k';
            supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
            supCritN = length(supCrits); % # suppression criteria
            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotSupDisp = 1; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe phase at suppression threshold? 0 = no; 1 = yes
            L90x = 12.5; % frequency where suppressor-evoked displacement falls below noise floor (kHz)
        case 'E'
            cond = 'RL_CF_L60';
            clr = 'r';
            supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
            supCritN = length(supCrits); % # suppression criteria
            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotSupDisp = 1; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe phase at suppression threshold? 0 = no; 1 = yes
            L90x = 13.5; % frequency where suppressor-evoked displacement falls below noise floor (kHz)
        case 'F'
            cond = 'RL_4kHz_L60';
            clr = [0 170 200]/255;   
            supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
            supCritN = length(supCrits); % # suppression criteria
            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotSupDisp = 1; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe phase at suppression threshold? 0 = no; 1 = yes
            L90x = 13.5; % frequency where suppressor-evoked displacement falls below noise floor (kHz)
    end
    
    d = data.(genvarname(cond)); % get data for that condition
    f1 = d.f1; % probe freq. (Hz)
    L1 = d.L1; % probe level (dB SPL)
    f1mag_ave = mean(mean(d.vib1.f1.mag)); % average unsuppressed probe response (nm)
    sup = d.supNal; % get suppression analysis
    
    %% SET UP FIGURE
    h = figure('units','normalized','position',[0 0 .28 1]);
    fillClr = [0.9 0.9 0.9];

    if plotSupThresh %% Suppression thresholds
        ax1 = axes('Position',[.12 .66 .8 pltH],'box','off', 'LineWidth', axlw, 'FontSize', fntSz,'Layer','top');
        hfill = rectangle('Position',[L90x,10,20,100],'FaceColor',fillClr,'EdgeColor',fillClr,'LineWidth',.1); hold on;
        ax1;
        
        plot([cf/1000 cf/1000], [0 100], 'k:', 'MarkerSize', 5, 'Linewidth', 1.2); hold on;
        plot([.1 100], [L1 L1], 'k--', 'LineWidth',1);
        plot(f1/1000, L1, 'ko','MarkerSize', mrkSz, 'LineWidth', 1.5);
        plot(cf/1000,34, 'k^','MarkerSize', mrkSz, 'Linewidth',1.5, 'MarkerFaceColor','k'); hold on;
        xlim(xlims); set(gca,'Xscale','log','Xtick',[0.5 1 2 5 10 15 20]); ax1.XRuler.TickLength = [0.025 0.025];
        ylim([30 90]); ax1.YRuler.TickLength = [0.025 0.025];
    end

    if plotSupDisp % Suppressor evoked displacement at suppression threshold
        ax2 = axes('Position',[.12 .305 .8 pltH],'box', 'off', 'LineWidth', axlw, 'FontSize', fntSz, 'Layer','top');
        hfill = rectangle('Position',[L90x,.000001,20,1000],'FaceColor',fillClr,'EdgeColor',fillClr,'LineWidth',.1); hold on;
        ax2;
        xlim(xlims); set(gca,'Xscale','log','Xtick',[0.5 1 2 5 10 15 20]); ax2.XRuler.TickLength = [0.025 0.025];
        ylim([0.07 70]); set(gca,'Yscale','log','Ytick',[0.1 1 10 50]); ax2.YRuler.TickLength = [0.025 0.025];

        plot([.1 100], [f1mag_ave f1mag_ave], 'k--', 'LineWidth', 1);  
        plot([cf/1000 cf/1000], [0.001 1000], 'k:','MarkerSize',5, 'Linewidth',1.2); hold on;
        plot(cf/1000, 0.1, 'k^','MarkerSize', mrkSz, 'Linewidth',.5, 'MarkerFaceColor','k'); hold on;
    end

    if plotProbePhase % Suppressor-evoked change in probe phase at suppression threshold
        ax3 = axes('Position',[.12 .05 .8 .22],'box','off', 'LineWidth', axlw, 'FontSize', fntSz,'Layer','top'); 
        hfill = rectangle('Position',[L90x,-20,20,1000],'FaceColor',fillClr,'EdgeColor',fillClr,'LineWidth',.1); hold on;
        ax3;
        
        xlim(xlims); set(gca,'Xscale','log','Xtick',[0.5 1 2 5 10 15 20]); ax3.XRuler.TickLength = [0.025 0.025];
        ylim([-0.2 0.2]); set(gca,'Ytick',-.3:.1:.3); ax3.YRuler.MinorTick = 'on'; ax3.YRuler.MinorTickValues = [-3:0.1:3]; ax3.YRuler.TickLength = [0.025 0.025];

        plot([cf/1000 cf/1000], [-5 5], 'k:','Linewidth',1.2); hold on;
        plot([0.5 100], [0 0], 'k--','LineWidth', 1); hold on;
        plot(cf/1000,-0.15, 'k^','MarkerSize', mrkSz, 'Linewidth',.5,'MarkerFaceColor','k'); hold on;
        plot(f1/1000, 0, 'k-','MarkerSize', mrkSz,'LineWidth',1.5);
    end
    
    %% Plot data for each suppression criterion
    for supCrit_i = 1:supCritN
        lw = 0.5+supCrit_i*0.25;
        
        %% Suppression threshold (dB SPL)
        if plotSupThresh
            supThreshL = sup.supThreshL(:, supCrit_i);
            axes(ax1);
            plot(d.f2s/1000, supThreshL, '-','LineWidth', lw, 'Color', clr); hold on;
        end
        
        %% Suppressor-evoked displacement at suppression threshold (nm)
        if plotSupDisp
            supDisp = sup.supDisp(:, supCrit_i);
            supNF = sup.supNF(:,supCrit_i);
            supDisp(supDisp<supNF) = NaN;
            
            axes(ax2);
            plot(d.f2s/1000, supDisp, '-','LineWidth', lw, 'Color', clr); hold on;
        end
        
        %% Suppressor-evoked change in probe phase at suppression threshold (cycles)
        if plotProbePhase
            probePhase =  sup.probePhase(:, supCrit_i);            
            axes(ax3);
            plot(d.f2s/1000, unwrap(probePhase)/(2*pi), '-','LineWidth', lw, 'Color', clr); hold on;
        end
    end
end