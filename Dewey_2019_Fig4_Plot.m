%% Plot Figure 4 from Dewey et al. 2019 J Neurosci paper
% Loads 'SUPPRESSION_DATA.mat' containing the 'allData' structure, with the following fields:

%      supF_re_CF: [1×101 double] % suppressor frequency re CF (octaves)
%        supCrits: [1.5000 3 4.5000 6 7.5000 9 10.5000 12] % suppression criteria (dB)
%             cfs: [9 9.5000 9.5000 9 9.5000 9.5000 9.5000 9 9 9.5000] % CFs (kHz)

% and these substructures containing suppression threshold data for either 'BM' or 'RL', 
% for different probe frequencies ('CF', '6kHz', etc.) and probe levels (e.g.,'L40' = 40 dB SPL):

%       BM_CF_L40: [1×1 struct]
%       BM_CF_L50: [1×1 struct]
%       BM_CF_L60: [1×1 struct]
%       RL_CF_L40: [1×1 struct]
%       RL_CF_L50: [1×1 struct]
%       RL_CF_L60: [1×1 struct]
%     RL_6kHz_L60: [1×1 struct]
%     RL_5kHz_L60: [1×1 struct]
%     RL_4kHz_L60: [1×1 struct]
%     RL_4kHz_L50: [1×1 struct]
%     RL_3kHz_L60: [1×1 struct]
%     RL_2kHz_L60: [1×1 struct]

%  Each suppression data substructure contains the following fields:
%     probeF_re_CFs: ... % probe frequency re CF (octaves)
%         probeDisp: [1×10 double] % unsuppressed probe displacement (nm)
%        supThreshL: [101×8×10 double] % suppression threshold (dB SPL)
%           supDisp: [101×8×10 double] % suppressor-evoked displacement at suppression threshold (nm)
%             supNF: [101×8×10 double] % noise floor (mean + 3SD) for suppressor-evoked displacement at suppression threshold (nm)
%     supDisp_Clean: [101×8×10 double] % "clean" suppressor-evoked displacement at suppression threshold (nm)
%        probePhase: [101×8×10 double] % suppressor-evoked change in phase at the probe frequency at suppression threshold (cycles)

% Data are organized by suppressor frequency x suppression criterion x mouse (101x8x10)

% Also loaded is the 'SUPPRESSION_PARAMS.mat' file which includes 'allParams', a structure that contains results of analyses on the
% suppression tuning curves, for a suppression criterion of 1.5 dB. 
% For each measurement condition, e.g., 'allParams.RL_CF_L60', the following were quantified:
%         tipF_re_CF: [10×1 double] % suppression threshold tip frequency (octaves re CF)
%               tipL: [10×1 double] % minimum suppression threshold (dB SPL)
%              q10dB: [10×1 double] % tuning sharpness (Q10dB)
%     int70dBF_re_CF: [10×1 double] % high-frequency 70 dB SPL intercept (octaves re CF)
    
clear;clc;close all;
%% Load data
load('SUPPRESSION_DATA.mat');
supF_re_CF = allData.supF_re_CF; % Suppressor frequency (octaves re CF)

%% Figure panels to plot
figNames = [{'Fig4A'}; {'Fig4B'}; {'Fig4C'}; {'Fig4D'}; {'Fig4E'}; {'Fig4F'}];
figN = length(figNames);

%% Plotting
axlw = 1.8; % axis line width
lw = 3; % figure line width 
selw = 1; % SEM line width
fntSz = 18; % font size
axClr = [0.35 0.35 0.35]; % color for lines indicating CF, probe level, etc.
pltH = 0.3; % plot height

%% Plot each figure panel
for fig_i = 1:figN
    figName = figNames{fig_i};
    panel = figName(end); % panel letter

    switch panel
        case 'A'
            conds = [{'BM_CF_L60'}]; % conditions to plot
            supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
            supCrit_lws = .5 + .3*(1:length(supCrits)); % line width based on suppression criteria
            strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
            plotSEM = 0; % plot SEM? 0 = no; 1 = yes
            minN = 5; % min # mice to include in average data

            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
            plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes

            plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes

        case 'B'
            conds = [{'RL_CF_L60'}]; % conditions to plot
            supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
            supCrit_lws = .5 + .3*(1:length(supCrits)); % line width based on suppression criteria
            strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
            plotSEM = 0; % plot SEM? 0 = no; 1 = yes
            minN = 7; % min # mice to include in average data

            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
            plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes

            plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
            
        case 'C'
            conds = [{'RL_4kHz_L60'}]; % conditions to plot
            supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
            supCrit_lws = .5 + .3*(1:length(supCrits)); % line width based on suppression criteria
            strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
            plotSEM = 0; % plot SEM? 0 = no; 1 = yes
            minN = 7; % min # mice to include in average data

            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
            plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes

            plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes

        case 'D'
            conds = [{'BM_CF_L60'};{'RL_CF_L60'};{'RL_4kHz_L60'}]; % conditions to plot
            cond_lws = [2.5 2.5 2.5];
            supCrits = 1.5; % suppression criterion (dB)
            strictInc = 1; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
            plotSEM = 1; % plot SEM? 0 = no; 1 = yes
            minN = 5; % min # mice to include in average data

            plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
            plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
            plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes

            plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
            
        case 'E'
            conds = [{'BM_CF_L60'};{'RL_CF_L60'};{'RL_4kHz_L60'}]; % conditions to plot
            cond_lws = [2.5 2.5 2.5];
            supCrits = 1.5; % suppression criterion (dB)
            strictInc = 1; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
            plotSEM = 1; % plot SEM? 0 = no; 1 = yes
            minN = 5; % min # mice to include in average data

            plotSupThresh = 0; % plot suppression threshold? 0 = no; 1 = yes
            plotProbeL = 0; % plot dashed line at probe level? 0 = no; 1 = yes
            plotProbeSym = 0; % plot symbol at probe frequency/level? 0 = no; 1 = yes

            plotSupDisp = 1; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
            plotProbeDisp = 1; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
            plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes

        case 'F'
            conds = [{'BM_CF_L60'};{'RL_CF_L60'};{'RL_4kHz_L60'}]; % conditions to plot
            supCrits = 1.5; % suppression criterion (dB)
            plotSEM = 1; % plot SEM? 0 = no; 1 = yes
            strictInc = 1; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
    end

    condN = length(conds); % # conditions
    supCritN = length(supCrits); % # suppression criteria

    if ~strcmp(panel,'F') %% PANELS A-E

        %% If comparing across within the same mice across conditions, get indices of relevant mice
        if condN > 1 && strictInc == 1
            tmp = [];
            for cond_i = 1:condN
                cond = conds{cond_i};
                d = allData.(genvarname(cond));
                tmp = [tmp; d.probeF_re_CFs];
            end 
            tmp = sum(tmp);
            [~,m_ix] = find(~isnan(tmp));
        else
            d = allData.(genvarname(conds{1}));
            m_ix = 1:length(d.probeF_re_CFs); % otherwise include all 
        end


        %% SET UP FIGURE
        h1=figure('units','normalized','position',[0 0 .25 1]);

        %% Suppressor Threshold
        if plotSupThresh
            ax1 = axes('Position',[.12 .66 .8 pltH],'box','off','LineWidth',axlw, 'FontSize', fntSz); hold all;set(gca,'Layer','top');
            xlim([-2.75 1]); set(gca,'Xtick',-3:1:3); ax1.XRuler.MinorTick = 'on'; ax1.XRuler.MinorTickValues = [-3:0.5:3]; ax1.XRuler.TickLength = [0.025 0.025]; %or 'off'
            ylim([30 90]); ax1.YRuler.TickLength = [0.025 0.025];
            plot([0 0], [0 100], ':', 'MarkerSize', 5, 'Linewidth',1.6,'Color',axClr); hold on;
        end

        %% Suppressor-evoked displacement at suppression threshold (nm)
        if plotSupDisp
            ax2 = axes('Position',[.12 .325 .8 pltH],'box','off','LineWidth',axlw,'FontSize',fntSz); hold all;set(gca,'Layer','top');
            xlim([-2.75 1]); set(gca,'Xtick',-3:1:3); ax2.XRuler.MinorTick = 'on'; ax2.XRuler.MinorTickValues = [-3:0.5:3]; ax2.XRuler.TickLength = [0.025 0.025];%or 'off'
            ylim([0.07 70]); set(gca,'Yscale','log','Ytick',[0.1 1 10 50]); ax2.YRuler.TickLength = [0.025 0.025];
            plot([0 0], [0.001 1000], ':','MarkerSize',5,'Linewidth',1.6,'Color',axClr); hold on;
        end

        %% Suppressor-evoked change in probe response phase at suppression threshold (cycles)
        if plotProbePhase  
            ax3 = axes('Position',[.12 .138 .8 .15],'box','off','LineWidth',axlw,'FontSize',fntSz); hold all;set(gca,'Layer','top'); 
            xlim([-2.75 1]); set(gca,'Xtick',-3:1:3); ax3.XRuler.MinorTick = 'on'; ax3.XRuler.MinorTickValues = [-3:0.5:3]; ax3.XRuler.TickLength = [0.025 0.025]; %or 'off'
            ylim([-0.2 0.2]); ax3.YRuler.MinorTick = 'on'; ax3.YRuler.MinorTickValues = [-3:0.1:3]; ax3.YRuler.TickLength = [0.025 0.025];
            plot([0 0], [-5 5], ':','Linewidth',1.6,'Color',axClr); hold on;
            plot([-100 100], [0 0], '--','LineWidth', 1.2, 'Color',axClr); hold on;
        end

        %% Average and plot data
        for cond_i = 1:condN
            cond = conds{cond_i}; % condition
            d = allData.(genvarname(cond)); % data for that condition

            if supCritN == 1
                lw = cond_lws(cond_i);
            end

            %% Parse condition string
            loc = cond(1:2); % location (BM or RL)
            seps = strfind(cond,'_');
            probeF_str = cond(seps(1)+1:seps(2)-1); % probe frequency string
            probeL_str = cond(seps(2)+1:end); % probe level string
            probeL = str2num(probeL_str(2:end)); % probe level (dB SPL)

            %% # Mice w/ data
            mouseN = sum(~isnan(d.probeF_re_CFs(m_ix))); % # mice

            %% Average probe frequency re CF (octaves)
            probeF_re_CF_ave = nanmean(d.probeF_re_CFs(m_ix));
            probeF_re_CF_se = nanstd(d.probeF_re_CFs(m_ix))/sqrt(mouseN); % SEM

            %% Average unsuppressed probe displacement (nm)
            probeDisp_ave = nanmean(d.probeDisp(m_ix)); % average probe displacement (nm)
            probeDisp_se = nanstd(d.probeDisp(m_ix))/sqrt(mouseN); % SEM

            % Specify line colors for each location/probe frequency
            switch loc
                case 'BM'
                    clr = 'k';
                    sym = 'o';
                    lSty = '-'; % line style
                case 'RL'
                    switch probeF_str
                        case 'CF'
                            clr = 'r';
                            sym = '^';
                            lSty = '-'; % line style
                        case '2kHz'
                            clr = [5 31 134]/255;
                            sym = 'v';
                            lSty = '-'; % line style
                        case '3kHz'
                            clr = [121 182 243]/255;
                            sym = 'v';
                            lSty = '-'; % line style
                        case '4kHz'
                            clr = [0 170 200]/255;
                            sym = 'v';
                            lSty = '-'; % line style
                        case '5kHz'
                            clr = [26 186 119]/255;;
                            sym = 'v';
                            lSty = '-'; % line style 
                        case '6kHz'
                            clr = [255 147 36]/255;
                            sym = 'v';
                            lSty = '-'; % line style
                    end
            end

            %% Plot data for each suppression criterion
            for supCrit_i = 1:supCritN
                supCrit = supCrits(supCrit_i);
                [~,supCrit_ix] = ismember(supCrit, allData.supCrits);

                if supCritN > 1
                    lw = supCrit_lws(supCrit_i); % line width if multiple suppression criteria
                end

                %% 1. Suppression Threshold
                if plotSupThresh
                    axes(ax1);
                    supThreshL = squeeze(d.supThreshL(:,supCrit_ix,m_ix));
                    supThreshL_ave = nanmean(supThreshL, 2);
                    supThreshL_sd = nanstd(supThreshL, 0, 2);
                    supThreshL_N = sum(~isnan(supThreshL), 2);
                    supThreshL_se = supThreshL_sd ./(sqrt(supThreshL_N));

                    supThreshL_ave(supThreshL_N < minN) = NaN;

                    % Plot average
                    plot(supF_re_CF, supThreshL_ave,'Color', clr,'LineWidth',lw,'LineStyle',lSty); hold on;
                    if plotSEM % plot SEM?
                        plot(supF_re_CF, supThreshL_ave + supThreshL_se,'Color', clr, 'LineWidth', selw, 'LineStyle','--'); hold on;
                        plot(supF_re_CF, supThreshL_ave - supThreshL_se,'Color', clr, 'LineWidth', selw, 'LineStyle','--'); hold on;
                    end

                    if cond_i == condN && supCrit_i == supCritN
                        if plotProbeL
                            plot([-3 1],[probeL probeL],'--','Color',[.7 .7 .7],'LineWidth',1.5);
                        end
                        if plotProbeSym
                            plot(probeF_re_CF_ave, probeL, 'ko', 'MarkerFaceColor','w','LineWidth',1.5, 'MarkerSize',9);
                        end
                    end
                end

                %% 2. Suppressor displacement
               if plotSupDisp
                    axes(ax2);
                    supDisp = squeeze(d.supDisp_Clean(:,supCrit_ix,m_ix));
                    supDisp_ave = nanmean(supDisp, 2);
                    supDisp_sd = nanstd(supDisp, 0, 2);
                    supDisp_N = sum(~isnan(supDisp), 2);
                    supDisp_se = supDisp_sd ./(sqrt(supDisp_N));

                    supDisp_ave(supDisp_N < minN) = NaN;

                    % Plot average
                    plot(supF_re_CF, supDisp_ave,'Color', clr,'LineWidth',lw,'LineStyle',lSty); hold on;
                    if plotSEM % plot SEM?
                        plot(supF_re_CF, supDisp_ave + supDisp_se,'Color', clr, 'LineWidth', selw, 'LineStyle','--'); hold on;
                        plot(supF_re_CF, supDisp_ave - supDisp_se,'Color', clr, 'LineWidth', selw, 'LineStyle','--'); hold on;
                    end

                    if supCrit_i == supCritN && plotProbeDisp
                        errorbar(probeF_re_CF_ave,probeDisp_ave, probeDisp_se, sym,'Color',clr,'LineWidth',1.5);
                        plot(probeF_re_CF_ave,probeDisp_ave, sym,'Color',clr,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',9);
                    end

                end

                %% 3. Probe phase
                if plotProbePhase
                    axes(ax3);
                    probePhase = squeeze(d.probePhase(:,supCrit_ix,m_ix));
                    probePhase_ave = nanmean(probePhase, 2);
                    probePhase_sd = nanstd(probePhase, 0, 2);
                    probePhase_N = sum(~isnan(probePhase), 2);
                    probePhase_se = probePhase_sd ./(sqrt(probePhase_N));

                    probePhase_ave(probePhase_N < minN) = NaN;

                    % Plot average
                    plot(supF_re_CF, probePhase_ave,'Color', clr,'LineWidth',lw,'LineStyle',lSty); hold on;
                    if plotSEM % plot SEM?
                        plot(supF_re_CF, probePhase_ave + probePhase_se,'Color', clr, 'LineWidth', selw, 'LineStyle','--'); hold on;
                        plot(supF_re_CF, probePhase_ave - probePhase_se,'Color', clr, 'LineWidth', selw, 'LineStyle','--'); hold on;
                    end
                end
            end
        end

    elseif strcmp(panel, 'F') % Fig 4F
        load('SUPPRESSION_PARAMS.mat');

        %% SET UP FIGURE
        h=figure('units','normalized','position',[0 0 .25 .35]);
        ax1 = axes('Position',[.15 .1 .7 .7],'box','off','LineWidth',axlw,'FontSize',fntSz,'XtickLabel',[]); hold all; set(gca,'Layer','top');
        xlim([0.5 3.5]); set(gca,'Xtick',-3:1:10); ax1.XRuler.TickLength = [0.025 0.025];%or 'off'
        ylim([0 1]); set(gca,'Ytick',0:0.2:1); ax1.YRuler.TickLength = [0.025 0.025]; ax1.YRuler.MinorTick = 'on'; ax1.YRuler.MinorTickValues = [-0.1:0.1:1];

        if condN > 1 && strictInc == 1
            tmp = [];
            for cond_i = 1:condN
                cond = conds{cond_i};
                d = allParams.(genvarname(cond));
                tmp = [tmp; d.int70dBF_re_CF'];
            end 
            tmp = sum(tmp);
            [~,m_ix] = find(~isnan(tmp));
        else
            d = allParams.(genvarname(conds{1}));
            m_ix = 1:length(d.int70dBF_re_CF); % otherwise include all 
        end    

        for cond_i = 1:condN
            cond = conds{cond_i}; % condition

            %% Parse condition string
            loc = cond(1:2); % location (BM or RL)
            seps = strfind(cond,'_');
            probeF_str = cond(seps(1)+1:seps(2)-1); % probe frequency string
            probeL_str = cond(seps(2)+1:end); % probe level string
            probeL = str2num(probeL_str(2:end)); % probe level (dB SPL)

            % Specify line colors for each location/probe frequency
            if strcmp(loc,'BM')
                clr = 'k';
            elseif strcmp(loc,'RL') && strcmp(probeF_str,'CF')
                clr = 'r';
            elseif strcmp(loc,'RL') && strcmp(probeF_str,'4kHz')
                clr = [0 170 200]/255;
            end

            int70dB = allParams.(genvarname(cond)).int70dBF_re_CF(m_ix); % 70 dB intercept frequency (re CF)
            n = sum(~isnan(int70dB));
            int70dB_ave = nanmean(int70dB);
            int70dB_sd = nanstd(int70dB);
            int70dB_se = int70dB_sd/sqrt(n);
            tmp = ones(1,n);
            jitter = rand(1,n)/4 + 0.03;
            plot(cond_i+jitter, int70dB, 'o', 'Color', clr,'LineWidth', 0.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w'); hold on;
            errorbar(cond_i, int70dB_ave, int70dB_se, 'o', 'Color', 'k', 'MarkerFaceColor', clr, 'LineWidth', 1.5,'MarkerSize', 12);
            plot(cond_i, int70dB_ave, 'o', 'Color', 'k', 'MarkerFaceColor', clr, 'LineWidth', 1.5,'MarkerSize', 12);
        end
    end
end

