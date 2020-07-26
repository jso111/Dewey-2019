%% Plot Figure 6 from Dewey et al. 2019 J Neurosci paper
% Loads 'SUPPRESSION_DATA.mat' and 'SUPPRESSION_PARAMS.mat'
% See 'Dewey_2019_Fig4_Plot.m' for explanation of data formatting.

clear;clc;close all;

%% Load data
load('SUPPRESSION_DATA.mat');
load('SUPPRESSION_PARAMS.mat');
supF_re_CF = allData.supF_re_CF; % Suppressor frequency (octaves re CF)

%% Figure panel to plot
figNames = [{'Fig6A'}; {'Fig6B'}; {'Fig6C'}; {'Fig6D'};{'Fig6E'}; {'Fig6F'}]; 
figNameN = length(figNames);

%% Plotting
lw_dflt = 2.5; % default line width
axlw = 1.8; % axis line width
selw = 1; % SE line width
fntSz = 18; % font size
axClr = [0.35 0.35 0.35]; % color for lines indicating CF, probe level, etc.
pltH=0.3; % plot height

%% Plot each figure panel
for fig_i = 1:figNameN
    figName = figNames{fig_i}; 
    panel = figName(end); % panel
    
    if ~strcmp(panel,'F')  
        switch panel
            case 'A'
                conds = [{'RL_2kHz_L60'}]; % conditions to plot
                supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
                supCrit_lws = .5 + .3*(1:length(supCrits)); % line widths based on suppression criteria
                plotSEM = 0; % plot SEM? 0 = no; 1 = yes
                strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
                plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
                plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes
                plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
                plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
                plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
                plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
                lSty = '-'; % line style
                
            case 'B'
                conds = [{'RL_3kHz_L60'}]; % conditions to plot
                supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
                supCrit_lws = .5 + .3*(1:length(supCrits)); % line widths based on suppression criteria
                plotSEM = 0; % plot SEM? 0 = no; 1 = yes
                strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
                plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
                plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes
                plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
                plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
                plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
                plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
                lSty = '-'; % line style
            case 'C'
                conds = [{'RL_5kHz_L60'}]; % conditions to plot
                supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
                supCrit_lws = .5 + .3*(1:length(supCrits)); % line widths based on suppression criteria
                plotSEM = 0; % plot SEM? 0 = no; 1 = yes
                strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
                plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
                plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes
                plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
                plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
                plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
                plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
                lSty = '-'; % line style
            case 'D'
                conds = [{'RL_6kHz_L60'}]; % conditions to plot
                supCrits = 1.5:1.5:12; % suppression criteria to plot (dB)
                supCrit_lws = .5 + .3*(1:length(supCrits)); % line widths based on suppression criteria
                plotSEM = 0; % plot SEM? 0 = no; 1 = yes
                strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
                plotProbeL = 1; % plot dashed line at probe level? 0 = no; 1 = yes
                plotProbeSym = 1; % plot symbol at probe frequency/level? 0 = no; 1 = yes
                plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
                plotSupDisp = 0;% plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
                plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
                plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
                lSty = '-'; % line style
            case 'E'
                conds = [{'RL_2kHz_L60'};{'RL_4kHz_L60'};{'RL_6kHz_L60'};{'RL_CF_L60'}]; % conditions to plot
                cond_lws = [2.5 2.5 2.5 2.5];
                supCrits = 1.5; % suppression criteria to plot (dB)
                plotSEM = 1; % plot SEM? 0 = no; 1 = yes
                strictInc = 0; % only include mice w/ data for all conditions? 0 = no (include data from all mice); 1 = yes
                plotProbeL = 0; % plot dashed line at probe level? 0 = no; 1 = yes
                plotProbeSym = 0; % plot symbol at probe frequency/level? 0 = no; 1 = yes
                plotSupThresh = 1; % plot suppression threshold? 0 = no; 1 = yes
                plotSupDisp = 0; % plot suppressor-evoked displacement at suppression threshold? 0 = no; 1 = yes
                plotProbeDisp = 0; % plot symbol for ave. unsuppressed probe-evoked displacement? 0 = no; 1 = yes
                plotProbePhase = 0; % plot change in probe response phase at suppression threshold? 0 = no; 1 = yes
                lSty = '-'; % line style
        end

        condN = length(conds); % # conditions
        supCritN = length(supCrits); % # suppression criteria

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
        h = figure('units','normalized','position',[0 0 .25 1]);

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
                if condN > 1
                    lw = cond_lws(cond_i); % if only one suppression criterion, set line widths based on condition
                else
                    lw = lw_dflt; % use default line width
                end
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

            % Specify line colors for each locure/probe frequency
            switch probeF_str

                case 'CF'
                   clr = 'r';
                   sym = 'o';
                   minN = 9;
                case '2kHz'
                   clr = [5 31 134]/255;
                   sym = 'o';
                   minN = 4;
                case '3kHz'
                   clr = [121 182 243]/255;
                   sym = 'o';
                   minN = 5;
                case '4kHz'
                   clr = [0 170 200]/255;
                   sym = 'o';
                   minN = 8;
                case '5kHz'
                   clr = [26 186 119]/255;
                   sym = 'o';
                   minN = 4;
                case '6kHz'
                   clr = [255 147 36]/255;
                   sym = 'o';
                   minN = 3;
            end


            %% Plot data for each suppression criterion
            for supCrit_i = 1:supCritN
                supCrit = supCrits(supCrit_i);
                [~,supCrit_ix] = ismember(supCrit, allData.supCrits);

                if supCritN > 1
                    lw = supCrit_lws(supCrit_i); % line width if multiple suppression criteria
                end

                %% 1. Suppressor Level
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
                        errorbar(probeF_re_CF_ave, probeDisp_ave, probeDisp_se, sym,'Color',clr,'LineWidth',1.5);
                        plot(probeF_re_CF_ave, probeDisp_ave, sym,'Color',clr,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',9);
                    end

                end

                %% 3. Change in probe phase at suppression threshold
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
        
    elseif strcmp(panel,'F')
        conds = [{'RL_2kHz_L60'};{'RL_3kHz_L60'};{'RL_4kHz_L60'};{'RL_5kHz_L60'};{'RL_6kHz_L60'};{'RL_CF_L60'}]; % conditions to plot
        condN = length(conds);

        %% SET UP FIGURE
        h = figure('units','normalized','position',[0 0 .25 1]);
        ax1 = axes('Position',[.12 .66 .8 pltH],'box','off','LineWidth',axlw,'FontSize',fntSz); hold all;set(gca,'Layer','top');
        xlim([-2.5 .2]); set(gca,'Xtick',-3:1:10); ax1.XRuler.MinorTick = 'on'; ax1.XRuler.MinorTickValues = [-3:0.5:3]; ax1.XRuler.TickLength = [0.025 0.025];%or 'off'
        ylim([0 0.8]); set(gca,'Ytick',0:0.2:1); ax1.YRuler.MinorTick = 'on'; ax1.YRuler.MinorTickValues = [-0.1:0.1:1]; ax1.YRuler.TickLength = [0.025 0.025];

        for cond_i = 1:condN
            cond = conds{cond_i}; % condition
            %% Parse condition string
            loc = cond(1:2); % location (BM or RL)
            seps = strfind(cond,'_');
            probeF_str = cond(seps(1)+1:seps(2)-1); % probe frequency string
            probeL_str = cond(seps(2)+1:end); % probe level string
            probeL = str2num(probeL_str(2:end)); % probe level (dB SPL)

            % Specify colors for each probe frequency
            switch probeF_str
                case 'CF'
                   clr = 'r';
                case '2kHz'
                   clr = [5 31 134]/255;
                case '3kHz'
                    clr = [121 182 243]/255;
                case '4kHz'
                    clr = [0 170 200]/255;
                case '5kHz'
                    clr = [26 186 119]/255;
                case '6kHz'
                    clr = [255 147 36]/255;
            end

            int70dBF_re_CF = allParams.(genvarname(cond)).int70dBF_re_CF(:,supCrit_i); % 70 dB intercepts (octaves re CF)
            probeF_re_CF = allData.(genvarname(cond)).probeF_re_CFs; % probe frequency (octaves re CF)
            probeF_re_CF_ave = nanmean(probeF_re_CF);
            
            n = sum(~isnan(int70dBF_re_CF)); % # mice w/ data
            int70dBF_re_CF_ave = nanmean(int70dBF_re_CF); % ave.
            int70dBF_re_CF_sd = nanstd(int70dBF_re_CF); % SD
            int70dBF_re_CF_se = int70dBF_re_CF_sd/sqrt(n); % SE

            jitter = rand(1,n)/4;
            plot(probeF_re_CF, int70dBF_re_CF, 'o', 'Color', clr,'LineWidth', 0.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w'); hold on;
            errorbar(probeF_re_CF_ave, int70dBF_re_CF_ave, int70dBF_re_CF_se, 'o', 'Color', 'k', 'MarkerFaceColor', clr, 'LineWidth', 1.5,'MarkerSize', 10);
        end
    end
end
