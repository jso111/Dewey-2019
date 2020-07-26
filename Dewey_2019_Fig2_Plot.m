%% Plot Figure 2 from Dewey et al. 2019 J Neurosci paper
% Plots single tone tuning curve data from the BM and RL in live and dead mice (C-H), individual/average RL-BM phase (I-J)
% Loads 'SINGLE_TONE_DATA.mat' file containing the structure 'allData' 
% 'allData' contains single-tone tuning curve data from the BM and RL, both live and dead, in 9 mice

% Fields include:
% allData.f1s; % stimulus frequency (Hz)
% allData.L1s; % stimulus level (dB SPL)

% and the following substructures containing vibration data at the stimulus frequency (f1):
% allData.LIVE.BM.f1; % Live BM data
% allData.LIVE.RL.f1; % Live RL data
% allDAta.DEAD.BM.f1; % Dead BM data
% allData.DEAD.RL.f1; % Dead RL data

% Each substructure contains the following:
%      mag: [29󭘹 double] % displacement magnitude (nm)
%       nf: [29󭘹 double] % mean displacement noise floor (nm)
%     nfsd: [29󭘹 double] % noise floor standard deviation (nm)
%      phi: [29󭘹 double] % displacement phase (radians)
%     magC: [29󭘹 double] % "clean" displacement magnitudes (nm); i.e., magnitudes falling > 3 SD above the mean noise floor
%     phiC: [29󭘹 double] % "clean" displacement phases (radians)

% 'RL' data also include the following fields:
%      phi_re_bm: [29󭘹 double] % RL phase re BM (radians)
%      phi_re_bmC: [29󭘹 double] % "clean" RL phase re BM (radians)

% Each field is organized as: freq x amp x mouse (29x9x9)

clear;clc;close all;
load('SINGLE_TONE_DATA.mat');

%% Plot Fig 2C-H: Representative BM and RL tuning curves from an individual mouse.
mi = 7; % index of representative data to plot
conds = [{'LIVE'}; {'DEAD'}];
condN = length(conds);
locs = [{'BM'};{'RL'}]; % locations = BM, RL
locN = length(locs); % # locations
clrs = [{'k'};{'r'}]; % plot colors
dead_clr = [0.7 0.7 0.7];

f1s = allData.f1s;
f1N = length(f1s);
L1s = allData.L1s;
L1N = length(L1s);

%% Plotting
smoothMag = 1; % 3-point smooth the magnitude data? 0 = no; 1 = yes
lws = 1 + (1:length(L1s))*.3; % line widths
live_L1s = 20:10:90; % levels to plot for live data (dB SPL)
dead_L1s = 60:10:90; % levels to plot for dead data (dB SPL)
axlw = 1.8; % axis line width
fntSz = 18; % font size
pltH=0.273; % plot height

for loc_i = 1:locN
    loc = locs{loc_i}; % location
    live_clr = clrs{loc_i}; % plot color
        
    %% SET UP FIGURE
    h = figure('units','normalized','position',[0 0 .25 1]);
    ax1 = axes('Position',[.12 .72 .8 pltH], 'box','off','FontSize',fntSz,'LineWidth',axlw,'Layer','top'); hold on; % Displacement magnitude (nm)
    ax2 = axes('Position',[.12 .485 .8 0.2],'box','off','FontSize',fntSz,'LineWidth',axlw,'Layer','top'); hold on; % Displacement gain (nm/Pa)
    ax3 = axes('Position',[.12 .04 .8 1.5*pltH],'box','off','FontSize',fntSz,'LineWidth',axlw,'Layer','top'); hold on; % Displacement phase (re EC phi)

    axes(ax1);
    xlim([0.8 14.2]); set(gca,'Xscale','log','Xtick',[1 2 5 10 15]);
    ylim([0.05 70]); set(gca,'Yscale','log', 'Ytick',[0.1 1 10 50]); ax1.YRuler.TickLength = [0.025 0.025];

    axes(ax2);
    xlim([0.8 14.2]); set(gca,'Xscale','log','Xtick',[1 2 5 10 15]);
    ylim([-5 0]); ax2.YRuler.MinorTick = 'on'; ax2.YRuler.MinorTickValues = [-6:1:0]; ax2.YRuler.TickLength = [0.025 0.025];

    axes(ax3);
    xlim([0.8 14.2]); set(gca,'Xscale','log','Xtick',[1 2 5 10 15]);
    ylim([0.2 7000]); set(gca,'Yscale','log','Ytick',[0.1 1 10 100 1000]); ax3.YRuler.TickLength = [0.025 0.025];

    
    %% Plot data for live and dead conditions
    for cond_i = 1:condN
        cond = conds{cond_i};
    
        d = allData.(genvarname(cond)).(genvarname(loc)).f1;
        
        if strcmp(cond,'LIVE')
            plot_L1s = live_L1s;
            clr = live_clr;
            
        elseif strcmp(cond,'DEAD')
            plot_L1s = dead_L1s;
            clr = dead_clr;
        end
        plot_L1N = length(plot_L1s);
        [~,plot_L1is] = ismember(plot_L1s,L1s);

        % Phase references for cleaning
        [~,phi_refL_i] = ismember(90,L1s); % align all curves to w/in 1 cycle of phase curves at 90 dB SPL
        phi_ref = unwrap(squeeze(d.phi(:, phi_refL_i, mi)))/(2*pi); % reference phase for adjusting phase curves
        [~,phi_refF_i] = ismember(8000,f1s); % use phase at 8000 Hz as reference
        phi_ref = phi_ref(phi_refF_i);

        %% Clean and plot data for all stimulus levels
        for L_i = plot_L1N:-1:1  
            L1i = plot_L1is(L_i); % stimulus level index
            L1 = L1s(L1i); % stimulus level (dB SPL)
            lw = lws(L1i); % line width
            LPa = 2e-5 * 10.^(L1/20); % stimulus pressure (Pa)

            mag = squeeze(d.mag(:,L1i,mi)); % displacement (nm)
            phi = squeeze(unwrap(d.phi(:,L1i,mi)))/(2*pi); % displacement phase (cycles)
            nf = squeeze(d.nf(:,L1i,mi)); % mean noise floor (nm)

            % Clean phase
            phi_tol = 0.15;
            phi_tmp = phi;
            for f_i = 2:f1N
                if phi_tmp(f_i) - phi_tmp(f_i-1) > phi_tol
                    phi_tmp = cat(1,phi_tmp(1:f_i-1),phi_tmp(f_i:end)-1); % subtract a whole cycle
                end
            end
            phi = phi_tmp; % replace orig. data
            phi = phi - round((phi(phi_refF_i)-phi_ref),0); % adjust by whole cycles to align w/ reference phase
            
            if smoothMag % 3-point smoothing of mag. and noise floor data
               magS = nan(size(mag));
               nfS = nan(size(nf));
               for f_i = 1:f1N
                    if f_i == 1
                        magS(f_i) = mean(mag(f_i:f_i+1));
                        nfS(f_i) = mean(nf(f_i:f_i+1));
                    elseif f_i == f1N
                        magS(f_i) = mean(mag(f_i-1:f_i));
                        nfS(f_i) = mean(nf(f_i-1:f_i));
                    else
                        magS(f_i) = mean(mag(f_i-1:f_i+1));
                        nfS(f_i) = mean(nf(f_i-1:f_i+1));
                    end
               end
               mag = magS;
               nf = nfS;            
            end
            
            nf = nf + 3 * squeeze(d.nfsd(:,L1i,mi)); % mean noise floor + 3 SDs (nm)
            magC = mag;
            magC(mag < nf) = NaN; % cleaned magnitudes
            phiC = phi;
            phiC(mag < nf) = NaN; % cleaned phase
            
            % interpolate non-continuous data for clarity
            magCx = magC(~isnan(magC));
            phiCx = phiC(~isnan(magC));
            f1sCx = f1s(~isnan(magC));
            for f_i = 2:f1N-1
                if isnan(magC(f_i)) && ~isnan(magC(f_i-1)) && ~isnan(magC(f_i+1)) % if neighbors are clean
                    magC(f_i) = interp1(f1sCx, magCx, f1s(f_i), 'linear');
                    phiC(f_i) = interp1(f1sCx, phiCx, f1s(f_i), 'linear');
                elseif ~isnan(magC(f_i)) && isnan(magC(f_i-1)) && ~isnan(magC(f_i+1)) && isnan(magC(f_i+2)) % eliminate if more than two points missing
                    magC(f_i) = NaN;
                    phiC(f_i) = NaN;
                end
            end
            
            gain = magC./LPa; % displacement (gain re stimulus Pa)

            if strcmp(cond, 'LIVE')
                axes(ax1); % Displacement (nm)
                plot(f1s/1000, magC, 'LineStyle','-','Color',clr,'LineWidth',lw); hold on;
                
                axes(ax2); % Phase (cycles)
                plot(f1s/1000, phiC, 'LineStyle','-','Color',clr,'LineWidth',lw); hold on;
            end
            
            axes(ax3); % Gain re stimulus pressure (nm/Pa)
            plot(f1s/1000, gain,'LineStyle','-','Color',clr,'LineWidth',lw); hold on;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Fig. 2I-J
% (I) Individual BM and RL phases at 70 dB SPL for live mice 
% (J) Individual/average RL - BM phase difference in live and dead mice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mouseN = 9; % # mice
cfs = 1000 * [9; 9.5; 9.5; 9; 9.5; 9.5; 9.5; 9; 9]; % CFs for each mouse
f1is = -3.5:.1:1; % frequency re CF interpolated vector
Lnal = 70; % plot BM and RL phase for 70 dB SPL data
[~,Lnal_i] = ismember(Lnal,L1s);

%% Plotting
axlw = 1.8; % axis line width
lw = 2.8; % figure line width
mrkSz = 5; % marker size
fntSz = 18; % font size

%% Set up figure
h1=figure('units','normalized','position',[0 0 .6 .25]);
ax1 = axes('Position',[.05 .15 .33336 0.7],'box','off','LineWidth',axlw, 'FontSize',fntSz,'Layer','top'); hold on;
ax2 = axes('Position',[.51  .15 .252 0.7],'box','off','LineWidth',axlw, 'FontSize',fntSz,'Layer','top'); hold on;

axes(ax1); % BM and RL phases (cycles)
xlim([-3.5 .6]);  ax1.XRuler.MinorTick = 'on'; ax1.XRuler.MinorTickValues = [-3:0.5:3]; ax1.XRuler.TickLength = [0.025 0.025]; %or 'off'
ylim([-5 0]); ax1.YRuler.MinorTick = 'on'; ax1.YRuler.MinorTickValues = [-6:1:0]; ax1.YRuler.TickLength = [0.025 0.025];

axes(ax2); % RL - BM phase difference (cycles)
xlim([-2.75 .5]); ax2.XRuler.MinorTick = 'on'; ax2.XRuler.MinorTickValues = [-3:0.5:3]; ax2.XRuler.TickLength = [0.025 0.025]; %or 'off'
ylim([-.7 0.1]);  set(gca,'Ytick',[-1:.2:1]); ax2.YRuler.MinorTick = 'on'; ax2.YRuler.MinorTickValues = [-3:0.1:3]; ax2.YRuler.TickLength = [0.025 0.025];
plot([-4 1], [0 0], 'k:','LineWidth',1.5); hold on;  % plot 0 phase difference line

%% Plot BM and RL phase in live mice, RL-BM phase difference for live and dead mice
for cond_i = condN:-1:1
    cond = conds{cond_i};
    
    for loc_i = 1:locN
        loc = locs{loc_i};    
        switch loc
            case 'BM'
                clr = 'k'; % plot color
            case 'RL'
                clr = 'r'; % plot color
        end    
        
        for m_i = 1:mouseN % plot phase re normalized frequency axis for each mouse
            f1_re_cf = log2(f1s/cfs(m_i)); % normalized frequency axis (octaves re CF)
            
            if strcmp(cond,'LIVE')
                phi = allData.LIVE.(genvarname(loc)).f1.phiC(:,Lnal_i,m_i)/(2*pi); % phase in cycles 
                axes(ax1); % Plot phase
                plot(f1_re_cf, phi,'Color', clr, 'LineWidth', 1.2)
            end
            
            if strcmp(loc,'RL')
                if strcmp(cond,'LIVE')
                    diff_clr = [255 203 195]/255;
                else
                    diff_clr = [.75 .75 .75];
                end
                phi_diff = allData.(genvarname(cond)).RL.f1.phi_re_bmC(:,Lnal_i,m_i)/(2*pi); % RL - BM phase (cycles)
                
                % Plot phase difference
                axes(ax2);
                plot(f1_re_cf, phi_diff,'Color', diff_clr, 'LineWidth', 1.2);
                
                % Store interpolated phase difference for averaging
                % purposes
                for L_i = 1:L1N
                    allData.(genvarname(cond)).(genvarname(loc)).f1.phiC_int(:,L_i,m_i) = interp1(f1_re_cf, allData.(genvarname(cond)).(genvarname(loc)).f1.phiC(:,L_i,m_i), f1is, 'linear');
                    allData.(genvarname(cond)).(genvarname(loc)).f1.phi_re_bmC_int(:,L_i,m_i) = interp1(f1_re_cf, allData.(genvarname(cond)).(genvarname(loc)).f1.phi_re_bmC(:,L_i,m_i), f1is, 'linear');
                end
            end
        end
    end
end

%% Plot average RL - BM phase difference for live and dead mice
for cond_i = 1:condN
    cond = conds{cond_i};
    
    if strcmp(cond, 'LIVE')
        clr = 'r';
    else
        clr = [.5 .5 .5];
    end
    
    phi_re_bm = squeeze(allData.(genvarname(cond)).RL.f1.phi_re_bmC_int(:,Lnal_i,:))/(2*pi); % use clean phase difference
    phi_re_bm_ave = nanmean(phi_re_bm,2);
    phi_re_bm_sd = nanstd(phi_re_bm,0,2);
    phi_re_bm_n = sum(~isnan(phi_re_bm),2);
    phi_re_bm_se = phi_re_bm_sd ./ sqrt(phi_re_bm_n);
    
    phi_re_bm_ave(phi_re_bm_n < 5) = NaN; % only show average when clean data are available from >= 5 mice
    
    %% Plot Ave. +/- SE phase difference
    axes(ax2);
    plot(f1is, phi_re_bm_ave, '-','Color',clr, 'LineWidth',2.8);
    plot(f1is, phi_re_bm_ave + phi_re_bm_se, '--', 'Color', clr, 'LineWidth',1);
    plot(f1is, phi_re_bm_ave - phi_re_bm_se, '--', 'Color', clr,'LineWidth',1);
end