%% 
% 
% Script used to plot the ERP data (grand averages, topoplots) in:
% Pomiechowska, B., & Gliga, T.: Nonverbal category knowledge limits the amount of information encoded in object representations
%
% Loop through the data stored in EEGDATA folder, plot and
% save grand average Nc plots and topolots.
%
% Requires EEGLAB.
%
%%

clear; clc; close all;

curr_path = pwd;

%Get list of all EEG files
name = '*.mat';
files    = dir(strcat(curr_path, '/EEGDATA/', name));
fileslist = cell(1,length(files));
for k = 1:length(files)
    fileslist{1,k} = files(k).name;
end

%Loop through EEG files, make topopolots, make ERP plots
for i = 1:length(fileslist)
    
    %Load data.
    curr_dataset = char(fileslist{1,i});
    curr_dataset = curr_dataset(1, 1:end-4);
    indat_file = strcat(curr_path,'/EEGDATA/', curr_dataset);
    load([indat_file '.mat'])
    plotname = curr_dataset(1,1:5);
    
    %Extract file name from full path.
    indat_path_delim = strsplit(indat_file, '\');
    file_hdr = indat_path_delim{end};
    
    %Define sample rate.
    sr = samplingRate;
    
    %Define samping interval (ms).
    sr_ms = 1000/sr;
    
    %Make EEG structure.
    EEG = eeg_emptyset;
    
    %Add data to EEG structure.
    EEG.data(:,:,1) = a1changeNO_Average_Combine; %condition1.
    EEG.data(:,:,2) = a2changeWITHIN_Average_Combine; %condition2.
    EEG.data(:,:,3) = a3changeACROSS_Average_Combine; %condition3.
    
    %Define times.
    times = 0:sr_ms:(length(EEG.data)-1)*sr_ms;
    times=times-200;
    
    %Add times to EEG structure.
    EEG.times = times;
    
    %Redraw EEG structure.
    EEG = eeg_checkset(EEG);
    
    %Remove channel no 125 (reference).
    %EEG = pop_select(EEG, 'nochannel', 125);
    
    %Define channel locations.
    EEG = pop_chanedit(EEG, 'lookup',...
        strcat(curr_path,'/EEGDATA/chans_splines/GSN-HydroCel-129.sfp'),...
        'load',{strcat(curr_path,'/EEGDATA/chans_splines/GSN-HydroCel-129.sfp') 'filetype' 'autodetect'});
    
    %Define time windows boundaries (in ms).
    time_wins = -200:200:1200;
    
    %Get dimensions of topoplot.
    topo_num=length(time_wins)-1;
    topo_dim = [size(EEG.data,3) topo_num];
    
    fig_h = figure; hold on; set(gcf, 'Position', [0         0        1050         1000]);
    
    % Get current axes object (just plotted on) and its position
    ax1 = gca;
    axPos = ax1.Position;
    % Change the position of ax1 to make room for extra axes
    % format is [left bottom width height], so moving up and making shorter here...
    ax1.Position = axPos + [0 0.3 0 -0.3];
    % Exactly the same as for plots (above), axes LineWidth can be changed inline or after
    ax1.LineWidth = 2;
    % Add two more axes objects, with small multiplier for height, and offset for bottom
    ax2 = axes('position', (axPos .* [1 1 1 1e-3]) + [0 0.15 0 0], 'color', 'none', 'linewidth', 2);
    ax3 = axes('position', (axPos .* [1 1 1 1e-3]) + [-.01 0.1 0.02 0], 'color', 'none', 'linewidth', 2);
    % You can change the limits of the new axes using XLim
    ax2.XLim = [0 10];
    ax3.XLim = [-200 1200];
    % You can label the axes using XLabel.String
    ax1.XLabel.String = 'Lambda [nm]';
    ax2.XLabel.String = 'Velocity [m/s]';
    ax3.XLabel.String = 'Time (ms)';
    set(gca, 'fontsize', 50, 'linewidth', 1);
    xticks([-100,100,300,500,700, 900, 1100])
    
    %Loop through conditions.
    for c = 1:size(EEG.data,3)
        
        %Loop through time windows.
        for t = 1:size(time_wins,2)-1
            
            %Extract time window bounds for current condition.
            twin_bounds = time_wins(t:t+1);
            
            %Get indices to EEG samples closest to time window boundaries.
            [~, twin_start_i] = min(abs(EEG.times-twin_bounds(1)));
            [~, twin_stop_i] = min(abs(EEG.times-twin_bounds(2)));
            
            %Extract data for time window.
            topodat_win = EEG.data(:,twin_start_i:twin_stop_i,c); %- EEG.data(:,twin_start_i:twin_stop_i,2);
            
            %Now get mean across window (average across time).
            topodat_winmean = squeeze(mean(topodat_win,2));
            
            %Now plot topographies of window mean for each condition.
            fig_h; subplot(topo_dim(1)+3, topo_dim(2), t+(topo_dim(2)*(c-1))); hold on; topoplot(topodat_winmean, EEG.chanlocs, ...
                'electrodes', 'off', 'style', 'map', 'maplimits', [-22,22],'headrad', (.5), 'plotrad', (.55));  %maplimits to change scaling.
            
            if c == 3
                %title([num2str(twin_bounds(1)+100)], 'FontSize', 36,'FontWeight','Normal'); % ' - ' num2str(twin_bounds(2)) 'ms']
            end
            
        end  %End time window loop.
        
    end %End condition loop.
    
    
    %ce=cbar;
    
    %Loop
    for a = 4:5
        %Get the difference data.
        for t = 1:size(time_wins,2)-1
            
            %Extract time window bounds for current condition.
            twin_bounds = time_wins(t:t+1);
            
            %Get indices to EEG samples closest to time window boundaries.
            [~, twin_start_i] = min(abs(EEG.times-twin_bounds(1)));
            [~, twin_stop_i] = min(abs(EEG.times-twin_bounds(2)));
            
            %Extract data for time window.
            topodat_win1 = EEG.data(:,twin_start_i:twin_stop_i,a-2) - EEG.data(:,twin_start_i:twin_stop_i,1);
            
            %Now get mean across window (average across time).
            topodat_winmean = squeeze(mean(topodat_win1,2));
            
            %Now plot topographies of window mean for each condition.
            
            fig_h; subplot(topo_dim(1)+3, topo_dim(2), t+(topo_dim(2)*(a-1))); hold on; topoplot(topodat_winmean, EEG.chanlocs, ...
                'electrodes', 'off', 'style', 'map','maplimits',[-12,12], 'headrad', (.5), 'plotrad', (.55)); %maplimits to change scaling.
            
        end  %End time window loop.
        
    end % End
        
    set(fig_h, 'color', 'w');
    
    saveas(fig_h,['Topo' plotname '.png']);
    
    
    %% ERP PLOT
    
    %Now make ERP plot.
    
    fig_erp = figure; hold on; set(gcf, 'Position', [0 0 800 1000], 'color', 'w');
    
    EEG.times = times(1:850);
    
    r = rectangle('Position',[276 -19.95 522 30]);
    r.FaceColor = [0.95 0.95 0.95]; r.EdgeColor = [0.95 0.95 0.95];
    
    roi_dat = squeeze(mean(EEG.data([3, 4, 5, 6, 7, 10, 11, 12, 13, 16, 18, 19, 20, 23, 24, 28, 29, 30, 35, 36, 104, 105, 106, 110, 111, 112, 117, 118, 124],1:850,:))); %roi plot.
    
    %% Smoothing
    roi_dat(:,1) = smoothdata(roi_dat(:,1),'loess',20);
    roi_dat(:,2) = smoothdata(roi_dat(:,2),'loess',20);
    roi_dat(:,3) = smoothdata(roi_dat(:,3),'loess',20);
    
    plot(EEG.times, roi_dat(:,1), 'b', EEG.times, roi_dat(:,3),'black', EEG.times, roi_dat(:,2), 'r', 'linewidth', 6);
    set(gca, 'Layer', 'top')
    line([-200,1400],[0,0],'Color',[.3 .3 .3])
    line([0,0],[-20,15],'Color',[.3 .3 .3])
    ylim([-20,10])
    xlim([-200,1400])
    xticks([0 500 1000])
    
    xlabel('Time (ms)');
    ylabel('Amplitude (\muv)');
    set(gca, 'fontsize', 50, 'linewidth', 1); % 'fontweight', 'bold',
    
    set(fig_erp, 'color', 'w');
    
    saveas(fig_erp,['ERP_' plotname '.png']);
    
end
