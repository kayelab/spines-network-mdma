%% ── Load data ─────────────────────────────────────────────────────────────
load([pwd '\data\cue_traces.mat']);
load([pwd '\data\baselined_freezing.mat']);
all_traces = [cue_traces.Saline; cue_traces.MDMA];   % Saline rows 1-8, MDMA rows 9-16
all_freeze = [pre_cue_corrected_freezing.Saline; pre_cue_corrected_freezing.MDMA];

%% ── Color palette ────────────────────────────────────────────────────────
MDMA_color   = [173,  84,  50] / 255;
Saline_color = [ 88,  88,  90] / 255;
fontSz       = 14;

cluster_cmap = [
    0.4510    0.3765    0.5216;   % cluster 1 – purple-gray
    0.6196    0.7922    0.8392;   % cluster 2 – blue
    0.9569    0.7020    0.4431;   % cluster 3 – light orange
    0.4353    0.5608    0.4471;   % cluster 4 – dark green
];

%% ── Parameters ───────────────────────────────────────────────────────────
treatments        = {'Saline'; 'MDMA'};
treatment_idx     = struct('Saline', 1:8, 'MDMA', 9:16);
n_mice_per_group  = 8;
total_number_of_mice       = 16;

frame_rate        = 30;    % Hz
time_before_cue   = 30;    % s
cue_time          = 30;    % s

days                   = 2:5;
n_days_ext             = numel(days);
n_days_total           = 5;
adj_day_pairs          = [2 3; 3 4; 4 5]; 
n_adj                  = size(adj_day_pairs,1);
% Tuning-curve analysis parameters (panels B–D)
first_day_for_analysis = 2;
n_sec_onset            = 7;
cue_onset_frame        = frame_rate * time_before_cue;
cue_offset_frame       = cue_onset_frame + frame_rate * cue_time;
shock_start_frame      = cue_onset_frame + frame_rate * (cue_time - 1) - 1;
onset_frames           = cue_onset_frame : cue_onset_frame + frame_rate * n_sec_onset - 1;
n_onset_frames         = length(onset_frames);
cues_per_day           = [1, 6, 8, 10, 10];
n_cues                 = sum(cues_per_day);
cue_day_ID_mat(:,1)    = [1 1 2 3 4 5 6 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10];
cue_day_ID_mat(:,2)    = [1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5];

% Panel H parameters — window within each cue for cluster correlation analysis
% h_win_sec    = [30 60];   % seconds from trial start
% nonempty_traces = all_traces(~cellfun(@isempty, all_traces));
% n_t          = size(nonempty_traces{1}, 2);   % time points per trial (robust to empty day-1 cells)
% h_win_idx    = [max(1, round(h_win_sec(1)*frame_rate)), min(n_t, round(h_win_sec(2)*frame_rate))];
% h_win        = h_win_idx(1):h_win_idx(2);
% h_smooth_win = max(1, round(1 * frame_rate));   % 1-s smoothing window
% h_bin        = 1;                               % no downsampling (1 sample/bin)

%% ── Analysis: tuning curve correlations (panels B–D) ─────────────────────
clear per_mouse_mean_corr_per_time_diff shuff_per_mouse_mean_corr_per_time_diff ...
      within_day_corr_per_mouse within_day_shuff_corr_per_mouse

for run_treatments = 1:2
    curr_treatment = treatments{run_treatments};

    for run_mice = 1:n_mice_per_group
        disp([curr_treatment ' mouse ' num2str(run_mice)])

        curr_mouse_data = cue_traces.(curr_treatment)(run_mice, :);
        n_cells         = size(curr_mouse_data{1}, 1);

        % Build per-mouse activity matrices (cells × frames × cues)
        onset_mat = [];
        for run_days = 1:n_days_total
            curr_day_act = curr_mouse_data{run_days};
            if isempty(curr_day_act)
                day_onset_mat = nan(n_cells, n_onset_frames, cues_per_day(run_days));
            else
                for run_cues = 1:cues_per_day(run_days)
                    curr_cue_mat = curr_day_act(:, :, run_cues);
                    curr_day_act(:, :, run_cues) = (zscore(curr_cue_mat'))';
                end
                day_onset_mat = curr_day_act(:, onset_frames, :);
            end
            onset_mat = cat(3, onset_mat, day_onset_mat);
        end

        % Replace cue 1 with post-shock activity from day 1
        day1_act       = curr_mouse_data{1};
        post_shock_mat = day1_act(:, shock_start_frame : (shock_start_frame + n_onset_frames - 1));
        onset_mat(:, :, 1) = post_shock_mat;

        % Temporal downsample to 0.5 Hz
        onset_target_frames = n_onset_frames / (frame_rate / 2);
        win_onset           = floor(n_onset_frames / onset_target_frames);
        sm_onset            = movmean(onset_mat, win_onset, 2, 'omitnan', 'Endpoints', 'shrink');
        centers_onset       = round(linspace(ceil(win_onset/2), n_onset_frames - floor(win_onset/2), onset_target_frames));
        ds_onset_mat        = sm_onset(:, centers_onset, :);   % cells × ds_frames × cues

        % Build cue × cue correlation matrix averaged over cells
        cue_corr_by_cell   = nan(n_cues, n_cues, n_cells);
        shuff_corr_by_cell = nan(n_cues, n_cues, n_cells);

        for run_cells = 1:n_cells
            cell_data     = squeeze(ds_onset_mat(run_cells, :, :));   % ds_frames × cues
            cell_corr_mat = corr(cell_data, 'Rows', 'pairwise');
            cell_corr_mat(logical(eye(size(cell_corr_mat)))) = NaN;
            cue_corr_by_cell(:, :, run_cells) = cell_corr_mat;

            for run_cues1 = 1:n_cues-1
                cue1_act = cell_data(:, run_cues1);
                if sum(isnan(cue1_act)) < 1
                    for run_cues2 = run_cues1:n_cues
                        cue2_mat     = ds_onset_mat(:, :, run_cues2);
                        active_cells = find(~isnan(sum(cue2_mat, 2)));
                        if ~isempty(active_cells)
                            rand_cell_act = cue2_mat(active_cells(randi(numel(active_cells))), :);
                            shuff_corr    = corr(cue1_act(:), rand_cell_act(:));
                            shuff_corr_by_cell(run_cues1, run_cues2, run_cells) = shuff_corr;
                            shuff_corr_by_cell(run_cues2, run_cues1, run_cells) = shuff_corr;
                        end
                    end
                end
            end
        end

        curr_mouse_corr_mat = mean(cue_corr_by_cell, 3, 'omitnan');
        curr_mouse_corr_mat(logical(eye(size(curr_mouse_corr_mat)))) = NaN;
        shuff_corr_mat = mean(shuff_corr_by_cell, 3, 'omitnan');
        all_mice_shuff_corr_mats.(curr_treatment)(:, :, run_mice) = shuff_corr_mat;

        % Fill lower triangle with day-pair means
        day_vec  = cue_day_ID_mat(:, 2);
        day_mean = nan(n_days_total, n_days_total);
        for run_day1 = 1:n_days_total
            rows_day1 = find(day_vec == run_day1);
            for run_day2 = 1:n_days_total
                cols_day2 = find(day_vec == run_day2);
                temp_sub  = curr_mouse_corr_mat(rows_day1, cols_day2);
                if run_day1 == run_day2
                    kk = min(size(temp_sub,1), size(temp_sub,2));
                    temp_sub(1:kk+1:kk*kk) = NaN;
                end
                day_mean(run_day1, run_day2) = mean(temp_sub(:), 'omitnan');
            end
        end

        corr_mat_filled = curr_mouse_corr_mat;
        for run_cue_row = 1:n_cues
            for run_cue_col = 1:(run_cue_row-1)
                day_of_row = day_vec(run_cue_row);
                day_of_col = day_vec(run_cue_col);
                val = day_mean(day_of_row, day_of_col);
                if isnan(val), val = day_mean(day_of_col, day_of_row); end
                corr_mat_filled(run_cue_row, run_cue_col) = val;
            end
        end
        all_mice_corr_mats.(curr_treatment)(:, :, run_mice) = corr_mat_filled;

        % Bin correlations by day distance and within-day
        if first_day_for_analysis == 2
            first_cue_idx = 2;
        else
            first_cue_idx = 8;
        end

        time_diff_corr_arr        = cell(1, 5);
        shuff_time_diff_corr_arr  = cell(1, 5);
        within_day_corr_arr       = cell(1, 5);
        shuff_within_day_corr_arr = cell(1, 5);

        for run_cues1 = first_cue_idx:n_cues-1
            cue1_day = cue_day_ID_mat(run_cues1, 2);
            for run_cues2 = run_cues1:n_cues
                cue2_day = cue_day_ID_mat(run_cues2, 2);
                day_dist = cue2_day - cue1_day;
                time_diff_corr_arr{day_dist+1}       = [time_diff_corr_arr{day_dist+1}       curr_mouse_corr_mat(run_cues1, run_cues2)];
                shuff_time_diff_corr_arr{day_dist+1} = [shuff_time_diff_corr_arr{day_dist+1} shuff_corr_mat(run_cues1, run_cues2)];
                if cue1_day == cue2_day
                    within_day_corr_arr{cue1_day}       = [within_day_corr_arr{cue1_day}       curr_mouse_corr_mat(run_cues1, run_cues2)];
                    shuff_within_day_corr_arr{cue1_day} = [shuff_within_day_corr_arr{cue1_day} shuff_corr_mat(run_cues1, run_cues2)];
                end
            end
        end

        for run_diff = 1:4
            per_mouse_mean_corr_per_time_diff.(curr_treatment)(run_mice, run_diff)       = mean(time_diff_corr_arr{run_diff},       'omitnan');
            shuff_per_mouse_mean_corr_per_time_diff.(curr_treatment)(run_mice, run_diff) = mean(shuff_time_diff_corr_arr{run_diff}, 'omitnan');
        end
        for run_days = 2:5
            within_day_corr_per_mouse.(curr_treatment)(run_mice, run_days-1)       = mean(within_day_corr_arr{run_days},       'omitnan');
            within_day_shuff_corr_per_mouse.(curr_treatment)(run_mice, run_days-1) = mean(shuff_within_day_corr_arr{run_days}, 'omitnan');
        end
    end
end

%% ── Panel B: cue-evoked correlation matrices ─────────────────────────────
mean_corr_mat.Saline = mean(all_mice_corr_mats.Saline(2:end, 2:end, :), 3, 'omitnan');
mean_corr_mat.MDMA   = mean(all_mice_corr_mats.MDMA(2:end,   2:end, :), 3, 'omitnan');

min_val = min([mean_corr_mat.Saline(:); mean_corr_mat.MDMA(:)]);
max_val = max([mean_corr_mat.Saline(:); mean_corr_mat.MDMA(:)]);

figure;
set(gcf, 'Position', get(0, 'Screensize'), 'Renderer', 'painters');
[ax1, ax2, cb_pos] = heatmap_axes_layout();
imagesc(ax1, mean_corr_mat.Saline);
title(ax1, 'Saline mice', 'FontSize', fontSz, 'FontWeight', 'normal', 'Color', Saline_color);
imagesc(ax2, mean_corr_mat.MDMA);
title(ax2, 'MDMA mice',   'FontSize', fontSz, 'FontWeight', 'normal', 'Color', MDMA_color);
cb = colorbar(ax2);
for ax = [ax1, ax2]
    axis(ax, 'square'); box(ax, 'on');
    ax.XTick = [3.5 10.5 19.5 29.5]; ax.YTick = ax.XTick;
    ax.XTickLabel = num2cell(2:5)';
    xlabel(ax, 'Day', 'FontSize', fontSz);
    if ax == ax1
        ylabel(ax, 'Day', 'FontSize', fontSz);
        ax.YTickLabel = ax.XTickLabel;
    else
        ax.YTickLabel = [];
    end
    clim(ax, [min_val max_val]); colormap(ax, jet);
    hold(ax, 'on');
    for run_days = 2:4
        sep = find(cue_day_ID_mat(:,2) == run_days, 1, 'last') - 1 + 0.5;
        plot(ax, [sep sep],        [0.5 n_cues+0.5], '-w', 'LineWidth', 2);
        plot(ax, [0.5 n_cues+0.5], [sep sep],        '-w', 'LineWidth', 2);
    end
    style_axes(ax, fontSz);
end
set(cb, 'Position', cb_pos); cb.FontSize = fontSz;
cb.Label.String = 'Tuning curve correlation'; cb.Label.FontSize = fontSz;
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = 'bottom';

%% ── Panel C: correlation vs day difference + within-day ──────────────────
figure;
set(gcf, 'Position', get(0, 'Screensize'), 'Renderer', 'painters');

subplot(1, 2, 1); hold on;
errorbar(1:4, mean(per_mouse_mean_corr_per_time_diff.Saline, 'omitnan'), std(per_mouse_mean_corr_per_time_diff.Saline, 'omitnan')/sqrt(n_mice_per_group), '-*',  'Color', Saline_color, 'LineWidth', 2);
errorbar(1:4, mean(per_mouse_mean_corr_per_time_diff.MDMA,   'omitnan'), std(per_mouse_mean_corr_per_time_diff.MDMA,   'omitnan')/sqrt(n_mice_per_group), '-*',  'Color', MDMA_color,   'LineWidth', 2);
errorbar(1:4, mean(shuff_per_mouse_mean_corr_per_time_diff.Saline, 'omitnan'), std(shuff_per_mouse_mean_corr_per_time_diff.Saline, 'omitnan')/sqrt(n_mice_per_group), '--*', 'Color', Saline_color, 'LineWidth', 2);
errorbar(1:4, mean(shuff_per_mouse_mean_corr_per_time_diff.MDMA,   'omitnan'), std(shuff_per_mouse_mean_corr_per_time_diff.MDMA,   'omitnan')/sqrt(n_mice_per_group), '--*', 'Color', MDMA_color,   'LineWidth', 2);
ylabel('Cue evoked ensemble activity correlation', 'FontSize', fontSz);
xlabel('Days difference', 'FontSize', fontSz);
legend({'Saline','MDMA','Saline - Shuffle','MDMA - Shuffle'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', fontSz);
xlim([0.95 4.05]); ylim([-0.005 0.101]);
ax = gca; ax.XTick = 1:4; ax.XTickLabel = {'0','1','2','3'};
axis square; box off;
style_axes(ax, fontSz);

subplot(1, 2, 2); hold on;
errorbar(1:4, mean(within_day_corr_per_mouse.Saline,        'omitnan'), std(within_day_corr_per_mouse.Saline,        'omitnan')/sqrt(n_mice_per_group), '-*',  'Color', Saline_color, 'LineWidth', 2);
errorbar(1:4, mean(within_day_corr_per_mouse.MDMA,          'omitnan'), std(within_day_corr_per_mouse.MDMA,          'omitnan')/sqrt(n_mice_per_group), '-*',  'Color', MDMA_color,   'LineWidth', 2);
errorbar(1:4, mean(within_day_shuff_corr_per_mouse.Saline,  'omitnan'), std(within_day_shuff_corr_per_mouse.Saline,  'omitnan')/sqrt(n_mice_per_group), '--*', 'Color', Saline_color, 'LineWidth', 2);
errorbar(1:4, mean(within_day_shuff_corr_per_mouse.MDMA,    'omitnan'), std(within_day_shuff_corr_per_mouse.MDMA,    'omitnan')/sqrt(n_mice_per_group), '--*', 'Color', MDMA_color,   'LineWidth', 2);
ylabel('Within-day cue evoked ensemble activity correlation', 'FontSize', fontSz);
xlabel('Day', 'FontSize', fontSz);
legend({'Saline','MDMA','Saline - Shuffle','MDMA - Shuffle'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', fontSz);
xlim([0.95 4.05]);
ax = gca; ax.XTick = 1:4; ax.XTickLabel = num2cell(2:5)';
axis square; box off;
style_axes(ax, fontSz);

%% ── Panel D: shuffle correlation matrices ────────────────────────────────
mean_shuff_mat.Saline = mean(all_mice_shuff_corr_mats.Saline(2:end, 2:end, :), 3, 'omitnan');
mean_shuff_mat.MDMA   = mean(all_mice_shuff_corr_mats.MDMA(2:end,   2:end, :), 3, 'omitnan');
min_val = min([mean_shuff_mat.Saline(:); mean_shuff_mat.MDMA(:)]);
max_val = max([mean_shuff_mat.Saline(:); mean_shuff_mat.MDMA(:)]);

figure;
set(gcf, 'Position', get(0, 'Screensize'), 'Renderer', 'painters');
[ax1, ax2, cb_pos] = heatmap_axes_layout();
imagesc(ax1, mean_shuff_mat.Saline);  title(ax1, 'Saline - Shuffle', 'FontWeight', 'bold', 'FontSize', fontSz);
imagesc(ax2, mean_shuff_mat.MDMA);    title(ax2, 'MDMA - Shuffle',   'FontWeight', 'bold', 'FontSize', fontSz);
cb = colorbar(ax2);
for ax = [ax1, ax2]
    axis(ax, 'square'); box(ax, 'on');
    ax.XTick = [3.5 10.5 19.5 29.5]; ax.YTick = ax.XTick;
    ax.XTickLabel = num2cell(2:5)';   ax.YTickLabel = ax.XTickLabel;
    xlabel(ax, 'Day', 'FontSize', fontSz); ylabel(ax, 'Day', 'FontSize', fontSz);
    clim(ax, [min_val max_val]); colormap(ax, jet);
    hold(ax, 'on');
    for run_days = 2:4
        sep = find(cue_day_ID_mat(:,2) == run_days, 1, 'last') - 1 + 0.5;
        plot(ax, [sep sep],        [0.5 n_cues+0.5], '-w', 'LineWidth', 2);
        plot(ax, [0.5 n_cues+0.5], [sep sep],        '-w', 'LineWidth', 2);
    end
    style_axes(ax, fontSz);
end
set(cb, 'Position', cb_pos); cb.FontSize = fontSz;
cb.Label.String = 'Cue evoked correlation'; cb.Label.FontSize = fontSz;
cb.Label.Rotation = 270; cb.Label.VerticalAlignment = 'bottom';

%% ── Analysis: hierarchical clustering ────────────────────────────────────
[cluster_info] = hierarchical_clustering(all_traces, frame_rate);
t_clust  = (0:2700) / frame_rate;
k_clust  = max(cluster_info.clust_id0);

%% ── Figure E: dendrogram, colorstrip, and heatmap ────────────────────────
figure;
set(gcf, 'Position', get(0, 'Screensize'));

ax_dend = subplot(1, 3, 1);
[H, ~, leaf_order_display] = dendrogram(cluster_info.Z0, 0, 'Orientation', 'right');
set(ax_dend, 'YDir', 'reverse');
set(H, 'LineWidth', 1.5);
box on; grid off; yticks([]);
for run_seg = 1:numel(H)   % warp x-axis for visual clarity
    xd = get(H(run_seg), 'XData');
    set(H(run_seg), 'XData', xd.^4);
end
Xplot      = cluster_info.X0(leaf_order_display, :);
clust_plot = cluster_info.clust_id0(leaf_order_display);
nLeaves    = numel(leaf_order_display);
for run_seg = 1:numel(H)   % color dendrogram segments by cluster
    yd  = get(H(run_seg), 'YData');
    ylo = max(1, floor(min(yd)));
    yhi = min(nLeaves, ceil(max(yd)));
    if yhi < ylo, continue; end
    seg_clusters = unique(clust_plot(ylo:yhi));
    if numel(seg_clusters) == 1
        set(H(run_seg), 'Color', cluster_cmap(seg_clusters(1), :));
    else
        set(H(run_seg), 'Color', 'k');
    end
end
style_axes(ax_dend, fontSz);

ax_strip = subplot(1, 3, 2);
imagesc(clust_plot);
colormap(ax_strip, cluster_cmap);
set(ax_strip, 'YDir', 'reverse');
box off; xticks([]); yticks([]);
style_axes(ax_strip, fontSz);

ax_heat = subplot(1, 3, 3);
imagesc(t_clust, 1:size(Xplot, 1), Xplot);
colormap(ax_heat, flipud(gray));
set(ax_heat, 'YDir', 'reverse');
axis tight; clim([-2 2]); colorbar;
xlabel('Time (s)', 'FontSize', fontSz);
ylabel('Neuron',   'FontSize', fontSz);
style_axes(ax_heat, fontSz);

%% ── Figure E: cluster mean traces ────────────────────────────────────────
figure;
set(gcf, 'Position', get(0, 'Screensize'));
for run_clusters = 1:k_clust
    subplot(1, k_clust, run_clusters); hold on;
    Xc  = cluster_info.X0(cluster_info.clust_id0 == run_clusters, :);
    m   = mean(Xc, 1, 'omitnan');
    sem = std(Xc, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Xc), 1));
    fill([t_clust fliplr(t_clust)], [m-sem fliplr(m+sem)], cluster_cmap(run_clusters,:), ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t_clust, m, 'Color', cluster_cmap(run_clusters,:), 'LineWidth', 2);
    xline([30 60], 'k--', 'LineWidth', 1);
    xlim([t_clust(1) t_clust(end)]); ylim([-2 2]);
    xlabel('Time (s)', 'FontSize', fontSz);
    if run_clusters == 1, ylabel('\DeltaF/F (z-score)', 'FontSize', fontSz); end
    title(sprintf('Cluster %d', run_clusters), 'FontSize', fontSz);
    axis square; box off;
    style_axes(gca, fontSz);
end

%% ── Figure F: cluster proportions per day (pie charts) ───────────────────
figure;
set(gcf, 'Position', get(0, 'Screensize'));
colormap(cluster_cmap);
for run_days_ext = 1:n_days_ext
    curr_day = days(run_days_ext);

    subplot(2, n_days_ext, run_days_ext); hold on;
    idx_sa    = (cluster_info.day_of_neuron_g == curr_day) & ismember(cluster_info.mouse_nums_g, treatment_idx.Saline);
    counts_sa = arrayfun(@(c) sum(cluster_info.clust_id0(idx_sa) == c), 1:k_clust);
    pct_sa    = 100 * counts_sa / sum(counts_sa);
    p = pie(counts_sa);
    txt = findobj(p, 'Type', 'Text');
    for run_txt = 1:numel(txt), txt(run_txt).String = sprintf('%.1f%%', pct_sa(run_txt)); end
    set(findobj(p, 'Type', 'Patch'), 'EdgeColor', 'none');
    axis square;
    title(sprintf('Saline  Day %d', curr_day), 'FontSize', fontSz);

    subplot(2, n_days_ext, run_days_ext + n_days_ext); hold on;
    idx_md    = (cluster_info.day_of_neuron_g == curr_day) & ismember(cluster_info.mouse_nums_g, treatment_idx.MDMA);
    counts_md = arrayfun(@(c) sum(cluster_info.clust_id0(idx_md) == c), 1:k_clust);
    pct_md    = 100 * counts_md / sum(counts_md);
    p = pie(counts_md);
    txt = findobj(p, 'Type', 'Text');
    for run_txt = 1:numel(txt), txt(run_txt).String = sprintf('%.1f%%', pct_md(run_txt)); end
    set(findobj(p, 'Type', 'Patch'), 'EdgeColor', 'none');
    axis square;
    title(sprintf('MDMA  Day %d', curr_day), 'FontSize', fontSz);
end

%% ── Analysis: cluster identity stability matrix ───────────────────────────
maintained_matrix_shared = nan(n_days_ext, n_days_ext, k_clust, 2);  % (dayA, dayB, cluster, cohort)
cluster_order = [1, 2, 4, 3];

for run_treatments = 1:2
    mice = treatment_idx.(treatments{run_treatments});
    for run_clusters = cluster_order
        per_mouse_maintaineds = nan(n_days_ext, n_days_ext, numel(mice));
        for run_mice_in_group = 1:numel(mice)
            curr_mouse = mice(run_mice_in_group);
            for run_days_A = 1:n_days_ext
                for run_days_B = 1:n_days_ext
                    dayA = days(run_days_A);  dayB = days(run_days_B);
                    neurons_A = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==curr_mouse & cluster_info.day_of_neuron_g==dayA);
                    neurons_B = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==curr_mouse & cluster_info.day_of_neuron_g==dayB);
                    shared    = intersect(neurons_A, neurons_B);
                    if isempty(shared)
                        per_mouse_maintaineds(run_days_A, run_days_B, run_mice_in_group) = NaN;
                        continue;
                    end
                    idx_A = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==curr_mouse & cluster_info.day_of_neuron_g==dayA & cluster_info.clust_id0==run_clusters & ismember(cluster_info.neuron_idx_g, shared));
                    idx_B = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==curr_mouse & cluster_info.day_of_neuron_g==dayB & cluster_info.clust_id0==run_clusters & ismember(cluster_info.neuron_idx_g, shared));
                    if isempty(idx_A)
                        per_mouse_maintaineds(run_days_A, run_days_B, run_mice_in_group) = NaN;
                    else
                        per_mouse_maintaineds(run_days_A, run_days_B, run_mice_in_group) = mean(ismember(idx_A, idx_B));
                    end
                end
            end
        end
        maintained_matrix_shared(:, :, run_clusters, run_treatments) = mean(per_mouse_maintaineds, 3, 'omitnan');
    end
end

%% ── Analysis: day-pair cluster identity overlap ───────────────────────────
day_pairs   = [2 3; 3 4; 4 5];
pair_labels = {'2–3', '3–4', '4–5'};
n_pairs     = size(day_pairs, 1);
pair_maintained_shared = nan(n_pairs, k_clust, total_number_of_mice);

for run_mice = 1:total_number_of_mice
    for run_clusters = cluster_order
        for run_pairs = 1:n_pairs
            dayA = day_pairs(run_pairs, 1);  dayB = day_pairs(run_pairs, 2);
            if ~ismember(dayA, days) || ~ismember(dayB, days), continue; end
            neurons_A = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==run_mice & cluster_info.day_of_neuron_g==dayA);
            neurons_B = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==run_mice & cluster_info.day_of_neuron_g==dayB);
            shared    = intersect(neurons_A, neurons_B);
            if isempty(shared), continue; end
            idx_A = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==run_mice & cluster_info.day_of_neuron_g==dayA & cluster_info.clust_id0==run_clusters & ismember(cluster_info.neuron_idx_g, shared));
            idx_B = cluster_info.neuron_idx_g(cluster_info.mouse_nums_g==run_mice & cluster_info.day_of_neuron_g==dayB & cluster_info.clust_id0==run_clusters & ismember(cluster_info.neuron_idx_g, shared));
            if ~isempty(idx_A)
                pair_maintained_shared(run_pairs, run_clusters, run_mice) = mean(ismember(idx_A, idx_B));
            end
        end
    end
end

%% ── Figure G: cluster identity stability across day pairs ─────────────────
figure;
set(gcf, 'Position', get(0, 'Screensize'));
for run_cluster_order_idx = 1:k_clust
    curr_cluster = cluster_order(run_cluster_order_idx);
    subplot(1, k_clust, run_cluster_order_idx); hold on;
    title(sprintf('Cluster %d', curr_cluster), 'FontSize', fontSz);
    for run_pairs = 1:n_pairs
        saline_vals = squeeze(pair_maintained_shared(run_pairs, curr_cluster, treatment_idx.Saline));
        mdma_vals   = squeeze(pair_maintained_shared(run_pairs, curr_cluster, treatment_idx.MDMA));
        errorbar(run_pairs-0.15, mean(saline_vals,'omitnan'), std(saline_vals,'omitnan')/sqrt(sum(~isnan(saline_vals))), ...
            'o', 'Color', Saline_color, 'MarkerFaceColor', Saline_color, 'CapSize', 6, 'LineWidth', 1.5, 'MarkerSize', 7);
        errorbar(run_pairs+0.15, mean(mdma_vals,  'omitnan'), std(mdma_vals,  'omitnan')/sqrt(sum(~isnan(mdma_vals))), ...
            'o', 'Color', MDMA_color,   'MarkerFaceColor', MDMA_color,   'CapSize', 6, 'LineWidth', 1.5, 'MarkerSize', 7);
    end
    set(gca, 'XTick', 1:n_pairs, 'XTickLabel', pair_labels, 'FontSize', fontSz);
    xlim([0.5 3.5]); ylim([0 1]);
    if run_cluster_order_idx == 1
        ylabel('% maintained (shared)', 'FontSize', fontSz);
        legend({'Saline', 'MDMA'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', fontSz);
    end
    axis square; box off;
    style_axes(gca, fontSz);
end

%% ── Statistics G: cluster identity stability (mixed-effects model) ────────
for run_clusters = 1:k_clust
    mouse_id_vec = [];  pair_id_vec = [];  cohort_vec = {};  value_vec = [];
    for run_mice = 1:total_number_of_mice
        for run_pairs = 1:n_pairs
            v = pair_maintained_shared(run_pairs, run_clusters, run_mice);
            if isnan(v), continue; end
            mouse_id_vec(end+1, 1) = run_mice;
            pair_id_vec(end+1,  1) = run_pairs;
            value_vec(end+1,    1) = v;
            if ismember(run_mice, treatment_idx.MDMA)
                cohort_vec{end+1, 1} = 'MDMA';
            else
                cohort_vec{end+1, 1} = 'Saline';
            end
        end
    end
    tbl = table(categorical(mouse_id_vec), categorical(pair_id_vec), categorical(cohort_vec), value_vec, ...
        'VariableNames', {'Mouse','Pair','Cohort','Value'});
    lme = fitlme(tbl, 'Value ~ Pair*Cohort + (1|Mouse)');
    disp(['Cluster ' num2str(run_clusters)]);
    disp(anova(lme));
end

%% ── Analysis: cluster tuning-curve correlations by day difference (panel H)
% Uses the same filtered cluster_info fields as panels F and G
% (mouse_nums_g, day_of_neuron_g, clust_id0, neuron_idx_g) — these only
% contain neurons that passed the variance filter and have valid cluster IDs.
number_of_clusters = cluster_info.n_clusters;
cue_duration_frames = cue_onset_frame : cue_offset_frame;
cluster_delta = cell(number_of_clusters,total_number_of_mice,4);
cluster_adj = cell(number_of_clusters,total_number_of_mice,n_adj);

for cluster = 1:number_of_clusters
    for run_mice = 1:total_number_of_mice
        cue_traces = []; % bins x cues
        cue_days = []; % 1 x cues
        for run_days = 1:n_days_total
            % for run_days = 1:n_days
            curr_day_traces = all_traces{run_mice,run_days};
            if isempty(curr_day_traces), continue; end
            this_cluster_cell_idx = (cluster_info.mouse_id==run_mice) & (cluster_info.day==run_days) & (cluster_info.cluster_full==cluster);
            neurons = cluster_info.neuron_index(this_cluster_cell_idx); % neurons in this cluster
            if isempty(neurons), continue; end % check


            for run_cues = 1:size(curr_day_traces,3)   % loop through each cue in this day
                pop_mat = []; % will be neurons x bins 
                for run_IDs = neurons(:)'
                    curr_cue_traces = squeeze(curr_day_traces(run_IDs, cue_duration_frames, run_cues))';
                    ds_cue_traces = downsample_trace(curr_cue_traces, frame_rate, 1);
                    if isempty(ds_cue_traces), continue; end
                    pop_mat = [pop_mat; ds_cue_traces];
                end

                if size(pop_mat,1) >= 2
                    pop_vec = mean(pop_mat,1);   % take the mean trace of all neurons in this cluster/mouse/day
                    cue_traces = [cue_traces, pop_vec(:)];
                    cue_days   = [cue_days, run_days];
                end
            end

        end
        if size(cue_traces,2) < 2, continue; end

        % -- calculate all correlations --
        R = corr(cue_traces,'rows','pairwise');  % pairwise correlation matrix for all cues and days (cuexday  x cuexday)
        R(eye(size(R))==1) = NaN; % set the diagonal to NaN so that the within day stuff isnt inflated by same to same

        % -- get Rs across delta day pairs --
        for i = 1:numel(cue_days)
            for j = i+1:numel(cue_days)
                delta = cue_days(j) - cue_days(i);
                if delta >= 0 && delta <= 3
                    cluster_delta{cluster,run_mice,delta+1} = [cluster_delta{cluster,run_mice,delta+1}, R(i,j)];
                end
            end
        end


    end
end

%% ── Figure H: cluster tuning-curve correlations vs day difference ─────────
figure;
set(gcf, 'Position', get(0, 'Screensize'));
for run_cluster_order_idx = 1:k_clust
    curr_cluster = cluster_order(run_cluster_order_idx);
    subplot(1, k_clust, run_cluster_order_idx); hold on;

    saline_vals = nan(numel(treatment_idx.Saline), 4);
    mdma_vals   = nan(numel(treatment_idx.MDMA),   4);

    for run_mice_in_group = 1:numel(treatment_idx.Saline)
        curr_mouse = treatment_idx.Saline(run_mice_in_group);
        for run_delta_bins = 1:4
            saline_vals(run_mice_in_group, run_delta_bins) = mean(cluster_delta{curr_cluster, curr_mouse, run_delta_bins}, 'omitnan');
        end
    end
    for run_mice_in_group = 1:numel(treatment_idx.MDMA)
        curr_mouse = treatment_idx.MDMA(run_mice_in_group);
        for run_delta_bins = 1:4
            mdma_vals(run_mice_in_group, run_delta_bins) = mean(cluster_delta{curr_cluster, curr_mouse, run_delta_bins}, 'omitnan');
        end
    end

    errorbar(0:3, mean(saline_vals,'omitnan'), std(saline_vals,'omitnan')/sqrt(size(saline_vals,1)), ...
        '-o', 'Color', Saline_color, 'LineWidth', 2, 'MarkerFaceColor', Saline_color);
    errorbar(0:3, mean(mdma_vals,  'omitnan'), std(mdma_vals,  'omitnan')/sqrt(size(mdma_vals,1)), ...
        '-o', 'Color', MDMA_color,   'LineWidth', 2, 'MarkerFaceColor', MDMA_color);

    xticks(0:3);
    xlabel('Days difference', 'FontSize', fontSz);
    ylim([-0.1 0.4]);
    title(sprintf('Cluster %d', curr_cluster), 'FontSize', fontSz);
    if run_cluster_order_idx == 1
        ylabel('Tuning curve correlation', 'FontSize', fontSz);
        legend({'Saline','MDMA'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', fontSz);
    end
    axis square; box off;
    style_axes(gca, fontSz);
end

%% ── Statistics H: t-tests per cluster and delta-day ─────────────────────
disp('── Panel H: t-tests between cohorts per delta-day value ──');
for run_clusters = 1:k_clust
    fprintf('\n--- Cluster %d ---\n', run_clusters);
    pvals = nan(1, 4);
    for run_delta_bins = 1:4
        saline_vec = nan(numel(treatment_idx.Saline), 1);
        mdma_vec   = nan(numel(treatment_idx.MDMA),   1);
        for run_mice_in_group = 1:numel(treatment_idx.Saline)
            saline_vec(run_mice_in_group) = mean(cluster_delta{run_clusters, treatment_idx.Saline(run_mice_in_group), run_delta_bins}, 'omitnan');
        end
        for run_mice_in_group = 1:numel(treatment_idx.MDMA)
            mdma_vec(run_mice_in_group) = mean(cluster_delta{run_clusters, treatment_idx.MDMA(run_mice_in_group), run_delta_bins}, 'omitnan');
        end
        [~, pvals(run_delta_bins)] = ttest2(mdma_vec, saline_vec);
    end
    fprintf('Delta-day p-values (delta 0–3):\n');
    disp(pvals);
end

%% ── Statistics H: LME (cohort × delta interaction) per cluster ───────────
disp('── Panel H: LME (cohort × delta interaction) ──');
for run_clusters = 1:k_clust
    fprintf('\n--- Cluster %d ---\n', run_clusters);
    y_vec = [];  delta_vec = [];  cohort_vec = {};  mouse_vec = [];

    for run_delta_bins = 1:4
        day_delta = run_delta_bins - 1;
        for run_mice_in_group = 1:numel(treatment_idx.Saline)
            curr_mouse = treatment_idx.Saline(run_mice_in_group);
            val = mean(cluster_delta{run_clusters, curr_mouse, run_delta_bins}, 'omitnan');
            if isnan(val), continue; end
            y_vec(end+1,     1) = val;
            delta_vec(end+1, 1) = day_delta;
            cohort_vec{end+1,1} = 'Saline';
            mouse_vec(end+1, 1) = curr_mouse;
        end
        for run_mice_in_group = 1:numel(treatment_idx.MDMA)
            curr_mouse = treatment_idx.MDMA(run_mice_in_group);
            val = mean(cluster_delta{run_clusters, curr_mouse, run_delta_bins}, 'omitnan');
            if isnan(val), continue; end
            y_vec(end+1,     1) = val;
            delta_vec(end+1, 1) = day_delta;
            cohort_vec{end+1,1} = 'MDMA';
            mouse_vec(end+1, 1) = curr_mouse;
        end
    end

    tbl = table(y_vec, categorical(delta_vec), categorical(cohort_vec), categorical(mouse_vec), ...
        'VariableNames', {'y','delta','cohort','mouse'});
    lme = fitlme(tbl, 'y ~ cohort*delta + (1|mouse)');
    disp(anova(lme));
end

%% ── Helper functions ─────────────────────────────────────────────────────
function style_axes(ax, font_sz)
% STYLE_AXES  Apply standard publication-quality axis formatting.
ax.FontSize        = font_sz;
ax.LineWidth       = 1;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
ax.XAxis.Color     = 'k';
ax.YAxis.Color     = 'k';
end

function [ax1, ax2, cb_pos] = heatmap_axes_layout()
% HEATMAP_AXES_LAYOUT  Create two side-by-side axes with shared colorbar position.
left_margin = 0.06;  right_margin = 0.04;
bottom      = 0.12;  top_margin   = 0.08;
gap_between = 0.06;  cb_width     = 0.015;  cb_gap = 0.01;

avail_width = 1 - left_margin - right_margin - cb_width - gap_between - 2*cb_gap;
main_width  = avail_width / 2;
main_height = 1 - bottom - top_margin;

ax1_pos = [left_margin,                            bottom, main_width, main_height];
ax2_pos = [left_margin + main_width + gap_between, bottom, main_width, main_height];
cb_pos  = [ax2_pos(1) + ax2_pos(3) + cb_gap,       bottom + 0.02, cb_width, main_height - 0.04];

ax1 = axes('Position', ax1_pos);
ax2 = axes('Position', ax2_pos);
end
