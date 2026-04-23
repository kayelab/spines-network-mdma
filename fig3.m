%% ── Load data ─────────────────────────────────────────────────────────────
load([pwd '\data\untethered_freezing.mat']);
load([pwd '\data\cue_traces.mat']);
load([pwd '\data\ITI_traces.mat']);
load([pwd '\data\baselined_freezing.mat']);
all_traces = [cue_traces.MDMA; cue_traces.Saline];  % concatenate cohorts
all_freeze = [pre_cue_corrected_freezing.MDMA;pre_cue_corrected_freezing.Saline];
all_ITI_vecs = vertcat(ITI_traces.Saline, ITI_traces.MDMA);
all_untethered_freezing = [untethered_freezing.Saline; untethered_freezing.MDMA];
%% ── Color palette ────────────────────────────────────────────────────────
MDMA_color   = [173,  84,  50] / 255;
Saline_color = [ 88,  88,  90] / 255;
ITI_color   = [197, 166, 112] / 255;
fontSz     = 14;

%% set paramenters
treatments = {'Saline'; 'MDMA'};
cohort_idx = struct('MDMA',1:8, 'Saline',9:16);
n_mice_per_group = 8;
days = 2:5;
n_days = numel(days);

%
frame_rate      = 30;   % Hz
time_before_cue = 30;   % s
que_time        = 30;   % s (cue duration)
time_after_cue_onset  = 60;   % s

bin_sec = 2;
bin_frames = round(bin_sec * frame_rate);
clip_trials = true; clip_trial_idx = 1:10;
smooth_win = round(0.8 * frame_rate);

n_frames_total = frame_rate * (time_before_cue +  time_after_cue_onset)+1;
cue_onset      = frame_rate * time_before_cue;
cue_offset     = cue_onset  + frame_rate * que_time;

xtick_diff = 10;
xticks_vec = ( - time_before_cue) : xtick_diff : ( + time_after_cue_onset);

AUC_frames              = (cue_onset + frame_rate*7) : (cue_offset + frame_rate*7);
initial_response_frames =  cue_onset : (cue_onset + 7*frame_rate);

%%
for run_treatments = 1:2
    curr_treatment = treatments{run_treatments};
    mean_BL_reduced_freezing.(curr_treatment)= nan(n_mice_per_group,n_frames_total,n_days);
    mean_offset_freezing.(curr_treatment) = nan(n_mice_per_group,n_days);
    mean_gcamp_activity_mats.(curr_treatment) = nan(n_mice_per_group,n_frames_total,n_days);
    AUC.(curr_treatment) = nan(n_mice_per_group,n_days);
    for run_mice = 1:n_mice_per_group
        curr_mouse_behavior_traces = pre_cue_corrected_freezing.(curr_treatment)(run_mice,days);
        curr_mouse_gcamp_traces = cue_traces.(curr_treatment)(run_mice,days);
        for run_days = 1:n_days
            curr_freezing_vec = mean(curr_mouse_behavior_traces{run_days},3,'omitnan');
            mean_freezing_mats.(curr_treatment)(run_mice,:,run_days) = curr_freezing_vec;
            mean_offset_freezing.(curr_treatment)(run_mice,run_days) = mean(curr_freezing_vec(cue_offset-10*frame_rate:cue_offset+10*frame_rate));
            curr_gcamp_activity_vec = mean(mean(curr_mouse_gcamp_traces{run_days},3,'omitnan'),'omitnan');
            mean_gcamp_activity_mats.(curr_treatment)(run_mice,:,run_days) = curr_gcamp_activity_vec;
            if ~all(isnan(curr_gcamp_activity_vec))
                AUC.(curr_treatment)(run_mice,run_days) = trapz(curr_gcamp_activity_vec(AUC_frames));
            end
        end
    end
end

ITI_mean_activity_vecs = nan(n_mice_per_group*2,n_frames_total,n_days);
AUC.ITI = nan(n_mice_per_group*2,n_days);
for run_mice = 1:2*n_mice_per_group
    curr_mouse_ITI_data = all_ITI_vecs(run_mice,days);
    for run_days = 1:n_days
        curr_gcamp_activity_vec = mean(mean(curr_mouse_ITI_data{run_days},3,'omitnan'),'omitnan');
        ITI_mean_activity_vecs(run_mice,:,run_days) = curr_gcamp_activity_vec;
        if ~all(isnan(curr_gcamp_activity_vec))
            AUC.ITI(run_mice,run_days) = trapz(curr_gcamp_activity_vec(AUC_frames));
        end
    end
end

%% plot panel B
n_mice_per_group_u = size(untethered_freezing.Saline,1);
number_of_days = size(untethered_freezing.Saline,2);
figure
set(gcf, 'Position', get(0, 'Screensize'));
hold on
S_mean = mean(untethered_freezing.Saline,'omitnan');
S_sem = std(untethered_freezing.Saline)/sqrt(n_mice_per_group_u);
M_mean = mean(untethered_freezing.MDMA,'omitnan');
M_sem = std(untethered_freezing.MDMA)/sqrt(n_mice_per_group_u);

plot(days,S_mean,'.-','Color',Saline_color,'LineWidth',2,'MarkerSize',35)
plot(days,M_mean,'.-','Color',MDMA_color,'LineWidth',2,'MarkerSize',35)

EB_S = errorbar(days,S_mean,S_sem,'Color',Saline_color,'LineWidth',2);
EB_M = errorbar(days,M_mean,M_sem,'Color',MDMA_color,'LineWidth',2);

legend({'Saline';'MDMA'}, 'Location', 'northeast', 'FontSize', 24, 'Box', 'off');

xlim([1.9 5.1])
ylim([0 0.7])
xlabel('Day')
ylabel('Cue evoked freezing (%)')

ax = gca;
style_axes(ax, 25);

%  Repeated-Measures ANOVA for untethered freezing 
dayLabels_B = {'Day2','Day3','Day4','Day5'};
n_mice_per_group_u = size(untethered_freezing.Saline, 1);  % use actual size (12)
treatment_B = [repmat({'Saline'}, n_mice_per_group_u, 1); repmat({'MDMA'}, n_mice_per_group_u, 1)];
mouseID_B   = (1:2*n_mice_per_group_u)';

T_B = array2table(all_untethered_freezing, 'VariableNames', dayLabels_B);
T_B.Treatment = categorical(treatment_B);
T_B.MouseID   = mouseID_B;

withinDesign_B = table((1:number_of_days)', 'VariableNames', {'Day'});
withinDesign_B.Day = categorical(withinDesign_B.Day);

rm_B = fitrm(T_B, 'Day2,Day3,Day4,Day5 ~ Treatment', 'WithinDesign', withinDesign_B);
ranovatbl_B = ranova(rm_B, 'WithinModel', 'Day');
disp('=== Repeated-Measures ANOVA Table for untethered freezing ===');
disp(ranovatbl_B);

%% plot panels F-G
figure;
set(gcf, 'Position', get(0, 'Screensize'));
for run_days = 1:n_days % panel F
    subplot(1,n_days+1,run_days)
    hold on
    S_mean = smooth_vector_data(mean(mean_freezing_mats.Saline(:,:,run_days),'omitnan'),smooth_win);
    S_sem = smooth_vector_data(std(mean_freezing_mats.Saline(:,:,run_days),'omitnan'),smooth_win)/sqrt(n_mice_per_group);
    M_mean = smooth_vector_data(mean(mean_freezing_mats.MDMA(:,:,run_days),'omitnan'),smooth_win);
    M_sem = smooth_vector_data(std(mean_freezing_mats.MDMA(:,:,run_days),'omitnan'),smooth_win)/sqrt(n_mice_per_group);

    plot(1:n_frames_total, S_mean,'Color',Saline_color)
    plot(1:n_frames_total, M_mean,'Color',MDMA_color)

    EB_S   = shadedErrorBar(1:n_frames_total, S_mean,   S_sem,   'k', 1);
    EB_M   = shadedErrorBar(1:n_frames_total, M_mean,   M_sem,   'k', 1);

    EB_S.mainLine.Color   = Saline_color;  EB_S.patch.FaceColor   = Saline_color;
    EB_M.mainLine.Color   = MDMA_color;    EB_M.patch.FaceColor   = MDMA_color;

    ylim([-0.4 0.7])
    YL = ylim;
    plot([1 1]*cue_onset,                        YL, 'k--', 'LineWidth', 1.1)
    plot([1 1]*cue_offset, YL, 'k--', 'LineWidth', 1.1)
    ylim(YL);

    if run_days == 1
        bg = patch([400 2400 2400 400], [YL(1) YL(1) YL(2) YL(2)], ...
            [241 179 153]/255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end

    x1 = cue_offset - 10*frame_rate;
    x2 = cue_offset + 10*frame_rate;
    rect = patch([x1 x2 x2 x1], [YL(1) YL(1) YL(2) YL(2)], ...
        [141 215 242]/255, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    uistack(rect, 'bottom');  % blue area around offset
    if run_days == 1
        uistack(bg, 'bottom');  % orange day 2 - drug
    end

    title(['Day ' num2str(run_days+1)]);
    xlabel('Time from cue onset');
    box off;

    ax = gca;
    ax.XTick      = 0 : frame_rate*xtick_diff : n_frames_total;
    ax.XTickLabel = xticks_vec;
    xlim([400 2400]);
    ylabel('\DeltaFreezing (% from BL)');
    style_axes(ax, 14);

end
% panel G
subplot (1,n_days+1,run_days+1)
hold on;

mean_sal = mean(mean_offset_freezing.Saline,'omitnan');
sem_sal = std(mean_offset_freezing.Saline,'omitnan')/sqrt(n_mice_per_group);
mean_mdma = mean(mean_offset_freezing.MDMA,'omitnan');
sem_mdma = std(mean_offset_freezing.MDMA,'omitnan')/sqrt(n_mice_per_group);


plot(1:4,         mean_sal,   '.', 'Color', Saline_color, 'MarkerSize', 45);
plot((1:4)-0.008, mean_mdma,  '.', 'Color', MDMA_color,   'MarkerSize', 45);
errorbar(1:4,         mean_sal,   sem_sal,   'Color', Saline_color, 'LineWidth', 2, 'CapSize', 15);
errorbar((1:4)-0.008, mean_mdma,  sem_mdma,  'Color', MDMA_color,   'LineWidth', 2, 'CapSize', 15);

xticks(1:4);  xticklabels(num2cell(days));
xlim([0.9 4.1]);
xlabel('Day',        'FontSize', fontSz);
ylabel('Freezing mean % from BL', 'FontSize', fontSz);
legend({'Saline','MDMA'}, 'Location', 'northeast', 'FontSize', 24, 'Box', 'off');
style_axes(gca, fontSz);

% Statistics panel G — RM-ANOVA on late-cue freezing 

dayLabels_G = {'Day2','Day3','Day4','Day5'};

data_G      = [mean_offset_freezing.Saline; mean_offset_freezing.MDMA];  % (2*n_mice) x n_days
treatment_G = [repmat({'Saline'}, n_mice_per_group, 1); ...
    repmat({'MDMA'},   n_mice_per_group, 1)];
mouseID_G   = (1 : 2*n_mice_per_group)';

T_G = array2table(data_G, 'VariableNames', dayLabels_G);
T_G.Treatment = categorical(treatment_G);
T_G.MouseID   = mouseID_G;

withinDesign_G     = table((1:number_of_days)', 'VariableNames', {'Day'});
withinDesign_G.Day = categorical(withinDesign_G.Day);

rm_G       = fitrm(T_G, 'Day2,Day3,Day4,Day5 ~ Treatment', 'WithinDesign', withinDesign_G);
ranova_G   = ranova(rm_G, 'WithinModel', 'Day');
between_G  = anova(rm_G);

disp('=== Panel G: RM-ANOVA on late-cue freezing (±10 s around cue offset) ===');
disp(ranova_G);

%% plots panel H-I
figure;
set(gcf, 'Position', get(0, 'Screensize'));
for run_days = 1:n_days % panel H
    subplot(1,n_days+1,run_days)
    hold on
    S_mean = smooth_vector_data(mean(mean_gcamp_activity_mats.Saline(:,:,run_days),'omitnan'),smooth_win);
    S_sem = smooth_vector_data(std(mean_gcamp_activity_mats.Saline(:,:,run_days),'omitnan'),smooth_win)/sqrt(n_mice_per_group);
    M_mean = smooth_vector_data(mean(mean_gcamp_activity_mats.MDMA(:,:,run_days),'omitnan'),smooth_win);
    M_sem = smooth_vector_data(std(mean_gcamp_activity_mats.MDMA(:,:,run_days),'omitnan'),smooth_win)/sqrt(n_mice_per_group);
    ITI_mean = smooth_vector_data(mean(ITI_mean_activity_vecs(:,:,run_days),'omitnan'),smooth_win);
    ITI_sem = smooth_vector_data(std(ITI_mean_activity_vecs(:,:,run_days),'omitnan'),smooth_win)/sqrt(n_mice_per_group*2);

    plot(1:n_frames_total, S_mean,'Color',Saline_color)
    plot(1:n_frames_total, M_mean,'Color',MDMA_color)
    plot(1:n_frames_total, ITI_mean,'Color',ITI_color)

    EB_S   = shadedErrorBar(1:n_frames_total, S_mean,   S_sem,   'k', 1);
    EB_M   = shadedErrorBar(1:n_frames_total, M_mean,   M_sem,   'k', 1);
    EB_ITI = shadedErrorBar(1:n_frames_total, ITI_mean, ITI_sem, 'k', 1);

    EB_S.mainLine.Color   = Saline_color;  EB_S.patch.FaceColor   = Saline_color;
    EB_M.mainLine.Color   = MDMA_color;    EB_M.patch.FaceColor   = MDMA_color;
    EB_ITI.mainLine.Color = ITI_color;    EB_ITI.patch.FaceColor = ITI_color;

    YL = ylim;
    plot([1 1]*cue_onset,                        YL, 'k--', 'LineWidth', 1.1)
    plot([1 1]*cue_offset, YL, 'k--', 'LineWidth', 1.1)
    ylim(YL);

    if run_days == 1
        bg = patch([400 2400 2400 400], [YL(1) YL(1) YL(2) YL(2)], ...
            [241 179 153]/255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end

    x1 = cue_offset - 10*frame_rate;
    x2 = cue_offset + 10*frame_rate;
    rect = patch([x1 x2 x2 x1], [YL(1) YL(1) YL(2) YL(2)], ...
        [141 215 242]/255, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    uistack(rect, 'bottom');  % blue just above orange
    if run_days == 1
        uistack(bg, 'bottom');  % orange goes to very bottom
    end

    title(['Day ' num2str(run_days+1)]);
    % % legend({'Saline','MDMA','ITI'}, 'Location', 'southoutside');
    xlabel('Time from cue onset');
    box off;

    ax = gca;
    ax.XTick      = 0 : frame_rate*xtick_diff : n_frames_total;
    ax.XTickLabel = xticks_vec;
    xlim([400 2400]);
    ylabel('dF/F0 (change from baseline)');
    style_axes(ax, 14);

end
% panel I
subplot(1,n_days+1,n_days+1)
hold on;

mean_sal = mean(AUC.Saline,'omitnan');
sem_sal = std(AUC.Saline,'omitnan')/sqrt(n_mice_per_group);
mean_mdma = mean(AUC.MDMA,'omitnan');
sem_mdma = std(AUC.MDMA,'omitnan')/sqrt(n_mice_per_group);
mean_ITI  = mean(AUC.ITI,'omitnan');
sem_ITI = std(AUC.ITI,'omitnan')/sqrt(n_mice_per_group*2);


plot(1:4,         mean_sal,   '.', 'Color', Saline_color, 'MarkerSize', 45);
plot((1:4)-0.008, mean_mdma,  '.', 'Color', MDMA_color,   'MarkerSize', 45);
plot((1:4)+0.008, mean_ITI, '.', 'Color', ITI_color,   'MarkerSize', 45);
errorbar(1:4,         mean_sal,   sem_sal,   'Color', Saline_color, 'LineWidth', 2, 'CapSize', 15);
errorbar((1:4)-0.008, mean_mdma,  sem_mdma,  'Color', MDMA_color,   'LineWidth', 2, 'CapSize', 15);
errorbar((1:4)+0.008, mean_ITI, sem_ITI, 'Color', ITI_color,   'LineWidth', 2, 'CapSize', 15);

xticks(1:4);  xticklabels(num2cell(days));
xlim([0.9 4.1]);
xlabel('Day',        'FontSize', fontSz);
ylabel('AUC (a.u.)', 'FontSize', fontSz);
legend({'Saline','MDMA','ITI'}, 'Location', 'northeast', 'FontSize', 14, 'Box', 'off');

style_axes(gca, fontSz);

%% plots panels J-K
% Analysis window: ±10 s around cue offset (the blue-shaded region in F-H).
% This matches cue_offset - 10*frame_rate : cue_offset + 10*frame_rate,
jk_win_frames = (cue_offset - 10*frame_rate) : (cue_offset + 10*frame_rate);

% -- compute per-trial window averages (activity and freezing) --
max_trials = 10;
window_activity = struct();
window_freezing = struct();
for run_treatments = 1:2
    curr_treatment = treatments{run_treatments};
    window_activity.(curr_treatment) = cell(n_mice_per_group, n_days, max_trials);
    window_freezing.(curr_treatment) = cell(n_mice_per_group, n_days, max_trials);
    mice = cohort_idx.(curr_treatment);
    for run_mice = 1:n_mice_per_group
        curr_mouse = mice(run_mice);
        for run_days = 1:n_days
            day_abs = days(run_days);           % absolute day index into all_traces
            X = all_traces{curr_mouse, day_abs};         % neurons × time × trials
            F = all_freeze{curr_mouse, day_abs};         % 1 × time × trials

            if isempty(X) || all(isnan(X(:))), continue; end

            n_trials  = size(X, 3);
            idx_win   = jk_win_frames(jk_win_frames <= size(X, 2));  % guard bounds
            trials_to_use = clip_trial_idx(clip_trial_idx <= n_trials);

            for run_trials = 1:numel(trials_to_use)
                curr_trial = trials_to_use(run_trials);
                pop_trace  = mean(X(:, :, curr_trial), 1, 'omitnan');
                freeze_vec = F(1, :, curr_trial);
                window_activity.(curr_treatment){run_mice, run_days, run_trials} = mean(pop_trace(idx_win),  'omitnan');
                window_freezing.(curr_treatment){run_mice, run_days, run_trials} = mean(freeze_vec(idx_win), 'omitnan');
            end
        end
    end
end

% -- correlate activity and freezing at the mouse level --
results_jk = struct();
for run_days = 1:n_days
    for run_treatments = 1:2
        curr_treatment = treatments{run_treatments};
        acts = nan(n_mice_per_group, 1);
        frzs = nan(n_mice_per_group, 1);
        for run_mice = 1:n_mice_per_group
            acts(run_mice) = mean(cell2mat(squeeze(window_activity.(curr_treatment)(run_mice, run_days, :))), 'omitnan');
            frzs(run_mice) = mean(cell2mat(squeeze(window_freezing.(curr_treatment)(run_mice, run_days, :))), 'omitnan');
        end
        results_jk(run_days).(curr_treatment).acts_x = acts;
        results_jk(run_days).(curr_treatment).frzs_y = frzs;

        valid  = ~isnan(acts) & ~isnan(frzs);
        acts_v = acts(valid);
        frzs_v = frzs(valid);
        results_jk(run_days).(curr_treatment).acts_x_valid = acts_v;
        results_jk(run_days).(curr_treatment).frzs_y_valid = frzs_v;

        mdl             = fitlm(acts_v, frzs_v);
        [r, p]          = corr(acts_v, frzs_v, 'rows', 'complete');
        ci              = coefCI(mdl);
        results_jk(run_days).(curr_treatment).mdl     = mdl;
        results_jk(run_days).(curr_treatment).r       = r;
        results_jk(run_days).(curr_treatment).p       = p;
        results_jk(run_days).(curr_treatment).beta    = mdl.Coefficients.Estimate(2);
        results_jk(run_days).(curr_treatment).beta_lo = ci(2, 1);
        results_jk(run_days).(curr_treatment).beta_hi = ci(2, 2);
    end
    all_acts_x = [results_jk(run_days).MDMA.acts_x; results_jk(run_days).Saline.acts_x];
    results_jk(run_days).xlims = [min(all_acts_x) max(all_acts_x)];
end

% -- panel J: scatter + regression fit per day --
figure;
set(gcf, 'Position', get(0, 'Screensize'));
for run_days = 1:n_days
    subplot(1, n_days+1, run_days); hold on;
    for run_treatments = 1:2
        curr_treatment = treatments{run_treatments};
        switch curr_treatment
            case 'MDMA'
                color    = MDMA_color;
            case 'Saline'
                color    = Saline_color;
        end
        acts_v = results_jk(run_days).(curr_treatment).acts_x_valid;
        frzs_v = results_jk(run_days).(curr_treatment).frzs_y_valid;
        mdl    = results_jk(run_days).(curr_treatment).mdl;
        xfit   = linspace(results_jk(run_days).xlims(1), results_jk(run_days).xlims(2), 100)';
        [yfit, yCI] = predict(mdl, xfit);
        scatter(acts_v, frzs_v, 30, color, 'filled');
        plot(xfit, yfit, 'Color', color, 'LineWidth', 1.5);
        fill([xfit; flipud(xfit)], [yCI(:,1); flipud(yCI(:,2))], color, ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    xlim(results_jk(run_days).xlims); ylim([-0.3 1]);
    YL = ylim;
    if run_days == 1
        bg = patch([results_jk(run_days).xlims(1) results_jk(run_days).xlims(2) ...
                    results_jk(run_days).xlims(2) results_jk(run_days).xlims(1)], ...
                   [YL(1) YL(1) YL(2) YL(2)], ...
                   [241 179 153]/255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        uistack(bg, 'bottom');  % orange day 2 - drug
        xlabel('Population activity (±10 s cue offset)', 'FontSize', fontSz);
        ylabel('Freezing (±10 s cue offset)',            'FontSize', fontSz);
    end
    title(sprintf('Day %d\nMDMA r=%.2f p=%.3f\nSaline r=%.2f p=%.3f', ...
        days(run_days), ...
        results_jk(run_days).MDMA.r,   results_jk(run_days).MDMA.p, ...
        results_jk(run_days).Saline.r, results_jk(run_days).Saline.p), ...
        'FontSize', fontSz);
    box off;
    style_axes(gca, fontSz);
end

% -- panel K: correlation coefficients and regression slopes across days --
r_mdma   = arrayfun(@(d) d.MDMA.r,   results_jk);
r_saline = arrayfun(@(d) d.Saline.r, results_jk);


subplot(1, n_days+1, n_days+1)
hold on
plot(days, r_mdma,   '-o', 'Color', MDMA_color,   'LineWidth', 2, 'MarkerFaceColor', MDMA_color);
plot(days, r_saline, '-o', 'Color', Saline_color, 'LineWidth', 2, 'MarkerFaceColor', Saline_color);
yline(0, '--k');
xlabel('Day',             'FontSize', fontSz);
ylabel('Correlation (r)', 'FontSize', fontSz);
title('Correlation coefficient', 'FontSize', fontSz);
legend({'MDMA','Saline'}, 'Location', 'northeast', 'FontSize', fontSz, 'Box', 'off');
xlim([days(1)-0.5 days(end)+0.5]); ylim([-1 1]);
box off;
style_axes(gca, fontSz);

figure
% subplot(1, 2, 2); hold on;
beta_mdma    = arrayfun(@(d) d.MDMA.beta,    results_jk);
beta_lo_mdma = arrayfun(@(d) d.MDMA.beta_lo, results_jk);
beta_hi_mdma = arrayfun(@(d) d.MDMA.beta_hi, results_jk);
beta_sal     = arrayfun(@(d) d.Saline.beta,    results_jk);
beta_lo_sal  = arrayfun(@(d) d.Saline.beta_lo, results_jk);
beta_hi_sal  = arrayfun(@(d) d.Saline.beta_hi, results_jk);

plot_shaded(days, beta_mdma, (beta_hi_mdma - beta_lo_mdma)/2, MDMA_color,   'FaceAlpha', 0.25, 'LineWidth', 1.5, 'Marker', 'o');
plot_shaded(days, beta_sal,  (beta_hi_sal  - beta_lo_sal )/2, Saline_color, 'FaceAlpha', 0.25, 'LineWidth', 1.5, 'Marker', 'o');
yline(0, '--k');
xlabel('Day',                      'FontSize', fontSz);
ylabel('Regression slope (\beta)', 'FontSize', fontSz);
title('Regression slopes, 95% CI', 'FontSize', fontSz);
xlim([days(1)-0.5 days(end)+0.5]);
box off;
style_axes(gca, fontSz);

%% statistics (J-K)
% Pool days 3–5 (run_days 2–4), i.e., the extinction testing days
% following the first MDMA-paired extinction session (day 2).
Acts_jk = []; Frzs_jk = []; Day_jk = []; Cohort_jk = {}; Mouse_jk = {};
for run_days = 2:n_days
    for run_treatments = 1:2
        curr_treatment = treatments{run_treatments};
        for run_mice = 1:n_mice_per_group
            a = mean(cell2mat(squeeze(window_activity.(curr_treatment)(run_mice, run_days, :))), 'omitnan');
            f = mean(cell2mat(squeeze(window_freezing.(curr_treatment)(run_mice, run_days, :))), 'omitnan');
            if ~isnan(a) && ~isnan(f)
                Acts_jk(end+1, 1)   = a;
                Frzs_jk(end+1, 1)   = f;
                Day_jk(end+1, 1)    = days(run_days);
                Cohort_jk{end+1, 1} = curr_treatment;
                Mouse_jk{end+1, 1}  = sprintf('%s_%d', curr_treatment, run_mice);
            end
        end
    end
end

tbl_jk = table(Acts_jk, Frzs_jk, Day_jk, categorical(Cohort_jk), categorical(Mouse_jk), ...
    'VariableNames', {'activity','freezing','day','cohort','mouse'});
lme_jk = fitlme(tbl_jk, 'freezing ~ activity*cohort*day + (1|mouse)');
disp(lme_jk); disp(anova(lme_jk));

%% ─────────────────────────────────────────────────────────────────────────
function style_axes(ax, font_sz)
% STYLE_AXES  Apply standard publication-quality axis formatting.
axis square
ax.FontSize        = font_sz;
ax.LineWidth       = 1;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
ax.XAxis.Color     = 'k';
ax.YAxis.Color     = 'k';
end