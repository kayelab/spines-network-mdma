function [cluster_info] = hierarchical_clustering(traces,smooth_win)

cluster_info = struct();

% all_traces = [traces.MDMA; traces.Saline]; % concatenate
[n_mice, n_days] = size(traces); % each all traces entry is neurons x time x trials

% --- create a 'megamouse' concatenated array of all mice and days ---
mega_traces_cell = cell(n_mice, n_days);
mouse_nums_cell = cell(n_mice, n_days);
day_cell = cell(n_mice, n_days);
neuron_idx_cell = cell(n_mice, n_days);
% filled_count = 0;
for run_days = 2:n_days
    for run_mice = 1:n_mice
        X = traces{run_mice, run_days};  % neurons × time × trials
        n_neurons = size(X,1);
        if isempty(X)
            continue;
        end
        Xmean = mean(X, 3); % actual data, neurons x time (avg over trials)
        Xmean = movmean(Xmean, smooth_win, 2);  % smooth
        mega_traces_cell{run_mice,run_days} = Xmean;
        mouse_nums_cell{run_mice,run_days} = run_mice * ones(n_neurons,1);
        day_cell{run_mice,run_days} = run_days * ones(n_neurons,1);
        neuron_idx_cell{run_mice,run_days} = (1:n_neurons)';
    end
end
mega_traces = vertcat(mega_traces_cell{:});
mouse_nums = vertcat(mouse_nums_cell{:});
day_of_neuron = vertcat(day_cell{:});
neuron_idx_within_mouse = vertcat(neuron_idx_cell{:});


% -- preprocess for clustering -- 
row_var = var(mega_traces,0,2);  % get rid of any zero-variance neurons
good = row_var > 0;
mega_traces_filt = mega_traces(good,:);
cluster_info.mouse_nums_g = mouse_nums(good);
cluster_info.day_of_neuron_g = day_of_neuron(good);
cluster_info.neuron_idx_g = neuron_idx_within_mouse(good);
days = unique(cluster_info.day_of_neuron_g);
days = days(~isnan(days));
days = sort(days(:)');
cluster_info.n_days_local = numel(days);

cluster_info.X0 = mega_traces_filt;




% ------- cluster --------
n_clust_desired = 4;
cluster_info.D0 = pdist(cluster_info.X0,'correlation');
cluster_info.Z0 = linkage(cluster_info.D0,'average');
perm = optimalleaforder(cluster_info.Z0, cluster_info.D0); % order
cluster_info.clust_id0 = cluster(cluster_info.Z0,'MaxClust',n_clust_desired); % cutoff
k = max(cluster_info.clust_id0); % number of clusters we got
% save workspace:
% clear trace_file freeze_file trace_data  all_traces  freeze_data  all_freeze;
%----------------

% save('hierarchical_clusters.mat', '-v7.3');

n_total = size(mega_traces,1);
cluster_full = nan(n_total,1);
cluster_full(good) = cluster_info.clust_id0;


cluster_info.cluster_full      = cluster_full;     % cluster ID per neuron (NaN = excluded)
cluster_info.good_mask         = good;             % logical mask used for clustering
cluster_info.mouse_id          = mouse_nums;       % original mouse index per neuron
cluster_info.day               = day_of_neuron;    % original day per neuron
cluster_info.neuron_index      = neuron_idx_within_mouse; % neuron ID within mouse
cluster_info.n_clusters        = k;
cluster_info.linkage_method    = 'average';
cluster_info.distance_metric   = 'correlation';
cluster_info.n_total_neurons   = n_total;

timestamp = datestr(now,'dd-mmm-yyyy');

cluster_info.note = [ ...
    'Cluster identities from hierarchical_v3_on ' timestamp '. ' ...
    'Clusters computed on filtered neurons (non-zero variance). ' ...
    'cluster_full is aligned to original mega_traces row order. ' ...
    'NaN entries indicate neurons excluded prior to clustering.'];

save('cluster_identities.mat','cluster_info');
