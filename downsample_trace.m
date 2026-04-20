%% downsample_trace.m
% 3/4/26
% helper function to reduce traces / downsample

function ds_trace = downsample_trace(x, smooth_win_frames, bin_size_frames)
% smooth/movmean & bin a single trace (1xtime -> 1xnBins) or trial matrix (trialxtime -> trials x nBins)

    if isempty(x) % check
        ds_trace = [];
        return;
    end

    smooth_win_frames = max(1, round(smooth_win_frames));
    bin_size_frames = max(1, round(bin_size_frames));

    is_vector = isvector(x); % is input a vec or matrx
    if is_vector
        ds_trace = reduceVector_(x, smooth_win_frames, bin_size_frames);
    else
        % Expect [trials x time]
        ds_trace = reduceTrials_(x, smooth_win_frames, bin_size_frames);
    end
end

% -- hopefully will work - local subfunctions? --

function out = reduceTrials_(trials_by_time, smooth_win, bin)
    tr = movmean(trials_by_time, smooth_win, 2); % smooth using movmean
    n_bins = floor(size(tr,2) / bin);  % bin data
    if n_bins < 2
        out = [];
        return;
    end
    tr = tr(:, 1:n_bins*bin);
    tr = reshape(tr, size(tr,1), bin, n_bins);
    out = squeeze(mean(tr, 2));
end

function out = reduceVector_(x, smooth_win, bin)
    x = x(:)'; % force row
    x = movmean(x, smooth_win); % smoothie movmean
    n_bins = floor(numel(x) / bin); % bin data
    if n_bins < 2
        out = [];
        return;
    end
    x = x(1:n_bins*bin);
    x = reshape(x, bin, n_bins);
    out = mean(x, 1);
end