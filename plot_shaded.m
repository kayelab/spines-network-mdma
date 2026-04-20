function h = plot_shaded(t, mean_vec, err_vec, color_vec, varargin)
    % -------- defaults --------
    p = inputParser;
    addParameter(p, 'FaceAlpha', 0.25);
    addParameter(p, 'LineWidth', 2);
    addParameter(p, 'EdgeColor', 'none');
    addParameter(p, 'Marker', 'none');
    addParameter(p, 'MarkerSize', 6);
    addParameter(p, 'MarkerFaceColor', []);
    parse(p, varargin{:});
    opts = p.Results;

    % ensure row vectors
    t = t(:)';
    mean_vec = mean_vec(:)';
    err_vec = err_vec(:)';

    % build upper/lower bounds
    lower = mean_vec - err_vec;
    upper = mean_vec + err_vec;

    X = [t fliplr(t)];
    Y = [lower fliplr(upper)];

    h.patch = fill(X, Y, color_vec,'FaceAlpha', opts.FaceAlpha,'EdgeColor', opts.EdgeColor);
    hold on;

    if isempty(opts.MarkerFaceColor)
        mfc = color_vec;
    else
        mfc = opts.MarkerFaceColor;
    end

    h.line = plot(t, mean_vec,'Color', color_vec,'LineWidth', opts.LineWidth,'Marker', opts.Marker, 'MarkerSize', opts.MarkerSize,'MarkerFaceColor', mfc);
end
