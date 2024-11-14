function imf_timeseries_plot(imf_data, grpname, plot_saveloc)
% imf_data : timepoints in the concatenated time series x number of imf's
% grpname : Name of the diagnostic group whose IMF is to be plotted
% plot_Saveloc : Location where the generated plot is to be saved
ncols = 3;
nrows = ceil(size(imf_data,2)/ncols);
colors = distinguishable_colors(size(imf_data,2));
figure('units','normalized','outerposition',[0 0 1 1]);
t = tiledlayout(nrows,ncols,TileSpacing="compact",Padding="compact");
for k = 1:size(imf_data,2)
    ax(k) = nexttile(t);
    plot(1:size(imf_data, 1), imf_data(:,k)', 'Color', colors(k, :));
    %xlim([tm(1) tm(end)])
    title(strcat("IMF ", num2str(k)));
    xlabel("Time frame", 'FontSize', 20);
    set(gca, 'FontSize', 18);
end
if ~exist(plot_saveloc,'dir')
    mkdir(plot_saveloc);
end
plot_title = strcat('Empirical mode decomposition (',grpname,')');
title(t,plot_title, 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(plot_saveloc,strcat('IMF_timeseries_', grpname, '.png')));
close(gcf);
end
