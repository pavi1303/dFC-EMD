function imf_freq_dist_plot(imf_data, TR, tdim, grpname, plot_saveloc)
fs = 1/TR;
f = fs*(0:(tdim/2))/tdim; %Frequency range limited to [0, fs/2] based on Nyquist criterion
[~,~,~,f_ins, ~] = hht(imf_data, fs);
ncols = size(f_ins, 2);
%colors = lines(ncols);
colors = distinguishable_colors(ncols);
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for i = 1:ncols
    [pdf, xi] = ksdensity(f_ins(:, i), f);
    plot(xi, pdf, 'Color', colors(i, :), 'LineWidth',2.5);
end
hxl = xline(0.01, '--', '0.01 Hz', 'LabelVerticalAlignment','middle', 'LineWidth',2.5);% Including a line for the drift frequency limit
hxl.FontSize = 18;
xticks(0:0.01:0.18);
legend(cellstr(num2str((1:ncols)', 'IMF %d')), 'Location','eastoutside');
box on;
xlabel('Frequency (Hz)', 'FontSize', 20);
set(gca, 'FontSize', 18);
grid on;
if ~exist(plot_saveloc,'dir')
    mkdir(plot_saveloc);
end
title(grpname, 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(plot_saveloc,strcat('Freq_dist_', grpname, '.png')));
close(gcf);
end
% for i = 1:ncols
%     [pdf_all(:, i), xi_all(:,i)] = ksdensity(f_ins(:, i), f);
%     %plot(xi, pdf, 'Color', colors(i, :), 'LineWidth',2.5);
% end