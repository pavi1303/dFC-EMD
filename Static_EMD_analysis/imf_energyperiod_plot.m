function imf_energyperiod_plot(Ek, Tk, f_drift, grpname, plot_saveloc)
% Ek - Mean energy for each IMF in the log scale (y - axis)
% Tk - Mean period for each IMF in the log scale (x - axis)

colors = distinguishable_colors(size(Ek, 2));
figure('Renderer', 'painters', 'Position', [1000 1000 2000 2000]);
hold on;
for k = 1:size(Ek, 2)
    scatter(Tk(:, k),Ek(:, k),20,colors(k,:), 'filled');
end
hxl = xline(log(1/f_drift), '--', '0.01 Hz', 'LabelVerticalAlignment','middle', 'LineWidth',2.5);% Including a line for the drift frequency limit
hxl.FontSize = 18;
%legend(cellstr(num2str(1:size(hht_Ek, 2), 'IMF %d')));
legend(num2str((1:size(Ek, 2))', 'IMF %d'), 'FontSize',30, 'Location','eastoutside');
box on;
hold off;
xlabel('log(period)', 'FontSize', 20);
ylabel('log(energy)', 'FontSize', 20);
set(gca, 'FontSize', 20);
if ~exist(plot_saveloc,'dir')
    mkdir(plot_saveloc);
end
title(grpname, 'FontSize', 24, 'FontWeight', 'bold');
saveas(gcf, fullfile(plot_saveloc,strcat('Energy_period_', grpname, '.png')));
close(gcf);
end