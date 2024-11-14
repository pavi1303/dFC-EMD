function [hs, e_ins, t_ins, f, t, mean_Ek, mean_Tk, log_mean_Ek, log_mean_Tk] = calculate_hht(imf, fs)
%% Estimation of mean energy and period profiles using Hilbert Huang Transform
% The IMF is of the form N*n_imf where
% N : number of timepoints
% n_imf : Number of IMF's extracted during EMD decomposition
[hs,f,t,f_ins, e_ins] = hht(imf, fs);
f_ins = abs(f_ins);
t_ins = abs(1./f_ins);
t_ins(t_ins>((1/fs)*size(imf,1))) = 60; % if the period is greater then the entire TR*tdim;
%f_ins = f_ins./fs; % Matlab calculates the frequency in Hz as fs/2*pi
mean_Ek = (1/size(imf,1))*sum(e_ins,1); % Mean energy
mean_Ek = mean_Ek./sum(mean_Ek); % Normalized mean energy
%Mean period using the density function
for k = 1:size(f_ins, 2)
    vk = f_ins(:, k)';
    [hvk, x1] = ksdensity(vk, 'support', [0, fs/2],'NumPoints', 1000);
    [cdf, x2] = ecdf(vk);
    cdflow = find(cdf > 0.001);
    cdfhigh = find(cdf > 1-0.001);
    hvklow = x2(cdflow(1));
    hvkhigh = x2(cdfhigh(1));
    ix1 = find(x1 > hvklow & x1 < hvkhigh);
    mean_Tk(:, k) = trapz(x1(ix1),hvk(ix1)./x1(ix1))/trapz(x1(ix1),hvk(ix1));
end
log_mean_Ek = log(mean_Ek); % log(energy)
log_mean_Tk = real(log(mean_Tk)); % log(period)
end