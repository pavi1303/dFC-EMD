function [corrvec_all, td_in_s_all, ws_emd_used_all, dfc_kmeans_input] = dfc_emd(window_type, rsn_names_ordered, n_sub_grp, rsn_trs_grp_emd_sigrecon, imfs_to_choose)
%% Adaptive window size estimation based on EMD
%% Inputs
% window_type: str of 'Rectangular' or 'Nonrectangular'
% rsn_names_ordered: Cell array of RSN names used for dFC estimation
% rsn_trs_grp_emd_sigrecon: Signal reconstructed from chosen IMFs for each RSN
% imfs_to_choose: IMFs chosen for calculation of adaptive window
%% Outputs
% corrvec_all: Correlation vectors across windows for each group
% td_in_s_all: Adaptive window size across windows for each group
% ws_emd_used_all: Mean adaptive window size across groups
% dfc_kmeans_input: Input for k-means clustering for extraction of brain states
n_rsn_comp = (length(rsn_names_ordered)*(length(rsn_names_ordered)-1))/2; % Number of correlation comparisons
td_in_s_all = cell(length(n_sub_grp), 1);
corrvec_all = cell(length(n_sub_grp), 1);
ws_emd_used_all = cell(1, length(n_sub_grp));
for grp = 1:length(n_sub_grp)
    n_patients = n_sub_grp(grp);
    rsn_trs_grp_emd = rsn_trs_grp_emd_sigrecon{1, grp};
    %% Obtaining the imfs for each subject for each rsn
    imf_sub_data = cell(n_patients, length(rsn_names_ordered));
    for rsn = 1:length(rsn_names_ordered)
        imf_loc = fullfile(emd_results_loc, rsn_names_ordered{1, rsn}, 'IMF', dx_status{grp, 1});
        imf_filename = extractfield(dir(fullfile(imf_loc, '*.mat')), 'name')';
        imf = load(fullfile(imf_loc, imf_filename{1})).imf;
        imf = imf(imfs_to_choose, :)';
        imf_sub_idx = mat2cell(1:size(imf, 1), 1, repmat(tdim, 1, n_patients));
        imf_sub_data(:, rsn)= cellfun(@(rows) imf(rows, :), imf_sub_idx, 'UniformOutput', false)';
    end
    %% Calculation of the time dependent window size using instantaneous period
    td_in_s = cell(n_patients, 1);
    corrmat_grp = cell(n_patients, 1);
    corrvec_grp = cell(n_patients, 1);
    ws_emd_used_grp = zeros(n_patients, 2);
    for sub = 1:n_patients
        imf_oi= imf_sub_data(sub, :);
        ins_T = cell(n_patients, length(rsn_names_ordered)); % Cell array with the instantaneous
        mean_Ek = cell(n_patients, length(rsn_names_ordered));
        td_in_s_rsn = zeros(tdim, length(imf_oi));
        for i = 1:length(imf_oi)
            ins_T_oi = zeros(tdim, size(imf_oi{1, i}, 2)); mean_Ek_oi = zeros(1, size(imf_oi{1, i}, 2));
            [~, ~, ins_T_oi, ~, ~, mean_Ek_oi, ~, ~, ~] = calculate_hht(imf_oi{1, i}, fs);
            ins_T_oi = ins_T_oi(:, :);mean_Ek_oi = mean_Ek_oi(:, :);
            ins_T{sub, i} = ins_T_oi;mean_Ek{sub, i} = mean_Ek_oi;
            td_in_s_oi = sum(ins_T_oi.*mean_Ek_oi, 2);
            td_in_s_oi(td_in_s_oi>tdim) = floor(tdim);
            td_in_s_rsn(:, i) = td_in_s_oi;
        end
        td_in_s{sub, 1} = td_in_s_rsn;
        %% Calculate the window size using the td (td_in_s_rsn)
        window_size_rsn = cell(size(td_in_s_rsn,2), size(td_in_s_rsn,2));
        for rsn1 = 1:size(td_in_s_rsn,2)
            td_in_s_rsn1 = td_in_s_rsn(:, rsn1);
            for rsn2 = 1:rsn1
                td_in_s_rsn2 = td_in_s_rsn(:, rsn2);
                td_in_s_temp = [td_in_s_rsn1'; td_in_s_rsn2'];
                td = (max(td_in_s_temp)./tr)';
                window_size_rsn{rsn1, rsn2} = td;
            end
        end
        %% Calculation of the sliding window correlation coefficient using the window size above
        rsn_trs_sub = rsn_trs_grp_emd{sub, 1};
        corrmat_t = NaN(length(rsn_names_ordered), length(rsn_names_ordered));
        corrmat_t_concat = cell(1, tdim-1);
        corrvec_concat= NaN(n_rsn_comp, tdim-1);
        for t = 1:tdim-1
            for rsnid1 = 1:length(rsn_names_ordered)
                for rsnid2 = 1:rsnid1
                    window_sizeoi = window_size_rsn{rsnid1, rsnid2};
                    x1 = rsn_trs_sub(:, rsnid1);x2 = rsn_trs_sub(:, rsnid2);
                    if round(t-window_sizeoi(t,1)/2) >= 1 && round(t+window_sizeoi(t,1)/2) <= tdim
                        if strcmp(window_type, 'Rectangular')
                            corrmat_t(rsnid1,rsnid2) = corr(x1(round(t-window_sizeoi(t,1)/2):round(t+window_sizeoi(t,1)/2)), ...
                                x2(round(t-window_sizeoi(t,1)/2):round(t+window_sizeoi(t,1)/2)));
                            for i  = 1:length(rsn_names_ordered)
                                corrmat_t(i, i) = 1; % Making the diagonal elements 1
                            end
                            corrvec_t = tril(corrmat_t, -1);
                            corrvec_concat(:, t) = corrvec_t(corrvec_t~=0);
                        else
                            x1_window = x1(round(t-window_sizeoi(t,1)/2):round(t+window_sizeoi(t,1)/2));
                            x2_window = x2(round(t-window_sizeoi(t,1)/2):round(t+window_sizeoi(t,1)/2));
                            Y = horzcat(x1_window, x2_window);
                            [dt, N] = size(Y);
                            %% Non rectangular window selection
                            nonrectangular_wt = ACgausswin(dt, 3); % Approximate confined Gaussian
                            %nonrectangular_wt = hann(dt); % Hanning window
                            %nonrectangular_wt = hamming(dt); % Hamming window
                            %nonrectangular_wt = parzenwin(dt); % Parzen window
                            %nonrectangular_wt = gausswin(dt); % Gaussian window
                            temp = Y - repmat(nonrectangular_wt' * Y, dt, 1); % Remove weighted mean
                            temp = temp' * (temp .* repmat(nonrectangular_wt, 1, N)); % Weighted covariance matrix
                            temp = 0.5 * (temp + temp'); % Making it symmetric
                            R = diag(temp); % Variances
                            R = temp ./ sqrt(R * R'); % Weighted Pearson correlation matrix
                            corrmat_t(rsnid1,rsnid2) = R(1,2);
                            for i  = 1:length(rsn_names_ordered)
                                corrmat_t(i, i) = 1; % Making the diagonal elements 1
                            end
                            mask = tril(true(size(corrmat_t)),-1);
                            corrvec_concat(:, t) = corrmat_t(mask);
                        end
                    end
                    corrmat_t_concat{1, t} = corrmat_t;
                end
            end
        end
        corrmat_grp{sub, :} = corrmat_t_concat;
        corrvec_grp{sub, :} = corrvec_concat;
        mask = tril(true(size(corrmat_t)),-1);
        ws_emd_used_grp(sub, :) = [mean(cell2mat(window_size_rsn(mask)), "all"), std(cell2mat(window_size_rsn(mask)),0,"all")];
    end
    % Temporal concatenation of correlation matrices across all subjects within a group
    corrvec_all{grp, 1} = cell2mat(corrvec_grp');
    td_in_s_all{grp, 1} = td_in_s;
    % Finding the mean td (s) in group
    ws_emd_used_all{1, grp} = ws_emd_used_grp;
end
dfc_kmeans_input = cell2mat(corrvec_all'); % Temporal concatenation across groups
end