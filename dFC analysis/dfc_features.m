function dfc_features(dfc_saveloc, window_type, nica, tdim, dx_status, n_sub_grp)
%%  Estimation of the dFC features
temp_list = dir(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica))));
folder_list = temp_list([temp_list.isdir] & ~ismember({temp_list.name}, {'.','..'}));
folder_list = extractfield(folder_list, 'name');
for ws_opt = 1:length(folder_list)
    cluster_list = dir(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, ws_opt}));
    cluster_list = cluster_list([cluster_list.isdir] & ~ismember({cluster_list.name}, {'.','..'}));
    cluster_list = extractfield(cluster_list, 'name');
    for cluster_opt = 1:length(cluster_list)
        cluster_resfile = extractfield(dir(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), ...
            folder_list{1, ws_opt}, cluster_list{1, cluster_opt},'*.mat')), 'name');
        cluster_res = load(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), ...
            folder_list{1, ws_opt}, cluster_list{1, cluster_opt}, cluster_resfile{1,1}));
        n_clusters = cluster_res.n_clusters;
        state_trans_vec = cluster_res.idx;
        state_trans_vec_grp = cell(size(dx_status, 1), 1);
        state_trans_vec_grp{1,1} = state_trans_vec(1:n_sub_grp(1,1)*(tdim-1));
        state_trans_vec_grp{2,1} = state_trans_vec(n_sub_grp(1,1)*(tdim-1)+1: ...
            n_sub_grp(1,1)*(tdim-1)+n_sub_grp(1,2)*(tdim-1));
        state_trans_vec_grp{3,1} = state_trans_vec(n_sub_grp(1,1)*(tdim-1)+n_sub_grp(1,2)*(tdim-1)+1:end);
        fraction_time_all = cell(size(dx_status, 1), 1);
        freq_state_change_all = cell(size(dx_status, 1), 1);
        mdt_all = cell(size(dx_status, 1), 1);
        fraction_time_mean = zeros(size(n_sub_grp, 2), n_clusters*2);
        freq_state_change_mean = zeros(size(n_sub_grp, 2), 2);
        mdt_mean = zeros(size(n_sub_grp, 2), n_clusters*2);
        state_trans_grp_all = cell(size(n_sub_grp, 2), 1);
        state_trans_grp_mean = cell(size(n_sub_grp, 2), 1);
        for grp_idx = 1:size(state_trans_vec_grp, 1)
            state_trans_vec_sub = mat2cell(state_trans_vec_grp{grp_idx,1}, diff([0:(tdim-1): ...
                numel(state_trans_vec_grp{grp_idx,1})-1,numel(state_trans_vec_grp{grp_idx,1})]));
            count_fraction_grp = zeros(n_sub_grp(1, grp_idx), n_clusters);
            freq_state_change_grp = zeros(n_sub_grp(1, grp_idx), 1);
            state_trans_grp_norm = zeros(n_clusters, n_clusters, n_sub_grp(1, grp_idx));
            mdt_grp = zeros(n_sub_grp(1, grp_idx), n_clusters);
            for sub = 1:size(state_trans_vec_sub, 1)
                sub_vec = state_trans_vec_sub{sub, 1};
                sub_vec=sub_vec(~isnan(sub_vec));
                count = histcounts(sub_vec, [1:n_clusters, inf]);
                %% 1. Fraction time
                count_fraction_grp(sub,:) = count./sum(count, "all");
                %count_fraction_grp(sub,:) = count./(tdim-1);
                %% 2. State change frequency
                freq_state_change_grp(sub, 1) = nnz(diff(sub_vec))/length(diff(sub_vec));
                state_trans_sub = zeros(n_clusters, n_clusters);
                for i = 1:length(sub_vec)-1
                    current_state = sub_vec(i);
                    next_state = sub_vec(i+1);
                    state_trans_sub(current_state, next_state) = state_trans_sub(current_state, next_state)+ 1;
                end
                %% 3. Transition probabilities
                state_trans_sub_norm = state_trans_sub./sum(state_trans_sub, 2);
                state_trans_grp_norm(:, :, sub) = state_trans_sub_norm;
                %% 4. Mean dwell time
                mdt_grp(sub, :) = calculate_mdt(sub_vec, n_clusters);
            end
            fraction_time_all{grp_idx, 1} = count_fraction_grp;
            freq_state_change_all{grp_idx, 1} = freq_state_change_grp;
            mdt_all{grp_idx, 1} = mdt_grp;
            state_trans_grp_all{grp_idx, 1} = state_trans_grp_norm;
            state_trans_grp_mean{grp_idx, 1} = nanmean(state_trans_grp_norm, 3);
            fraction_time_mean(grp_idx, :) = [mean(count_fraction_grp, 1), std(count_fraction_grp, 1)];
            freq_state_change_mean(grp_idx, :) = [mean(freq_state_change_grp); std(freq_state_change_grp)];
            mdt_mean(grp_idx, :) = [mean(mdt_grp, 1), std(mdt_grp, 1)];
        end
        %% Saving the results
        %% 1. Fraction time
        ft_results_saveloc = fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, ws_opt}, cluster_list{1, cluster_opt}, 'Fraction_time');
        if ~exist(ft_results_saveloc ,'dir')
            mkdir(ft_results_saveloc);
        end
        save(fullfile(ft_results_saveloc, ...
            sprintf(strcat('FT_results_', num2str(n_clusters), 'clusters.mat'))),'fraction_time_all', 'fraction_time_mean', ...
            'state_trans_vec_grp');
        %% 2. State change frequency
        scf_results_saveloc = fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, ws_opt}, cluster_list{1, cluster_opt}, 'State_change_frequency');
        if ~exist(scf_results_saveloc ,'dir')
            mkdir(scf_results_saveloc);
        end
        save(fullfile(scf_results_saveloc, ...
            sprintf(strcat('SCF_results_', num2str(n_clusters), 'clusters.mat'))),'freq_state_change_all', 'freq_state_change_mean', 'state_trans_vec_grp');
        %% 3. Probability of state transitions
        tp_results_saveloc = fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, ws_opt}, cluster_list{1, cluster_opt}, 'Transition_probabilities');
        if ~exist(tp_results_saveloc ,'dir')
            mkdir(tp_results_saveloc);
        end
        save(fullfile(tp_results_saveloc, ...
            sprintf(strcat('TP_results_', num2str(n_clusters), 'clusters.mat'))),'state_trans_grp_all', 'state_trans_grp_mean','state_trans_vec_grp');
        %% 4. Mean dwell time
        mdt_results_saveloc = fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, ws_opt}, cluster_list{1, cluster_opt}, 'Mean_dwell_time');
        if ~exist(mdt_results_saveloc ,'dir')
            mkdir(mdt_results_saveloc);
        end
        save(fullfile(mdt_results_saveloc, ...
            sprintf(strcat('MDT_results_', num2str(n_clusters), 'clusters.mat'))),'mdt_all', 'mdt_mean', 'state_trans_vec_grp');
    end
end
end

