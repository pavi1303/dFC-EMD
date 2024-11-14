function run_kmeans_clustering_dfc(dfc_saveloc, window_type, nica, ncluster_min, ncluster_max)
temp_list = dir(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica))));
folder_list = temp_list([temp_list.isdir] & ~ismember({temp_list.name}, {'.','..'}));
folder_list = extractfield(folder_list, 'name');
n_reps = 100;
distance_metric = 'sqeuclidean';
max_iter = 500;
for j = 1:length(folder_list)
    dfc_resfile = extractfield(dir(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, j},'*.mat')), 'name');
    dfc_res = load(fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1, j}, dfc_resfile{1,1}));
    dfc_kmeans_input = dfc_res.dfc_kmeans_input;
    for n_clusters = ncluster_min:ncluster_max
        %% Performing standard k-means clustering
        [idx, clust_mat, sumd, D] = kmeans(dfc_kmeans_input', n_clusters,  'Replicates',n_reps, 'Distance', distance_metric, 'MaxIter', max_iter);
        cluster_saveloc = fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1,j}, strcat(num2str(n_clusters), '_clusters'));
        if ~exist(cluster_saveloc ,'dir')
            mkdir(cluster_saveloc);
        end
        save(fullfile(cluster_saveloc, sprintf(strcat('kmeans_cluster_res_', num2str(n_clusters), 'clusters.mat'))),'n_reps', ...
            'distance_metric', 'max_iter', 'n_clusters', 'idx', 'clust_mat', 'sumd', ...
            'D', 'dfc_kmeans_input');
        %% Deriving the brain states
        brain_states = cell(n_clusters, 1);
        for bs = 1:n_clusters
            bsstate_mat = eye(nica, nica);
            idx_l = tril(true(nica, nica), -1);
            idx_u = triu(true(nica, nica), 1);
            bsstate_mat(idx_l) = clust_mat(bs,:);
            brain_states{bs, 1} = tril(bsstate_mat)+tril(bsstate_mat,-1)';
        end
        bs_saveloc = fullfile(dfc_saveloc, window_type, strcat('GICA', num2str(nica)), folder_list{1,j}, ...
            strcat(num2str(n_clusters), '_clusters'), 'Brain_states');
        if ~exist(bs_saveloc ,'dir')
            mkdir(bs_saveloc);
        end
        save(fullfile(bs_saveloc, sprintf(strcat('Brain_states_', num2str(n_clusters), 'clusters.mat'))),'n_reps', ...
            'distance_metric', 'max_iter', 'n_clusters', 'brain_states');
    end
end
end