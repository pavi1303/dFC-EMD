function class_res = bfs_criterion_v2(X, Y, cat_ind, fea_names, t, hyperparam_variable_info, n_hyperparam_runs, best_t_input, NLearn_input, tuning_flag)
%% Classification using random forest and estimation of different performance metrics
if tuning_flag ==1
    class_model = fitcensemble(X, Y, 'Method', 'Bag', 'Learners', t, 'CategoricalPredictors', cat_ind, 'FResample',0.75,...
        'OptimizeHyperparameters', hyperparam_variable_info, 'HyperparameterOptimizationOptions', ...
        struct('ShowPlots',0, 'MaxObjectiveEvaluations', n_hyperparam_runs));
    best_t = class_model.ModelParameters.LearnerTemplates{1,1};
    %% Classifier model with the trained hyperparameters
    class_model_best = fitcensemble(X, Y, 'Method', 'Bag', 'Learners', best_t, 'CategoricalPredictors', cat_ind, 'FResample',0.75,...
        'NumLearningCycles',class_model.ModelParameters.NLearn);
else
    class_model_best = fitcensemble(X, Y, 'Method', 'Bag', 'Learners', best_t_input, 'CategoricalPredictors', cat_ind, 'FResample',0.75,...
        'NumLearningCycles',NLearn_input);
end
perm_imp = oobPermutedPredictorImportance(class_model_best);
[Y_pred, scores] = oobPredict(class_model_best);
diffscore = scores(:, 3) - max(scores(:, 1), scores(:, 2)); % One vs rest with AD vs the rest
%[fpr,tpr, ~, auc] = perfcurve(Y, scores(:,2), 1, 'NBoot', 100); %Bootstrap replicas using the max of the scores
[fpr,tpr, ~, auc] = perfcurve(Y, diffscore, 2, 'NBoot', 100); %Bootstrap replicas using the difference scores
c_mat = confusionmat(Y, Y_pred); % Confusion matrix
TP = c_mat(2, 2);
FN = sum(c_mat(2,:)) - c_mat(2,2);
FP = sum(c_mat(:, 2)) - c_mat(2,2);
% c_mat = confusionmat(Y, Y_pred); % Confusion matrix
% TP = c_mat(3, 3);
% FN = sum(c_mat(3,:)) - c_mat(3,3);
% FP = sum(c_mat(:, 3)) - c_mat(3,3);
Total_pred = sum(c_mat, 'all');
missclass_rate = (FP + FN)/Total_pred; % Misclassification rate

class_res = struct('X', {}, 'Y', {}, 'categorical_idx', {}, 'feature_names', {}, 'hyperparam_variable_info', {}, 'n_bayesopt_runs', {}, 'class_model', {}, 'class_model_best', {}, ...
    'perm_imp', {}, 'Y_pred', {}, 'scores', {}, 'missclass_rate', {}, 'AUC', {}, 'FPR', {}, 'TPR', {}, 'confusion_matrix', {});
%% Saving the classification results as a structure
class_res(1).X = X;
class_res(1).Y = Y;
class_res(1).categorical_idx = cat_ind;
class_res(1).feature_names = fea_names;
class_res(1).hyperparam_variable_info = hyperparam_variable_info;
class_res(1).n_bayesopt_runs = n_hyperparam_runs;
if tuning_flag ==1
    class_res(1).class_model = class_model;
else
    class_res(1).class_model = [];
end
class_res(1).class_model_best = class_model_best;
class_res(1).perm_imp = perm_imp;
class_res(1).Y_pred= Y_pred;
class_res(1).scores = scores;
class_res(1).missclass_rate = missclass_rate;
class_res(1).AUC = auc;
class_res(1).FPR = fpr;
class_res(1).TPR = tpr;
class_res(1).confusion_matrix = c_mat;
end