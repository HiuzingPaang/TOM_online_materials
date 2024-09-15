function fitmodel_realdata(irun)

%% run data
fitmodel_mi_prepare(irun);

%% load data
load(['/seastor/a853898293/analysis/model_simulation/fitmodel_mi_prepa_run_',num2str(irun),'.mat'],'runmi');
runmi = mean(runmi,3);

%% fit the mi with all real data
ncolumn = 4:27;
all_bhv = [77.65, 67.39, 75.13, 70.52, 59.85, ...
           79.04, 66.53, 56.48, 74.99, 74.98, 69.48, ...
           78.05, 74.88, 73.34, 68.08, 75.10, ...
           72.89, 74.45, ...
           69.47, 69.54, ...
           67.75, 69.64, 63.84, 59.71]';
subjID=[ones(2,1); ones(3,1)*2; ones(1,1)*3; ones(2,1)*4; ones(3,1)*5; ones(5,1)*6; ones(2,1)*7; ...
    ones(2,1)*8; ones(2,1)*9; ones(2,1)*10];
subjID=categorical(subjID);

crit_glm = nan(size(runmi,1),6);  % ordinary adjusted AIC BIC p f
crit_glmm = nan(size(runmi,1),6);

for ilog = 1:size(runmi,1)
    
    A = runmi(ilog,ncolumn);
    modelA = A(:);
    tbl1 = table(modelA,subjID,all_bhv);

    % fixed effect model
    mdl1_fixed = fitglm(tbl1,'all_bhv ~ 1 + modelA');
    % extract values
    crit_glm(ilog,1) = mdl1_fixed.Rsquared.Ordinary;
    crit_glm(ilog,2) = mdl1_fixed.Rsquared.Adjusted;
    crit_glm(ilog,3) = mdl1_fixed.ModelCriterion.AIC;
    crit_glm(ilog,4) = mdl1_fixed.ModelCriterion.BIC;
    [p,f] = coefTest(mdl1_fixed);  % compute p and f
    crit_glm(ilog,5) = p;
    crit_glm(ilog,6) = f;

    % random effect model full_model
    mdl1_rdm = fitglme(tbl1,'all_bhv ~ 1 + modelA + (1+modelA|subjID)');
    % extract values
    crit_glmm(ilog,1) = mdl1_rdm.Rsquared.Ordinary;
    crit_glmm(ilog,2) = mdl1_rdm.Rsquared.Adjusted;
    crit_glmm(ilog,3) = mdl1_rdm.ModelCriterion.AIC;
    crit_glmm(ilog,4) = mdl1_rdm.ModelCriterion.BIC;
    [p,f] = coefTest(mdl1_rdm);  % compute p and f
    crit_glmm(ilog,5) = p;
    crit_glmm(ilog,6) = f;
    
    save(['/seastor/a853898293/analysis/model_simulation/fitmodel_prepa_realdata_run_',num2str(irun),'.mat'],...
        'crit_glm','crit_glmm');
end

end
