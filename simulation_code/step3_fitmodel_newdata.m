function step3_fitmodel_newdata(irun,ipart)

%% load data
load(['/seastor/a853898293/analysis/data_mi_run',num2str(irun),'.mat'],'runmi');

%% fit the mi with Pu's real data and the new one
ncolumn = 4:23;
all_bhv = [77.65, 67.39, 75.13, 70.52, 59.85, ...
           79.04, 66.53, 56.48, 74.99, 74.98, 69.48, ...
           78.05, 74.88, 73.34, 68.08, 75.10, ...
           72.89, 74.45, ...
           69.47, 69.54]';
subjID=[ones(2,1); ones(3,1)*2; ones(1,1)*3; ones(2,1)*4; ones(3,1)*5; ones(5,1)*6; ones(2,1)*7; ...
    ones(2,1)*8];
subjID=categorical(subjID);

nrep = size(runmi,3);
npart = 9;
nlong = size(runmi,1)/npart;

i0 = (ipart-1)*nlong + 1;
i1 = ipart*nlong;
runmi_ = runmi(i0:i1,:,:);
clear runmi

crit_glm = nan(nlong,nrep,6);  % ordinary adjusted AIC BIC p f
crit_glmm = nan(nlong,nrep,6);

for ilog = 1:nlong
    tic
    for irep = 1:nrep
        A = mean(runmi_(ilog,ncolumn,1:irep),3);
        modelA = A(:);
        tbl1 = table(modelA,subjID,all_bhv);

        % fixed effect model
        mdl1_fixed = fitglm(tbl1,'all_bhv ~ 1 + modelA');
        % extract values
        crit_glm(ilog,irep,1) = mdl1_fixed.Rsquared.Ordinary;
        crit_glm(ilog,irep,2) = mdl1_fixed.Rsquared.Adjusted;
        crit_glm(ilog,irep,3) = mdl1_fixed.ModelCriterion.AIC;
        crit_glm(ilog,irep,4) = mdl1_fixed.ModelCriterion.BIC;
        [p,f] = coefTest(mdl1_fixed);  % compute p and f
        crit_glm(ilog,irep,5) = p;
        crit_glm(ilog,irep,6) = f;

        % random effect model full_model
        mdl1_rdm = fitglme(tbl1,'all_bhv ~ 1 + modelA + (1+modelA|subjID)');
        % extract values
        crit_glmm(ilog,irep,1) = mdl1_rdm.Rsquared.Ordinary;
        crit_glmm(ilog,irep,2) = mdl1_rdm.Rsquared.Adjusted;
        crit_glmm(ilog,irep,3) = mdl1_rdm.ModelCriterion.AIC;
        crit_glmm(ilog,irep,4) = mdl1_rdm.ModelCriterion.BIC;
        [p,f] = coefTest(mdl1_rdm);  % compute p and f
        crit_glmm(ilog,irep,5) = p;
        crit_glmm(ilog,irep,6) = f;
    end
    toc
end

save(['/seastor/a853898293/analysis/data_FitNewData_run',num2str(irun),'_part',num2str(ipart),'.mat'],...
    'crit_glm','crit_glmm');

end
