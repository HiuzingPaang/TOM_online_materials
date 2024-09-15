function recovery_model_run1(irun_simulate)

irun = 1;

ncolumn = 4:27;
para_loop1 = 0.01:0.06:0.99;

subjID = [ones(2,1); ones(3,1)*2; ones(1,1)*3; ones(2,1)*4; ones(3,1)*5; ones(5,1)*6; ones(2,1)*7;...
    ones(2,1)*8; ones(2,1)*9; ones(2,1)*10];
subjID = categorical(subjID);

crit_glm_1 = nan(2500,length(para_loop1),length(para_loop1));  % ordinary R^2
crit_glm_2 = nan(2500,length(para_loop1),length(para_loop1));  % AIC
crit_glm_3 = nan(2500,length(para_loop1),length(para_loop1));  % BIC
crit_glmm_1 = nan(2500,length(para_loop1),length(para_loop1));  % ordinary R^2
crit_glmm_2 = nan(2500,length(para_loop1),length(para_loop1));  % AIC
crit_glmm_3 = nan(2500,length(para_loop1),length(para_loop1));  % BIC
load(['/seastor/a853898293/analysis/model_simulation/data_mi_run_',num2str(irun),'.mat'],'runmi');
for i1 = 1:length(para_loop1)
    for i2 = 1:length(para_loop1)
        tic
        
        % create the simulated mi
        rng(i1+i2+irun_simulate);
        if irun_simulate==1
            all_bhv = simulate_mi(irun_simulate,para_loop1(i1),para_loop1(i2));
        elseif irun_simulate==2
            all_bhv = simulate_mi(irun_simulate,para_loop1(i1),para_loop1(i2),50);
        end
        for ilog = 1:size(runmi,1)
            A = runmi(ilog,ncolumn);
            modelA = A(:);
            tbl1 = table(modelA,subjID,all_bhv);

            % fixed effect model
            mdl1_fixed = fitglm(tbl1,'all_bhv ~ 1 + modelA');
            % extract values
            crit_glm_1(ilog,i1,i2) = mdl1_fixed.Rsquared.Ordinary;
            crit_glm_2(ilog,i1,i2) = mdl1_fixed.ModelCriterion.AIC;
            crit_glm_3(ilog,i1,i2) = mdl1_fixed.ModelCriterion.BIC;

            % random effect model full_model
            mdl1_rdm = fitglme(tbl1,'all_bhv ~ 1 + modelA + (1+modelA|subjID)');
            % extract values
            crit_glmm_1(ilog,i1,i2) = mdl1_rdm.Rsquared.Ordinary;
            crit_glmm_2(ilog,i1,i2) = mdl1_rdm.ModelCriterion.AIC;
            crit_glmm_3(ilog,i1,i2) = mdl1_rdm.ModelCriterion.BIC;
        end
       
        toc
        save(['/seastor/a853898293/analysis/model_simulation/para_rec_run_',num2str(irun),'_isim_',num2str(irun_simulate),'_r.mat'],...
            'crit_glm_1','crit_glmm_1');
        save(['/seastor/a853898293/analysis/model_simulation/para_rec_run_',num2str(irun),'_isim_',num2str(irun_simulate),'_aic.mat'],...
            'crit_glm_2','crit_glmm_2');
        save(['/seastor/a853898293/analysis/model_simulation/para_rec_run_',num2str(irun),'_isim_',num2str(irun_simulate),'_bic.mat'],...
            'crit_glm_3','crit_glmm_3');
    end
end

end