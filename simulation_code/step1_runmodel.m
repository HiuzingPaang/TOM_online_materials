function runmi = step1_runmodel(irun)

%% for repeatability and set some parameters
rng('default');
nsub = 200;  % number of iteration
seq6 = repmat(1:6,1,6);
seq4 = repmat(1:4,1,9);
seq336 = repmat([1:3,1:3,1:6],1,3);
seq633 = repmat([1:6,1:3,1:3],1,3);
seq6no = 1:36;
seq_n1_1 = repmat([1,2,3,4,1,2,1,2,1,2,1,2],1,3);
seq_n1_2 = repmat([1,2,1,2,1,2,1,2,3,4,1,2],1,3);
    
%% run model
ind = 0;
runmi = nan(99*99,3+18+2,nsub);
if irun == 1
    %% context resetting model
    for p = 0.01:0.01:0.99
        tic
        for lambda = 0.01:0.01:0.99
            ind = ind + 1;
            subwholemi = nan(nsub,18+2);
            for isub = 1:nsub
                seq6_con = temreset(seq6,p,lambda);
                seq4_con = temreset(seq4,p,lambda);
                seq336_con = temreset(seq336,p,lambda);
                seq633_con = temreset(seq633,p,lambda);
                seq6no_con = temreset(seq6no,p,lambda);

                subwholemi(isub,1:5) = mical(seq6_con,'seq6');
                subwholemi(isub,6:11) = mical(seq4_con,'seq4');
                seq336_submi = mical(seq336_con,'seq336');
                seq633_submi = mical(seq633_con,'seq633');
                subwholemi(isub,12:16) = (seq336_submi+seq633_submi)/2;
                subwholemi(isub,17:18) = mical(seq6no_con,'seq6no');

                seq_n1_1_con = temreset(seq_n1_1,p,lambda);
                seq_n1_2_con = temreset(seq_n1_2,p,lambda);
                seq_n1_1_submi = mical(seq_n1_1_con,'seq1234');
                seq_n1_2_submi = mical(seq_n1_2_con,'seq1212');
                subwholemi(isub,19:20) = (seq_n1_1_submi+seq_n1_2_submi)/2;
            end
            runmi(ind,1,:) = irun;
            runmi(ind,2,:) = p;
            runmi(ind,3,:) = lambda;
            runmi(ind,4:end,:) = subwholemi';
        end
        toc
    end
elseif irun == 2
    %% position coding model, drift-like codes, 50% weight
    w = 50;
    for rho = 0.01:0.01:0.99
        tic
        for beta = 0.01:0.01:0.99
            ind = ind + 1;
            subwholemi = nan(nsub,18+2);
            for isub = 1:nsub
                seq6_con = temcon(seq6,rho,100-w);
                seq4_con = temcon(seq4,rho,100-w);
                seq336_con = temcon(seq336,rho,100-w);
                seq633_con = temcon(seq633,rho,100-w);
                seq6no_con = temcon(seq6no,rho,100-w);

                seq6_l = postem(seq6,beta,w);
                seq4_l = postem(seq4,beta,w);
                seq336_l = postem(seq336,beta,w);
                seq633_l = postem(seq633,beta,w);
                seq6no_l = postem(seq6no,beta,w);

                subwholemi(isub,1:5) = mical([seq6_con,seq6_l],'seq6');
                subwholemi(isub,6:11) = mical([seq4_con,seq4_l],'seq4');
                seq336_submi = mical([seq336_con,seq336_l],'seq336');
                seq633_submi = mical([seq633_con,seq633_l],'seq633');
                subwholemi(isub,12:16) = (seq336_submi+seq633_submi)/2;
                subwholemi(isub,17:18) = mical([seq6no_con,seq6no_l],'seq6no');

                seq_n1_1_con = temcon(seq_n1_1,rho,100-w);
                seq_n1_2_con = temcon(seq_n1_2,rho,100-w);
                seq_n1_1_l = postem(seq_n1_1,beta,w);
                seq_n1_2_l = postem(seq_n1_2,beta,w);
                seq_n1_1_submi = mical([seq_n1_1_con,seq_n1_1_l],'seq1234');
                seq_n1_2_submi = mical([seq_n1_2_con,seq_n1_2_l],'seq1212');
                subwholemi(isub,19:20) = (seq_n1_1_submi+seq_n1_2_submi)/2;
            end
            runmi(ind,1,:) = irun;
            runmi(ind,2,:) = rho;
            runmi(ind,3,:) = beta;
            runmi(ind,4:end,:) = subwholemi';
        end
        toc
    end
end

save(['/seastor/a853898293/analysis/data_mi_run',num2str(irun),'.mat'],'runmi');

end