function mi = simulate_mi(irun,para1,para2,w)

if irun==1
    p = para1;
    lambda = para2;
elseif irun==2
    rho = para1;
    beta = para2;
end

nsub = 1000;  % number of iteration
seq6 = repmat(1:6,1,6);
seq4 = repmat(1:4,1,9);
seq336 = repmat([1:3,1:3,1:6],1,3);
seq633 = repmat([1:6,1:3,1:3],1,3);
seq6no = 1:36;
seq_n1_1 = repmat([1,2,3,4,1,2,1,2,1,2,1,2],1,3);
seq_n1_2 = repmat([1,2,1,2,1,2,1,2,3,4,1,2],1,3);
seq_n2_1 = repmat([1,2,3,4,5,6,1,2,3,1,2,3,1,2,3,1,2,3],1,2);
seq_n2_2 = repmat([1,2,3,1,2,3,1,2,3,1,2,3,4,5,6,1,2,3],1,2);
subwholemi = nan(nsub,18+2+4);

if irun == 1
    % context resetting model    
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

        seq_n2_1_con = temreset(seq_n2_1,p,lambda);
        seq_n2_2_con = temreset(seq_n2_2,p,lambda);
        seq_n2_1_submi = mical(seq_n2_1_con,'seq123456');
        seq_n2_2_submi = mical(seq_n2_2_con,'seq123123');
        subwholemi(isub,21:24) = (seq_n2_1_submi+seq_n2_2_submi)/2;
    end
elseif irun == 2
    % position coding model, drift-like codes, different weight ratio
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

        seq_n2_1_con = temcon(seq_n2_1,rho,100-w);
        seq_n2_2_con = temcon(seq_n2_2,rho,100-w);
        seq_n2_1_l = postem(seq_n2_1,beta,w);
        seq_n2_2_l = postem(seq_n2_2,beta,w);
        seq_n2_1_submi = mical([seq_n2_1_con,seq_n2_1_l],'seq123456');
        seq_n2_2_submi = mical([seq_n2_2_con,seq_n2_2_l],'seq123123');
        subwholemi(isub,21:24) = (seq_n2_1_submi+seq_n2_2_submi)/2;
    end  
end
    mi = mean(subwholemi,1)';
    clear subwholemi
end