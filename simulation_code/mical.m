function submi = mical(wholevector,eventype)
    
    % calculate the similarity between all event vectors
    corrs = squareform(1-pdist(wholevector,'cosine'));
    
    % calculate the mi
    if strcmp(eventype,'seq6')
        
        % exp1_6_within, exp1_6_across, exp3_6_within(early), exp3_6_within(late), exp3_6_across
        stim1 = [2,8,14,20,26,33,  5,11,17,23,29,  2,14,26, 10,22,34, 5,11,23,29];
        stim2 = [6,12,18,24,30,36  9,15,21,27,33,  4,16,28, 12,24,36, 9,15,27,33];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:6)),mean(diff_stim(7:11)),mean(diff_stim(12:14)),...
            mean(diff_stim(15:17)),mean(diff_stim(18:21))];
        
    elseif strcmp(eventype,'seq4')
        
        % exp2_4_withinlag1, exp???_4_acrosslag3, exp???_4_acrosslag1,
        % exp3_4_withinlag1(red), exp3_4_withinlag1(green), exp3_4_acrosslag3
        stim1 =[2,6,10,14,18,22,26,30,34,  3,7,11,15,19,23,27,31,  3,7,11,15,19,23,27,31,  2,14,26, 10,22,34, 3,11,23,31];
        stim2 =[4,8,12,16,20,24,28,32,36, 7,11,15,19,23,27,31,35,  5,9,13,17,21,25,29,33,  4,16,28, 12,24,36, 7,15,27,35];

        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:9)),mean(diff_stim(10:17)),mean(diff_stim(18:25)),...
            mean(diff_stim(26:28)),mean(diff_stim(29:31)),mean(diff_stim(32:35))];
        
    elseif strcmp(eventype,'seq336')
        
        % exp4_336_withinlag1(early), exp4_336_withinlag1(late),
        % exp4_336_acrosslag3(early), exp4_336_acrosslag3(late),
        % exp4_336_grey
        stim1= [4,16,28,  10,22,34,   5,17,29,  11,23,  2,14,26,7,19];
        stim2= [6,18,30,  12,24,35,   9,21,33,  15,27,  8,20,32,13,25];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:3)),mean(diff_stim(4:6)),mean(diff_stim(7:9)),...
            mean(diff_stim(10:11)),mean(diff_stim(12:16))];
        
    elseif strcmp(eventype,'seq633')
        
        % exp4_633_withinlag1(early), exp4_633_withinlag1(late),
        % exp4_633_acrosslag3(early), exp4_633_acrosslag3(late),
        % exp4_633_grey
        stim1= [10,22,34,  4,16,28,   11,23,  5,17,29,  2,14,26,7,19];
        stim2= [12,24,35,  6,18,30,   15,27,  9,21,33,  8,20,32,13,25];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:3)),mean(diff_stim(4:6)),mean(diff_stim(7:8)),...
            mean(diff_stim(9:11)),mean(diff_stim(12:16))];
        
    elseif strcmp(eventype,'seq6no')
        stim1 = [2,8,14,20,26,33,  5,11,17,23,29];
        stim2 = [6,12,18,24,30,36  9,15,21,27,33];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:6)),mean(diff_stim(7:11))];
        
    elseif strcmp(eventype,'seq1234')
        % across one, across two
        stim1 = [2,14,26,   8,20,32];
        stim2 = [6,18,30,  12,24,36];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:3)),mean(diff_stim(4:6))];
        
    elseif strcmp(eventype,'seq1212')
        % across one, across two
        stim1 = [8,20,32,   2,14,26];
        stim2 = [12,24,36,  6,18,30];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:3)),mean(diff_stim(4:6))];
        
    elseif strcmp(eventype,'seq123456')
        % across one across two, within event across event
        stim1 = [2,20, 11,29, 3,21, 12,30];
        stim2 = [8,26, 17,35, 5,23, 14,32];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:2)),mean(diff_stim(3:4)),...
            mean(diff_stim(5:6)),mean(diff_stim(7:8))];
        
    elseif strcmp(eventype,'seq123123')
        % across one across two, within event across event
        stim1 = [11,29, 2,20, 12,30, 3,21];
        stim2 = [17,35, 8,26, 14,32, 5,23];
        
        n_stim = length(stim1);
        diff_stim = nan(n_stim,1);
        for s=1:n_stim
            sim_pre = corrs(1,stim1(s));
            sim_post = corrs(1,stim2(s));
            diff_stim(s,1) = (1-sim_post)-(1-sim_pre);
        end
        submi = [mean(diff_stim(1:2)),mean(diff_stim(3:4)),...
            mean(diff_stim(5:6)),mean(diff_stim(7:8))];
        
    end 
end