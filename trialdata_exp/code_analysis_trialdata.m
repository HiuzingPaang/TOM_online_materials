%% experiment 1
% for each column in the variable resp_test
% 1 = sequence id
% 2 = test trial id
% 3 = condition id (1=across one, 2=across two)
% 4 = accuracy (1=correct,0=incorrect)
% 5 = reaction time
% 6 = pressed button(70='f',74='j')

clear; clc; close all;
data = '.\trialdata_exp1';
subwise = nan(length(dir(data))-2,3);
for isub = 1:length(dir(data))-2
	load([data,'\sub',num2str(isub),'\sub',num2str(isub),'exp1_resp_test.mat'],'resp_test');
    subwise(isub,1) = isub;
    subwise(isub,2) = mean(resp_test(resp_test(:,3)==1,4));  % across one;
    subwise(isub,3) = mean(resp_test(resp_test(:,3)==2,4));  % across two;
end

[h,p,ci,stats] = ttest(subwise(:,2),subwise(:,3),'Tail','right');

%% experiment 2
% for each column in the variable resp_test
% 1 = sequence id
% 2 = test trial id
% 3 = condition id (1=across one, 2=across two, 3=within event, 4=across event)
% 4 = accuracy (1=correct,0=incorrect)
% 5 = reaction time
% 6 = pressed button(70='f',74='j')

clear; clc; close all;
data = '.\trialdata_exp2';
subwise = nan(length(dir(data))-2,5);
for isub = 1:length(dir(data))-2
	load([data,'\sub',num2str(isub),'\sub',num2str(isub),'exp2_resp_test.mat'],'resp_test');
    subwise(isub,1) = isub;
    subwise(isub,2) = mean(resp_test(resp_test(:,3)==1,4));  % across one
    subwise(isub,3) = mean(resp_test(resp_test(:,3)==2,4));  % across two
    subwise(isub,4) = mean(resp_test(resp_test(:,3)==3,4));  % within event
    subwise(isub,5) = mean(resp_test(resp_test(:,3)==4,4));  % across event
end

[~,p1,ci1,t1] = ttest(subwise(:,2),subwise(:,3),'Tail','right');
[~,p2,ci2,t2] = ttest(subwise(:,4),subwise(:,5),'Tail','right');
[~,p3,ci3,t3] = ttest(subwise(:,2),subwise(:,5),'Tail','right');
[~,p4,ci4,t4] = ttest(subwise(:,3),subwise(:,5),'Tail','right');