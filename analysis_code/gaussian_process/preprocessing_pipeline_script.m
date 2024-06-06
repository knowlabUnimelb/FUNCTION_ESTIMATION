%% Preprocessing pipeline for FUNCTION ESTIMATION Experiment
%  Subject numbers are listed in socfun_subject_numbers.txt


clear all
clc

subs = dlmread('socfun_subject_numbers.txt');

%% Step 1
for i = 1:numel(subs) 
%     Examine all fits, rotate and organize first & second attempts
        preprocess_1(subs(i), true, true)
%         preprocess_1(subs(i), true, false)
    
    pause
%     return
end

%% Step 2
for i = 1:numel(subs)
    disp(['Starting Subject #: ' num2str(subs(i))])
    preprocess_2_correct_multiple_responses(subs(i))
    disp(['Ending Subject #: ' num2str(subs(i))])
end

%% Step 3
for i = 1:numel(subs)
    preprocess_3(subs(i), true, true)
    pause(1);
end