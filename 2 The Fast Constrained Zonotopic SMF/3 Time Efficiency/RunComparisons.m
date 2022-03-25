%%% Run Comparisons

close all;
clear;
clc;


%%  Save
ExampleB_1_Comparison_DT_n10
save('experiment_comparison_n_10_20220322')

ExampleB_1_Comparison_DT_n15
save('experiment_comparison_n_15_20220322')

ExampleB_1_Comparison_DT_n20
save('experiment_comparison_n_20_20220322')


%%  Load
load('experiment_comparison_n_10_20220322')
load('experiment_comparison_n_15_20220322')
load('experiment_comparison_n_20_20220322')

averaged_running_time_description_OITCZSMF_n10
averaged_running_time_description_CZSMF_n10
ratio_n10

averaged_running_time_description_OITCZSMF_n15
averaged_running_time_description_CZSMF_n15
ratio_n15

averaged_running_time_description_OITCZSMF_n20
averaged_running_time_description_CZSMF_n20
ratio_n20