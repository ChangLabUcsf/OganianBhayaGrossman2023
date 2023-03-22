%% startup

addpath(genpath('utils'));
datapath = '../../data_anon/exp1/';

bef=20;
aft=50;

SIDs = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S9'}; %'S8'

Dvow = loadDD('dimex', bef, aft, SIDs, datapath);
Tvow = loadDD('timit', bef, aft, {}, datapath);

model = 'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2';

% load in all beta model versions
beta_model = 'curveFits_vowelftest_audF1F2binwtF0_spanish_anon.mat'; 
load([datapath '/beta_models/' beta_model])
