%% startup

addpath(genpath('utils'));
datapath = '../data/';

bef=20;
aft=50;

SIDs={'EC172', 'EC163', 'EC105', 'EC214', 'EC100', 'EC203', 'EC152', 'EC252'};

Dvow = loadDD('dimex', bef, aft, SIDs, datapath);
Tvow = loadDD('timit', bef, aft, {}, datapath);

model = 'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2';

% load in all beta model versions
beta_model = 'curveFits_vowelftest_audF1F2binwtF0_spanish.mat'; 
load([datapath '/beta_models/' beta_model])
