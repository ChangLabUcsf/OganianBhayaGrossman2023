%% startup

addpath(genpath('utils'));
datapath = '../../data for sharing/exp1/';

bef=20;
aft=50;

SIDs = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S9'}; %'S8'

Dvow = loadDD('dimex', bef, aft, SIDs, datapath);
Tvow = loadDD('timit', bef, aft, {}, datapath);

model = 'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2';

% load in all beta model versions
beta_model = 'curveFits_vowelftest_audF1F2binwtF0_spanish_anon.mat'; 
load([datapath '/beta_models/' beta_model])

%% ensure only repeat sentences exist in HGA out structures
% for s = SIDs
%     corpus = 'dimex';
%     corpus_path = fullfile(datapath, 'pt_data', s{1}, corpus);
% 
%     outfile = dir(fullfile(corpus_path,'*HilbAA_70to150_8band*out_resp_log.mat'));
%     if isempty(outfile)
%         warning('no out file for %s', fullfile(corpus_path, '*out_resp_log.mat'));
%         out = [];
%         return
%     end
%     out = getfield(load(fullfile(outfile(1).folder, outfile(1).name)), 'out');
%     
%     repeats = cellfun(@(x) length(x), {out.Trials});
%     out = out(repeats>2);
% 
%     save(fullfile(outfile(1).folder, outfile(1).name), 'out');
% end