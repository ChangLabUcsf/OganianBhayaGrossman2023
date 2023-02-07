function [out] = addNeuralResp(SID, corpus, datapath)
    % pred - whether you want strf predictions saved in the out structure
    % cstrf - a strf model structure, output of strf code
    % cout - out structure with data; predictions are calculated for every single sentence in this our structure
    % sentdet - outstructure with sentence features, must contain the feature fields that are used in the cstrf model. used to create feature matrix for prediction.

    corpus_path = fullfile(datapath, 'pt_data', SID, corpus);

    outfile = dir(fullfile(corpus_path,'*HilbAA_70to150_8band*out_resp_log.mat'));
    if isempty(outfile)
        warning('no out file for %s', fullfile(corpus_path, '*out_resp_log.mat'));
        out = [];
        return
    end
    out = getfield(load(fullfile(outfile(1).folder, outfile(1).name)), 'out');
end