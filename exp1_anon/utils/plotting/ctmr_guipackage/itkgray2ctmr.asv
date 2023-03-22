function result = itkgray2ctmr(seg_nii, hemisphere, subject)
%  Convert ITKGray segmentations for use with SPM's ctmr package.
%
% Usage:
%   result = itkgray2ctmr([seg_nii=prompt], [hemisphere='both'], [subject=prompt])
%
% Inputs:
%   
% seg_nii: either a loaded NII struct (using load_nii) or a path to a NIFTI
% file containing the ITKGray segmentation.
%
% hemisphere: either 'left', 'right', or 'both' (or 'l', 'r', or 'b').
% Specifies which hemisphere to convert. [Default: 'both']
%
% subject: subject name. Will save the data in a directory called
% ['data_' subject] relative to the current dir.
%
% Outputs:
%
%  Produces 2 files -- [subject]_white.nii and [subject]_grey.nii, which
%  are suitable for input into the ctmr package (starting w/ get_mask) 
%  in place of the spm5 segment files.
%
% result: struct with paths to the saved files, and a status flag (true if
% the code was successful, false otherwise).
%
% POSSIBLE TODO: if results of this aren't good, one parameter to change
% is the value of the WM and GM masks -- they may interact with the spm
% smoothing / thresholding code. 
%
% ras, 09-02-2010

if notDefined('seg_nii')
    seg_nii = mrvSelectFile('r', 'nii', 'Select a T1 segmentation file');
end

if notDefined('hemisphere')
    hemisphere = 'both';
end

if notDefined('subject')
    subject = input('Subject name: ', 's');
end

% allow hemisphere to be specified w/ just a letter
if strncmp(lower(hemisphere), 'l', 1)
    hemisphere = 'left';
elseif strncmp(lower(hemisphere), 'r', 1)
    hemisphere = 'right';
else
    hemisphere = 'both';
end

data_dir = fullfile(pwd, ['data_' subject]);
result.status = false;
result.grey_file = [subject '_' hemisphere '_grey.nii'];
result.white_file = [subject '_' hemisphere '_white.nii'];
result.grey_path = fullfile(data_dir, result.grey_file);
result.white_path = fullfile(data_dir, result.white_file);
result.mask_val = uint8(127);

try
    % make sure we have a data directory
    ensureDirExists(data_dir);
    cd(data_dir);
    
    % load the segmentation
    if ~isstruct(seg_nii)
        seg = load_nii(seg_nii);
    end
    
    % the WM and GM values we keep depend on the hemispheres we're
    % interested in:
    switch hemisphere
        case 'left',
            wm_vals = 3;
            gm_vals = 5;
        case 'right',
            wm_vals = 4;
            gm_vals = 6;
        case 'both',
            wm_vals = [3 4];
            gm_vals = [5 6];
    end
    
    % convert the white matter
    % (We use a hack here, for compatibility with MRIcron and probably
    % other packages: we use the SPM header format, hacked off the 
    % segmentation file, and spm's save code.)    
    wm = result.mask_val .* uint8(ismember(seg.img, wm_vals));
    hdr = spm_vol(seg_nii);
    hdr.fname = result.white_file;
    hdr.dim = size(wm);
    spm_write_vol(hdr, wm);
    fprintf('Saved white matter as %s.\n', result.white_path);
    
    % convert the grey matter
    gm = result.mask_val .* uint8(ismember(seg.img, gm_vals));
    hdr = spm_vol(seg_nii);
    hdr.fname = result.grey_file;
    hdr.dim = size(gm);
    spm_write_vol(hdr, gm);
    fprintf('Saved grey matter as %s.\n', result.grey_path);
    
    fprintf('[%s]: done! \t %s', mfilename, datestr(now));
    
    cd ..
    
catch
    disp('itkgray2ctmr failed!')
    rethrow(lasterror)
end

return


%% debug / fix if this doesn't work -- 
% resave the output of get_mask using spm_write_vol
genFile = mrLoad;  % select the output .img from get_mask
newhdr = genFile.hdr; %grab header info and place it in a new variable
newhdr.fname = 'new_mask20_35p5.img'; %Replace with new file name.img
tempIm = spm_vol('tempg.img'); % to get the header from an intermediate file
newhdr.descrip = 'Tissue class 2 - conv(10,10,10)'; %might be unnecessary
newhdr.mat = tempIm.mat; %place the intermediate file's .mat into the new header
newhdr.dim = tempIm.dim; %place teh intermediate file's .dim into the new header
spm_write_vol(newhdr, genFile.data); %write the new volume with the new header and original data

% hdr = spm_vol('tempw.img');  % or any of the intermediate files
% hdr.fname = 'resave_test.img'; % change to what you want
% spm_write_vol(hdr, sfc.data);  % somehow this works ... :P