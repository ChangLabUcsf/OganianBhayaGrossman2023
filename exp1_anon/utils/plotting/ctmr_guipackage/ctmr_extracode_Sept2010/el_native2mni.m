function el_mni = el_native2mni(el_native,norm_inv_sn,norm_anat)
% function converts electrode positions in native (subject)
% space to MNI space. 
% usage: el_mni = el_native2mni(el_native)
%
% input: el_native = electrode coordinates in native space 
%           (not voxel space, rows: electrodes, columns: xyz)
%        norm_inv_sn = inverse normalization parameters from SPM
%           unified segmentation
%        norm_anat = normalized anatomical for the same subject 
% 
% output: el_mni electrodes in MNI space  
%           (rows: electrodes, columns: xyz)
%
% Uses the function get_orig_coord5 from SPM5


% electrodes from native space to MNI voxel space coordinates
el_mni = get_orig_coord5(el_native,norm_inv_sn,norm_anat);

% convert MNI voxels space to MNI coordinates
brain_info=spm_vol(norm_anat); 
el_mni=(el_mni*brain_info.mat(1:3,1:3)')+...
    repmat(brain_info.mat(1:3,4),1,length(el_mni(:,1)))';
