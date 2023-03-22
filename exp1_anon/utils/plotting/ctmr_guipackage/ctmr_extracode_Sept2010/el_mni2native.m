function el_native = el_mni2native(el_mni,norm_sn,native_anat)
% function converts electrode positions in MNI space 2 native 
% (subject) space. 
% usage: el_native = el_mni2native(el_mni)
%
% input: el_mni = coordinates in MNI space 
%        norm_sn = normalization parameters from SPM
%           unified segmentation
%        native_anat = anatomical for the same subject 
% 
% output: el_native coordinates in native space  
%           (rows: electrodes, columns: xyz)
%
% Uses the function get_orig_coord5 from SPM5


% electrodes from MNI space to native voxel space coordinates
el_native = get_orig_coord5(el_mni,norm_sn,native_anat);

% convert MNI voxels space to MNI coordinates
brain_info=spm_vol(native_anat); 
el_native=(el_native*brain_info.mat(1:3,1:3)')+...
    repmat(brain_info.mat(1:3,4),1,length(el_native(:,1)))';
