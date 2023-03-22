
% el_native = electrodes localized on the individual subject CT
% norm inv = the inversed normalization parameters from the SPM unified
% segmentation procedure
% norm_anat = a normalized anatomical scan
% 
% el_mni = electrode coordinates in MNI space

el_mni = el_native2mni(el_native,norm_inv_sn,norm_anat);

% and reverse
el_native = el_mni2native(el_mni,norm_sn,native_anat);


