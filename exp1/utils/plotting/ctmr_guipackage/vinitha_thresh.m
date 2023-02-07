function vinitha_thresh




% load data in 
w1=spm_vol(spm_select(1,'image','select white matter image',''));
w=spm_read_vols(w1);

% rename to get new file
ff=w1.fname; ff((end-4):end)=[]; ff=[ff '_thresh.nii']; w1.fname=ff;

% thresholds data
w=w>.25;

% writes data
spm_write_vol(w1,w)





















