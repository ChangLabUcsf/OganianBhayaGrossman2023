function [imgall] = loadElecLoc(datapath)
    %% load data

    sSIDs = {'S5', 'S4', 'S2', 'S1', 'S7', 'S6', 'S3', 'S8', 'S9'};
    shems = {'rh',  'rh', 'rh', 'lh', 'lh', 'rh', 'rh', 'lh', 'rh'};
    sel={256, 256, 256, 256, 256, 256, 256, 128, 256};
    sz=[10, 10, 10, 75, 10, 10, 10, 10, 20];
    
    subjects=sSIDs;
    hems = shems;
    elects = sel;
    z = sz;
    
    %% load all imaging data
    ctr=1;
    for subj = subjects
        hemis = hems{ctr};
        SID=subj{1, 1};
    
        % imaging info should be the same for both VOT and Vowel
        imgPath = fullfile(datapath, SID, 'imaging');
    
        % native
        try
            imgfile = fullfile(imgPath, 'TDT_elecs_all.mat');
            if exist(imgfile, 'file'),             img = load(imgfile); end
    
        catch
            imgfile = fullfile(imgPath, 'temporal_grid.mat');
            if exist(imgfile, 'file'),             img = load(imgfile); end
            imgfile = fullfile(imgPath, 'frontal_grid.mat');
            if exist(imgfile, 'file')
                img2 = load(imgfile);
                img.elecmatrix(129:256,:) = img2.elecmatrix;
            end
        end
        imgfile = fullfile(imgPath, [SID '_'  hemis '_pial.mat']);
        if exist(imgfile, 'file'), img.cortex = getfield(load(imgfile, 'cortex'), 'cortex'); end
    
        % mni
        try
            img_mni = load(fullfile(imgPath, 'TDT_elecs_all_warped.mat'));
        catch
            warning('no warped electrode file for subject %s.' ,SID);
            img_mni = [];
        end
        img_mni.cortex = getfield(load(fullfile(datapath, 'mni_meshes', ...
            ['cvs_avg35_inMNI152_'  hemis '_pial.mat']), 'cortex'), 'cortex');
    
        if ~exist('img', 'var'), img = [];end
    
        imgall.(SID).img_native=img;
        imgall.(SID).imgNative=img;
        imgall.(SID).img_mni=img_mni;
        imgall.(SID).hemi=hemis;
        imgall.(SID).elect=elects{ctr};
        imgall.(SID).z=z(ctr);
        clear DD img img_*
        ctr=ctr+1;
    end 
end
