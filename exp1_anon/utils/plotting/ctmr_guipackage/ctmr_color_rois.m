subj = 'EC41';
hemi = 'rh';
roi = 'superiortemporal';
clr = [0 1 0.9];
elecsize = 80;
rootdir = '/Users/mattleonard/Documents/Research/data/MRI';

load([rootdir '/' subj '/Meshes/' subj '_' hemi '_pial.mat']);
load([rootdir '/' subj '/elecs/hd_grid.mat']);

dat = h5read([rootdir '/DM_data/elecs.hdf5'],['/' subj]);
dat(isnan(dat)) = 0;

ctab_fid = fopen([rootdir '/aparc.annot.ctab']);
ctab = textscan(ctab_fid,'%d%s%d%d%d%d');
roi_names = ctab{2};
roi_idx = find(strcmpi(roi_names,roi));

fid = fopen([rootdir '/' subj '/' hemi '.aparc.annot.dpv']);
vert = textscan(fid,'%d%d%d%d%d');
vert = vert{5};

roi_verts = find(vert == roi_idx);

figure;
ctmr_gauss_plot(cortex,[0 0 0],0,hemi);

kids = get(gca,'Children');
kids(2).FaceVertexCData(roi_verts,:) = 'r';
% kids(2).FaceAlpha = 'flat';
% kids(2).FaceVertexAlphaData(roi_verts) = 0.5;
colormap([[0.9412    0.9412    0.9412] ; clr]);
caxis([0 1]);

dat1 = dat / (max(dat) - min(dat));

if strcmpi(hemi,'lh')
    offset = -5;
else
    offset = 5;
end
cmap = cbrewer('seq','Reds',101);
for i = 1:size(elecmatrix,1)
    if dat1(i) == 0
        scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
            elecsize,'k');
    else
        scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
            elecsize,cmap(round(dat1(i)*100)+1,:),'filled','MarkerEdgeColor','k');
    end
    hold on;
end