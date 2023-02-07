function ctmr_vox_plot(cortex,xyz,weights,ssize,v)

% ctmr_vox_plot(cortex,xyz,weights,ssize)
% xyz and weights output from voxplot_func_gm
% ssize=2; % size of voxel squares on surface
% v='l'; % 'l' or 'r' left or right
% ask for load cortex

brain=cortex.vert;

if length(weights)~=length(xyz(:,1))
    error('you sent a different number of weights than xyz (perhaps a whole matrix instead of vector)')
end

c=zeros(length(cortex.vert),1);
for k=1:length(xyz(:,1))
    b_z=abs(brain(:,3)-xyz(k,3));
    b_y=abs(brain(:,2)-xyz(k,2));
    b_x=abs(brain(:,1)-xyz(k,1));
    
    d=b_z<ssize & b_y<ssize & b_x<ssize;
    d=d*weights(k); % no smoothing

    c=max(c,d); %overlap is going to maximum
end

a=tripatch(cortex, '', c);
shading interp;
a=get(gca);
%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
d=a.CLim;
set(gca,'CLim',[0 max(weights)]);

l=light;
%load in colormap
% cm=colormap('hot');
% cm=[.7 .7 .7;cm];
cm1=[repmat([1 0 0],100,1)];
cm1(1:40,1)=[0.7:(1-0.7)/39:1];
cm1(1:40,2)=[0.7:-0.7/39:0];
cm1(1:40,3)=[0.7:-0.7/39:0];
cm1(41:100,2)=[0:1/59:1];
cm=cm1;
colormap(cm);

lighting gouraud; %play with lighting...
% material dull;
material([.3 .8 .1 10 1]);
axis off
set(gcf,'Renderer', 'zbuffer')

if v=='l'
view(270, 0);
%view(-92,32);
set(l,'Position',[-1 0 1])        
elseif v=='r'
view(90, 0);
set(l,'Position',[1 0 1])        
end
