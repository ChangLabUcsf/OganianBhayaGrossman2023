function ctmr_med_plot(cortex,coords_dark)
% function [electrodes]=ctmr_gauss_plot(cortex,electrodes,weights)
% projects electrode locations onto their cortical spots in the 
% left hemisphere and plots about them using a gaussian kernel
% for only cortex use: 
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)

%     Copyright (C) 2009  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%   Version 1.1.0, released 26-11-2009



%load in colormap
load('loc_colormap')

brain=cortex.vert;
v='l';
% %view from which side?
% temp=1;
% while temp==1
%     disp('---------------------------------------')
%     disp('to view from right press ''r''')
%     disp('to view from left press ''l''');
%     v=input('','s');
%     if v=='l'      
%         temp=0;
%     elseif v=='r'      
%         temp=0;
%     else
%         disp('you didn''t press r, or l try again (is caps on?)')
%     end
% end

c=zeros(length(cortex(:,1)),1);
c(coords_dark)=1;

% c=(c/max(c));
a=tripatch(cortex, '', c');
shading interp;
a=get(gca);
%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
d=a.CLim;
set(gca,'CLim',[-max(abs(d)) max(abs(d))])
l=light;
colormap(cm)
lighting gouraud; %play with lighting...
% material dull;
material([.3 .8 .1 10 1]);
axis off
set(gcf,'Renderer', 'zbuffer')

if v=='l'
view(270, 0);
set(l,'Position',[-1 0 1])        
elseif v=='r'
view(90, 0);
set(l,'Position',[1 0 1])        
end
% %exportfig
% exportfig(gcf, strcat(cd,'\figout.png'), 'format', 'png', 'Renderer', 'painters', 'Color', 'cmyk', 'Resolution', 600, 'Width', 4, 'Height', 3);
% disp('figure saved as "figout"');