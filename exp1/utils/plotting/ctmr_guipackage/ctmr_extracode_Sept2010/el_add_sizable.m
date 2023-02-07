function el_add_sizable(els,r2,varargin)
% function el_add_withr2size(els,r2)
% 
% els: rows = electrodes, columns = xyz
% r2: for size and color of electrodes
% varargin:
% {1} : maximum for scale (if {''} absmax)
% example: el_add_withr2size(els,r2,1)

hold on
% % black circle around electrode:
% plot3(els(:,1),els(:,2),els(:,3),'.','Color', 'k','MarkerSize',msize+5)
% plot3(els(:,1),els(:,2),els(:,3),'.','Color', [0.7 0.7 0.7],'MarkerSize',msize)

%gray2green
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,2)=[0:0.7/9:0.7]';
% cm1(10:100,2)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,1)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,1)=[0:1/80:1]';
% cm1=[repmat([0 1 0],100,1)];

cm1(1:100,2)=[0:1/99:1]';
cm1=[repmat([0 0 0],100,1)];
cm1(1:10,1)=[0:0.7/9:0.7]';
cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
cm1(1:10,2)=[0]';
cm1(1:10,3)=[0]';
cm1(20:100,2)=[0:1/80:1]';


cm2=[repmat([0 0 0],100,1)];
cm2(1:10,3)=[0:0.7/9:0.7]';
cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
cm2(1:10,2)=[0]';
cm2(1:10,1)=[0]';
cm2(20:100,2)=[0:1/80:1]';

% cm2=cm2(end:-1:1,:);
% cm=[cm2; cm1];
elsize=[15:(45-15)/(100-1):45];

maxr2=round(max(r2));
if abs(round(min(r2)))>maxr2
    maxr2=abs(round(min(r2)));
end

if isempty(varargin)
    r2=r2/maxr2;% scaled to absmax
else
    r2=r2/varargin{1};% scaled to varargin{1}
end

% electrode with r2:
for k=1:length(els)
    if ~isnan(r2(k))
        if abs(r2(k))>0.01
            ind_color=abs(round(100*r2(k)));
            if ind_color>100
                ind_color=100;
            end
            elsize_r2=elsize(ind_color);
            if r2(k)>0.01
                elcol_r2=cm1(ind_color,:); 
                plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
                plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
            elseif r2(k)<0.01
                elcol_r2=cm2(ind_color,:);
                plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize_r2)
                plot3(els(k,1),els(k,2),els(k,3),'.','Color',elcol_r2,'MarkerSize',elsize_r2-5)
            end
        else
            plot3(els(k,1),els(k,2),els(k,3),'.','Color','k','MarkerSize',elsize(1))
        end
    end
end


%% plot color scale

 
% cm1(1:100,2)=[0:1/99:1]';
% cm1=[repmat([0 0 0],100,1)];
% cm1(1:10,1)=[0:0.7/9:0.7]';
% cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
% cm1(1:10,2)=[0]';
% cm1(1:10,3)=[0]';
% cm1(20:100,2)=[0:1/80:1]';
% 
% cm2=[repmat([0 0 0],100,1)];
% cm2(1:10,3)=[0:0.7/9:0.7]';
% cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
% cm2(1:10,2)=[0]';
% cm2(1:10,1)=[0]';
% cm2(20:100,2)=[0:1/80:1]';
% 
% % cm2=cm2(end:-1:1,:);
% % cm=[cm2; cm1];
% elsize=[15:(45-15)/(100-1):45];
% 
% figure('Color',[1 1 1],'Position',[30 50 50 300]),hold on
% for k=1:100
%     plot(1,k,'.','MarkerSize',elsize(k),'Color',cm1(k,:))
% end
%     
% for k=1:100
%     plot(1,-k,'.','MarkerSize',elsize(k),'Color',cm2(k,:))
% end
% 
% ylim([-120 120])
% set(gcf, 'PaperPositionMode', 'auto');
% print('-painters','-r300','-dpng',strcat(['./figures/ecog/colorscale_el_add_sizable']));
% print('-painters','-r300','-depsc',strcat(['./figures/ecog/colorscale_el_add_sizable']));

