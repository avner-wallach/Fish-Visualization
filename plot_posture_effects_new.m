% plot effects of posture on tuning map
pxl2mm=ops.pxlsize;
setenv('PXLSIZE',num2str(pxl2mm));
imgpath='Z:\mormyrid_data\fish_images\';
% load('tail_network.mat');
bgcol=0;
%% images
load([imgpath,'straight_image.mat']);
str_imgs=IMGS;
load([imgpath,'tail_images.mat']);
tail_img=IMGS;
N_v=[16,14:-1:6];
tail_imgs=IMGS(N_v);
N_tail=numel(N_v);
load([imgpath,'fin_images.mat']);
N_v=[2:2:14];
fin_imgs=IMGS(N_v);
N_fin=numel(N_v);    
%%
unitind=10;
segnum=6;
zcol=21;
% clim=[0.5 1.5]; %lfp
clim=[-1 1]; %lfp
wall=0;
pcol=6;
R=5;
maxa=.3;
ind_lim=database.ind_lim(unitind,:);
F1=figure;
p=linspace(0,1,(N_tail-1)*R+3);
p=p(2:end-1);
if(pcol==1)
    q=quantile(1./eod(segnum).data(:,pcol),p);
else
    q=quantile(eod(segnum).data(:,pcol),p);
end
K=[reshape(ones(R,1)*[1:N_tail-1],1,[]) N_tail];
clear N0 Mz;
for i=1:numel(p)
    k=K(i);
%     k=floor(i/R)+1; 
    g=mod(i-1,R)/(R);
    if(pcol==6)
        if(i==numel(p))
            img=255-(tail_imgs(k).cdata);
        else
            img=255-(tail_imgs(k).cdata*(1-g)+tail_imgs(k+1).cdata*g);
        end
    else
        img=255-(str_imgs.cdata);
    end
    if(wall)
        [Mz(:,:,i),N0(:,:,i)]=plot_tuning_polar(eod(segnum),file(segnum),zcol,'t_col',pcol,...
            't_val',q(i),'clim',clim,'image',img,'maxa',maxa,'ind_lim',ind_lim,'objidx',~wall,'upsamp',.1);
    else
        Mz(:,:,i)=plot_tuning2d(eod(segnum),file(segnum),zcol,'t_col',pcol,'t_val',q(i),...
            'clim',clim,'image',img,'maxa',maxa,'upsamp',2);
    end
    colorbar('off');
    set(gca,'Color',bgcol*[1 1 1],'XColor',bgcol*[1 1 1],'YColor',bgcol*[1 1 1]);
    set(gcf,'Color',bgcol*[1 1 1]);
    set(gca,'alim',[0 12]);
    I=findobj(gcf,'Type','Image');
    Im(i).cdata=I.CData;
    S=findobj(gcf,'Type','Surface');
    set(gcf,'Position',[1 41 1920 964]);
    frame1(i)=getframe;
end

% % get correlation map
% for i=1:size(Mz,1)
%     for j=1:size(Mz,2)
%         MM=Mz(i,j,:);
%         NN=N0(i,j,:);
%         if(~isnan(sum(MM)))
% %             MM=MM(~isnan(MM));
%             f=fit(p',MM(:),'poly1','Weights',NN(:));
%             P1(i,j)=f.p1;
%             P2(i,j)=f.p2;
%             [C(i,j),pval(i,j)]=corr(p',MM(:));
%         else
%             P1(i,j)=NaN;
%             P2(i,j)=NaN;
%             C(i,j)=0;
%             pval(i,j)=1;
%         end            
%     end
% end

% Mzz=[P2 P2(:,1)];
% AL=ones(size(Mzz));
% figure;
% if(isfield(params,'image'))
%     sX=size(params.image,2);
%     sY=size(params.image,1);
%     ix=([1:sX]-sX/2)*pxl2mm;
%     iy=(-[1:sY]+sY/3)*pxl2mm;
%     imagesc(ix,iy,params.image);
%     hold on;
% end

% S=polarplot3d(Mzz,AL,'AngularRange',[-pi pi],'plottype','surfa','RadialRange',[0 150]);
% view(0,90);
% colormap(modified_jet);
% set(gca,'Clim',[-1 1])
% set(gca,'Alim',[-1 0])
