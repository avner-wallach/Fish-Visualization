function [FF,ax_grid]=visualize_fish_modeling(vid,txt,num,ss,framenum,circle,bgimage)
    th=0.8;
    sfactor=3;
    k=10;
    ssize=25;
    rpoles=10;
    bgcol=[1 1 1];
    
    %model params
    N=200; %number of positive charges
    params.grid_M=1e2;        
    params.object_x=[0 0];
    params.object_R=0;
    params.object_c=0;
    p0=[circle(1)+circle(3)/2 circle(2)+circle(4)/2]; %tank center
    r=mean(circle(3:4))/2; %tank radius    
    
    SE1 = strel('disk',7,4);   
    SE2 = strel('disk',20,4);    
    imth=80;
%     [num,txt,raw] = xlsread([trackfile,'.csv']);
%     txt=txt(2,2:end);
    [fish,seg,coornames,data]=get_posture(txt,num);
    fnames1={'mouth','Trunk1','Trunk2','Tail1','Tail2','CaudalFork'};      
%     vid = VideoReader([vidname,'.avi']);    
    params.grid_center=p0;
    params.r_max=r;
    
    Rcirc=1.1;
    
    currentFrame=ss;
    F=figure;
    set(F,'Color',bgcol,'Position',[1 41 1920 964]);
%     COL=colormap('lines');
    COL=colormap(invert_map(brewermap(64,'BrBG')));
    A1=axes;
    set(A1,'Color',bgcol);

    for i=1:framenum
        axes(A1);
        frame = rgb2gray(read(vid, currentFrame));        
        if(i==1)
            ix=[1:size(frame,2)];
            iy=[1:size(frame,1)];
            [IX,IY]=meshgrid(ix,iy);
            imask=repmat(uint8(((IX-p0(1)).^2/r^2+(IY-p0(2)).^2/r^2)<Rcirc)*255,1,1,3);
        end
        
        ax_grid={ix,iy};   
        if(nargin==7)
            FR=(abs(double(frame)-double(bgimage(:,:,1)))>45);
        else
            FR=ones(size(frame));
        end
        I=imagesc(ix,iy,repmat((frame),[1 1 3]));
        set(I,'AlphaData',FR);
        hold on;
        imagesc(ix,iy,255*bgcol(1)-imask,'AlphaData',~imask(:,:,1))
        
        %get main axis points
        X=[]; Y=[];        
        for j=1:numel(fnames1)
            xy= data.(fnames1{j}).xy(currentFrame,:);
            c= data.(fnames1{j}).c(currentFrame);
            if(c>th)
              X=[X;xy(1)];
              Y=[Y;xy(2)];
            else
              X=[X;NaN];
              Y=[Y;NaN];
            end                
        end
        X=X*sfactor;
        Y=Y*sfactor;
        
        %find edges
        fr1=imdilate(imclose(imopen(frame<imth | frame>200,SE1),SE2),SE1);
        fr1(((IX-p0(1)).^2+(IY-p0(2)).^2)>r^2)=0;
        CC=bwconncomp(fr1);
        A=regionprops('table',CC,'Area')
        [marea,iarea]=max(A.Area);
        fr2=zeros(size(fr1));        
        fr2(CC.PixelIdxList{iarea})=1;
        [Ce,he]=contour(ix,iy,fr2,1);        
        delete(he);
        
        %circle edges
%         az=linspace(0,2*pi,40)';
%         Ce(1,:)=p0(1)+0.95*r*cos(az);
%         Ce(2,:)=p0(2)+0.95*r*sin(az);
                
        %resample points
        [X,Y]=resample_curve(X,Y,N);
        X=X(rpoles:end-rpoles);
        Y=Y(rpoles:end-rpoles);
        Q=[ones(numel(X)-1,1)/(numel(X)-1);-1];
        
        % get wall reflection
        d=normrows([X-p0(1) Y-p0(2)]); %distance of pole from tank center
        alpha=atan2(Y-p0(2),X-p0(1)); %azimouth of pole relative to center
        L=r^2./d;
        Xw=p0(1)+L.*cos(alpha);
        Yw=p0(2)+L.*sin(alpha);
        Qw=Q;%*r./d;

        %get potential map
        [ Z_V,Z_E,Xg,Yg ] = get_potential_field_map([X Y;Xw Yw],[Q;Qw],params)
        Z_V(((Xg-p0(1)).^2+(Yg-p0(2)).^2)>r^2)=nan;
        Z_E(((Xg-p0(1)).^2+(Yg-p0(2)).^2)>r^2)=nan;
%         [IN,ON]=inpolygon(Xg,Yg,Ce(1,:),Ce(2,:));
%         Z_E(IN & ~ON)=nan;

%         H=scatter(Xw(1:k:end),Yw(1:k:end),ssize,COL(1,:));        
%         set(H,'MarkerFaceAlpha',.5,'MarkerFaceColor',COL(1,:),'MarkerEdgeAlpha',0);
        c=quantile(Z_V(:),linspace(0.05,0.95,18));
        [C]=contour(Xg,Yg,Z_V,c,'LineWidth',3);
        S=streamline(Xg',Yg',Z_E(:,:,1)',Z_E(:,:,2)',Ce(1,1:2*k:end),Ce(2,1:2*k:end));
        S=[S;streamline(Xg',Yg',-Z_E(:,:,1)',-Z_E(:,:,2)',Ce(1,1:2*k:end),Ce(2,1:2*k:end))];
        set(S,'Color',.25*[1 1 1],'LineWidth',1);
        for s=1:numel(S)  
            ind=find(inpolygon(S(s).XData,S(s).YData,Ce(1,:)',Ce(2,:)'));
%             ind=[ind1(:);ind2(:)];
            S(s).XData(ind)=nan;
            S(s).YData(ind)=nan;
        end
        H=scatter(X(1:k:end-1),Y(1:k:end-1),ssize,COL(end,:));        
        set(H,'MarkerFaceAlpha',1,'MarkerFaceColor',COL(end,:),'MarkerEdgeAlpha',0);
        H=scatter(X(end),Y(end),ssize,COL(1,:));        
        set(H,'MarkerFaceAlpha',1,'MarkerFaceColor',COL(1,:),'MarkerEdgeAlpha',0);        
        
        

        hold off;
        axis('image');
        xs=size(frame,2); ys=size(frame,1);
%         set(gca,'Xlim',[-xs/4 xs*1.24],'Ylim',[-ys/4 ys*1.25]);
        set(gca,'XColor',bgcol,'YColor',bgcol);
%         set(gca,'Xlim',[-size(frame,2) 2*size(frame,2)],'Ylim',[-size(frame,1) 2*size(frame,1)],'XColor',[0 0 0],'YColor',[0 0 0]);
        set(gca,'Color',bgcol);
        R=rectangle('Position',circle,'Curvature',[1 1]);
        
        drawnow;
        pause(0.01);
        FF(i)=getframe(gca);
        
%         writeVideo(vidout,FF(i));
        currentFrame=currentFrame+1;        
    end
    
end

                

