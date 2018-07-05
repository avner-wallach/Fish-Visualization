function plot_pcaspace(frame,eod,file)
%
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
sesspath=[datapath,'\',sdate,'\'];
%%
    selfile=1;
    filenum=num2str(selfile);
    vid = VideoReader([sesspath,'video_',filenum,'.avi']);
    % get data file
    datafile=[sesspath,'data_',filenum];
    s=load(datafile);        
    data=s.data;        
    %get tracking file
    trackfile=[sesspath,'trackedFeaturesRaw_',filenum];
    [num,txt,raw] = xlsread([trackfile,'.csv']);    

% get neutral fish image
    presentpath='Z:\AW\presentations\';    
    selected_imgs=[];   
%     set(F,'Position',[0 0 1000 1000]);
    int=0.025;
    M=  [0 int;
        0.25-int/2 0.25+int/2;
        0.5-int/2 0.5+int/2;
        0.75-int/2 0.75+int/2;
        1-int 1]
    pcvec=[14:20];
    for p=1:numel(pcvec)
        zcol=pcvec(p);
        ctlcol=setdiff(pcvec,zcol);  
        F=figure;   
        imgs=get_posture_images(frame,file,zcol,ctlcol,M,'imgnum',10,'track','off',...
            'selfile',selfile,...
            'vid',vid,'txt',txt,'data',data,'num',num);
        close(F);
        plot_imgs_dlg;
    end
        
    save([presentpath,'presentation_imgs'],'selected_imgs');
    
    function plot_imgs_dlg
        F=figure;
        k=1;
        for i=1:size(imgs,1)
            for j=1:size(imgs,2)
                subplot(size(imgs,1),size(imgs,2),k);
                imshow(imgs(i,j).cdata);
                k=k+1;
            end
            prompt{i} = ['Line ',num2str(i),':'];
            defans{i}='1';
        end
        
        dlg_title = 'Choose images';
        num_lines = size(imgs,1);        
        opt.WindowStyle='normal';
        answer = inputdlg(prompt,dlg_title,num_lines,defans,opt);
        for i=1:size(imgs,1)
            selected_imgs=[selected_imgs imgs(i,str2num(answer{i}))];
        end
        close(F);
    end
end