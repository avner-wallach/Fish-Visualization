function make_video_from_imgs(imgs,vidname,mode,N)
    presentpath='Z:\presentations\';    
%     presentpath='C:\Users\Avner\iCloudDrive\Documents\Mormyrid_Data\presentations\group meeting 111518\';
    rep=1;
%     load([presentpath,filename]);
    
    for i=1:size(imgs,1)
        i
        vidout=VideoWriter([presentpath,vidname,num2str(N),'.avi']);
        vidout.FrameRate=10;
        open(vidout);
        for k=1:rep
            for j=1:size(imgs,2)            
                writeVideo(vidout,imgs(i,j));
            end
            if(strcmp(mode,'freq'))
                for k=1:5
                    writeVideo(vidout,imgs(i,j));
                end
            elseif(strcmp(mode,'posture'))
                for j=(size(imgs,2)-1):-1:2
                    writeVideo(vidout,imgs(i,j));
                end
            end

        end
        close(vidout);
%         F=figure;
%         imshowpair(imgs(i,1).cdata,imgs(i,end).cdata);
        
    end

end
   
