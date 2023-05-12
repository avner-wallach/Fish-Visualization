function [ database ] = plot_one_electrode( database,ops,segnum,fileind,lfpnum,unitnums,eod,file,IMGS,wind,expnum)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
interactive=str2num(getenv('INTERACTIVE'));


i0=numel(database.tail.ex.p1);
I=[]; k=1;
for i=1:numel(unitnums)
    uname{i}=[num2str(expnum),'_',num2str(segnum),'_',num2str(lfpnum),'_',num2str(unitnums(i))];
    iu=find(cellfun(@(x) (strcmp(x,uname{i})),database.unitname));
    if(numel(iu))
        I=[I iu];
    else
        I=[I i0+k];
        k=k+1;
    end
end

% I=i0+(1:numel(unitnums));
database.lfp.ratio(I)=plot_lfp(eod,ops,segnum,lfpnum,fileind);
[database.psth.median(I),database.psth.peak(I),database.psth.latency(I),utype]=...
    plot_trace_spikes(ops,segnum,fileind,lfpnum,unitnums,wind,expnum);
if(interactive)
    database.unit_type(I)=utype;
end

for i=1:numel(unitnums)
    database.unitname{I(i)}=uname{i};
    unitnum=unitnums(i);
    unitname=['sc',num2str(unitnum)];
    lfpname=['lfp',num2str(lfpnum)];
    unitcol=find(cellfun(@(x) (strcmp(x,unitname)),eod(segnum).fnames));
    lfpcol=find(cellfun(@(x) (strcmp(x,lfpname)),eod(segnum).fnames));
    if(interactive)
        F=figure, plot(smooth(eod(segnum).data(:,unitcol),100));
        op.WindowStyle='normal';
        prompt = {'ind start:','ind end:'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'1','inf'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans,op);    
        ind_lim=[str2num(answer{1}) str2num(answer{2})];
        database.ind_lim(I(i),:)=ind_lim;
        close(F);
    else
        ind_lim=database.ind_lim(I(i),:);
%         ind_lim=[1 inf];
    end

    t_col=6;
    [database.tail.ex.p1(I(i)),database.tail.ex.p2(I(i)),...
        database.tail.re.p1(I(i)),database.tail.re.p2(I(i)),...
        database.tail.ex.rsq(I(i)),database.tail.re.rsq(I(i)),...
        database.tail.ex.vr(I(i)),database.tail.re.vr(I(i))]=...
        plot_ex_re(eod,file,ops,segnum,unitname...
        ,lfpname,t_col,IMGS,ind_lim);

    t_col=1;
    [database.freq.ex.p1(I(i)),database.freq.ex.p2(I(i)),...
        database.freq.re.p1(I(i)),database.freq.re.p2(I(i)),...
        database.freq.ex.rsq(I(i)),database.freq.re.rsq(I(i)),...
        database.freq.ex.vr(I(i)),database.freq.re.vr(I(i))]=...
        plot_ex_re(eod,file,ops,segnum,unitname...
        ,lfpname,t_col,IMGS,ind_lim);
    
    t_col=find(cellfun(@(x) (strcmp(x,lfpname)),eod(segnum).fnames));
    [database.io.ex.p1(I(i)),database.io.ex.p2(I(i)),...
        database.io.re.s(I(i)),database.io.re.m(I(i)),...
        database.io.ex.rsq(I(i)),database.io.re.rsq(I(i))]=...
        plot_ex_re(eod,file,ops,segnum,unitname...
        ,lfpname,t_col,IMGS,ind_lim);

    %input-output rasters/psths  
    alllfp=find(cellfun(@(x) numel(strfind(x,'lfp')),eod(segnum).fnames));    
    rastnum=unitcol-alllfp(end);
    plot_sorted_raster(eod(segnum).raster{rastnum},eod(segnum).data(:,lfpcol),ops,database.ind_lim(I(i),:),1e3,7,.5);
    
    %info
    [database.info.Izt(I(i)),database.info.Izt_xy(I(i)),...
        database.info.Izxy(I(i)),database.info.Izxy_t(I(i))]=...
        compute_info(eod(segnum),file(segnum),unitcol,'t_col',lfpcol);
end

end

