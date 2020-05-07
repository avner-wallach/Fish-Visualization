% plot raw data example
d='20190624';
fnum='35';
samplerate=3e4;
chan=[12];
channum=16;
LFPc=[236,157 118]/255;
Uc=[210 149 0;74 133 34]/255;
% [t,amp]=read_bonsai_binary(['Z:\mormyrid_data\',d,'\amp_',fnum],30e3,channum,256,chans,'adc');
% load(['Z:\mormyrid_data\',d,'\data_',fnum]);
ind=[1:5e4]+1e4;
t0=t(ind(1));
figure;
H=plot(t(ind)-t0,amp(ind));
set(H,'Color',0.5*[1 1 1]);
hold on;

%% LFP
idx=find(inrange(data.EOD.t+data.FILE.offset,[t(ind(1)),t(ind(end))]));
for i=1:numel(idx)
    T=data.EOD.t(idx(i))+([1:size(data.EOD.traces,2)]-2)/samplerate+data.FILE.offset+1e-3;
    H=plot(T-t0,data.EOD.traces(idx(i),:,chan));
    set(H,'Color',LFPc);
end

%% units
rast=data.SPIKES.raster;
idx=find(inrange(rast(:,1)+data.FILE.offset,[t(ind(1)),t(ind(end))]));
V=150;
for i=1:numel(idx)
    T=rast(idx(i),1)+data.FILE.offset;
    H=plot(T-t0,V,'.');
    set(H,'Color',Uc(rast(idx(i),2),:),'MarkerSize',12);
end



