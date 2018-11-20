function SaveDownstSpectra(handles)

n=5; %number of windows
windT=0.8; %duration of time windows, in seconds.
sp=handles.parentHandles.sp;
if length(handles.downstates.fromI)>=2
    from1=handles.downstates.fromI(2);
else
    warning('only 1 upstate')
end
from=from1+(0:(n-1))*floor((length(handles.DS)-from1)/n);
to=from+round(windT/sp);

params.tapers=[5 9];
params.Fs=1/sp;
params.fpass=[0 1000];

for i=1:n
    [pw(i,:),fr(i,:)]=mtspectrumc(handles.DS(from(i):to(i)),params);
end
tmpOut=zeros(n+1,length(fr(1,:)));
tmpOut(1,:)=fr(1,:);
tmpOut(2:n+1,:)=pw;

pw=pw';
fmt='';
for i=1:n+1
    fmt=[fmt '%12.10f ']; %one column for each spectrum (window)
end
fmt=[fmt '\n'];
currCh=num2str(handles.parentHandles.currentCh);
currTr=num2str(handles.parentHandles.currentTrial);
fid=fopen([handles.parentHandles.dir_in 'DSspectra_' num2str(n) '_windows_ch' currCh '_tr' currTr '.dat'],'w');
fprintf(fid,fmt,tmpOut);
fclose(fid);
