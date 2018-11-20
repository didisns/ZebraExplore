function [ output_args ] = plotPWhisto(hObject, handles)

% plotPWhisto
%   Detailed explanation goes here

handles = zebraguidata(hObject);
% what is the visualised BP?
if get(handles.analyzePWdata.BP1display,'Value'), BP = 1;
end
if get(handles.analyzePWdata.BP2display,'Value'), BP = 2;
end
if get(handles.analyzePWdata.BP3display,'Value'), BP = 3;
end
if get(handles.analyzePWdata.HPdisplay,'Value'), BP = 4;
end

% what is the visualized channel?
nChDisplay = get(handles.analyzePWdata.selectCh, 'Value');
if (nChDisplay == 0), nChDisplay = 1;
end
axes (handles.analyzePWdata.PWdistribution);
cla
tmpx = squeeze(handles.centers{BP}(nChDisplay,1,:));
tmpy = squeeze(handles.freqPW{BP}(nChDisplay,1,:));
bar(tmpx,tmpy);
axis ([handles.distFrom(BP) handles.distTo(BP) 0 inf], 'auto y');

hold on
amp = handles.mainModeAmp{BP}(nChDisplay,1);
mn = handles.mainModeMean{BP}(nChDisplay,1);
sd = handles.mainModeStd{BP}(nChDisplay,1);
secMode = handles.secondaryMode{BP}(nChDisplay,1,:);

mainModeGauss = amp * exp(-(tmpx-mn).^2/(2*(sd^2)));
plot (tmpx,mainModeGauss,'red','LineWidth',2);
plot (tmpx,squeeze(secMode(1,1,:)),'green','LineWidth',2);

% Fill the table with the proper data corresponding to selected channel,
% trial and band pass.

for i=1:handles.fileInNum
    handles.dtaPWlocalSummary (i, 1:15) = handles.dtaPWsummary(i, nChDisplay, 1, BP, 1:15); 
end

set(handles.analyzePWdata.distributionLogPowerTable,'Data',handles.dtaPWlocalSummary);
zebraguidata(hObject,handles);

end

