function varargout = trialRejection(varargin)
% TRIALREJECTION MATLAB code for trialRejection.fig
%      TRIALREJECTION, by itself, creates a new TRIALREJECTION or raises the existing
%      singleton*.
%
%      H = TRIALREJECTION returns the handle to a new TRIALREJECTION or the handle to
%      the existing singleton*.
%
%      TRIALREJECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIALREJECTION.M with the given input arguments.
%
%      TRIALREJECTION('Property','Value',...) creates a new TRIALREJECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trialRejection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trialRejection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trialRejection

% Last Modified by GUIDE v2.5 21-May-2018 15:48:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trialRejection_OpeningFcn, ...
                   'gui_OutputFcn',  @trialRejection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


% End initialization code - DO NOT EDIT


% --- Executes just before trialRejection is made visible.
function trialRejection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trialRejection (see VARARGIN)

% Choose default command line output for trialRejection
handles.output = hObject;

%use synchroFunction now to get zebraMain handles
synchroFunction(hObject,handles)
handles=zebraguidata(hObject);

%initialize these variables with the default value in the GUI
handles.bins=str2double(handles.bins_edit.String);
handles.segm=str2double(handles.segm_edit.String);
handles.iterFraction = str2double(handles.iterationsFrac_edit.String);
handles.iterN = str2double(handles.iterationsN_edit.String);

%initialize handles.mode (the work modality) with the default selection in
%the GUI
mode_rb_SelectionChangedFcn(hObject,[], handles)
handles=zebraguidata(hObject);

%BPcheckA works as a logical vector; it stores the selected BP for the Analysis; LP index is 1, BPn index is n+1,
%HP index is 5, no-filter index is 6
handles.BPcheckA=zeros(1,6);
BPa_bg_SelectionChangedFcn(hObject, [], handles);
handles=zebraguidata(hObject);

handles.BPcheckE=zeros(1,6);
BPe_bg_SelectionChangedFcn(hObject, [], handles);
handles=zebraguidata(hObject);

% Update handles structure
zebraguidata(hObject, handles);

% UIWAIT makes trialRejection wait for user response (see UIRESUME)
% uiwait(handles.trialRejection);


% --- Outputs from this function are returned to the command line.
function varargout = trialRejection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in compPlot_btn.
function compPlot_btn_Callback(hObject, eventdata, handles)
% hObject    handle to compPlot_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%newLFP is the signal we work on here. It is imported from zebraMain
%according to the selected BandPass.
handles.newLFP=[];
selBP=find(handles.BPcheckA);
if ~(selBP==6)
    handles.newLFP(:,:)=handles.zebraMainData.bandPassed_LFP(selBP,handles.currentCh,:,:);
else
    handles.newLFP(:,:)=handles.zebraMainData.LFP(handles.currentCh,:,:);
end

if handles.nTr == 1
    handles.newLFP = handles.newLFP' ;
end
    
%dissect each trial trace in "segm" segments; the obtained segmLFP variable
%is a matrix in the form trial*segment*time
handles.segmLFP=[];
for i=1:handles.nTr
    for j=1:handles.segm
        startI=floor(handles.trLen/handles.segm)*(j-1)+1;
        endI=floor(startI+handles.trLen/handles.segm)-1;
        %endI=floor(startI+handles.trLen/handles.segm);
        %handles.segmLFP(i,j,:)= handles.zebraMainData.LFP(handles.currentCh,i,startI:endI);
        handles.segmLFP(i,j,:)= handles.newLFP(i,startI:endI);
    end
end

segmDim=size(handles.segmLFP);
handles.segmLen=segmDim(3);


%different modes make different computation!
switch handles.mode
    case 'e'
        computeNrgDistr(hObject,handles)
    
    case 'ht'
        compute_ht(hObject,handles)
        
        handles=zebraguidata(hObject);
        plot_result(hObject,handles)
        
    case 'manual'
        compute_manual(hObject,handles)
        
        handles=zebraguidata(hObject);
        plot_result(hObject,handles)
    case 'ssr'
        compute_ssr(hObject,handles)
end

function compute_ssr(hObject,handles)
%SSR = sum of the squared residuals. Actually, the rms of the residuals is
%calculated.

handles.segmOn=1;

if handles.discFrac_rb.Value
    handles.nIter = floor(handles.nTr/handles.iterFraction);
else
    handles.nIter = handles.iterN;
end

%lavoro sempre su newLFP ma uso una matrice di logicals.
handles.TrOk = true(handles.nTr,handles.nIter);
%handles.newLFP = handles.newLFP - mean(handles.newLFP,2);
handles.newLFP = squeeze(handles.segmLFP(:,1,:));

i=1;
% figure %devo decidere dove plottare
% hold on
% title('medians rss: i-th vs  (i-1)-th')

handles.med = []; %is a structure iterations*time, whose i-th row is the median of the trials at the i-th iteration.
handles.res = zeros(handles.nIter,handles.nTr,length(handles.newLFP)); %matrix of the residuals. its dimensioins are iteration*trial*time. 
handles.ssr = zeros(handles.nIter,handles.nTr); %matrix of the rms of the residuals (ssr). its dimensioins are iteration*trial.
handles.MaxMed = []; %stores the max rss for every iteration.
%medSsr = []; %stores the rss between a median and the previous one.


while i <= handles.nIter
    handles.med(i,:) = median(handles.newLFP(handles.TrOk(:,i),:),1);
    tmp = handles.newLFP(handles.TrOk(:,i),:);
    dim=size(tmp);
    handles.res(i,:,:) = handles.newLFP - handles.med(i,:);
    %handles.ssr(i,:) = rms(squeeze(handles.res(i,:,:)),2); %ssr for each trial
    handles.ssr(i,:) = sum(squeeze(handles.res(i,:,:)).^2,2); %ssr per provare
    handles.res(i,~handles.TrOk(:,i),:)=NaN;
    handles.ssr(i,~handles.TrOk(:,i))=NaN;
    [handles.MaxMed(i),iMm] = nanmax(handles.ssr(i,:)); %Save the highest ssr between a trial and the median, and the index of that trial (to be eliminated).

    handles.TrOk(:,i+1) = handles.TrOk(:,i);
    handles.TrOk(iMm,i+1) = false;

%     % plot the rss between a median and the previous one    
%     if i>1
%         medSsr(i) = rms(med(i,:)-med(i-1,:)); %ssr between a median and the previous one
%         %medSsr(i) = sum((med(i,:)-med(i-1,:)).^2);
%         plot(i,medSsr(i),'bo');
%         drawnow
%         if i>2
%             plot([i-1 i],[medSsr(i-1) medSsr(i)],'b-')
%         end
%     end
    
    i=i+1;
end

handles.med(i,:) = median(handles.newLFP(handles.TrOk(:,i),:));

handles.legal = handles.TrOk(:,end);

plot_ssr(hObject,handles);
plot_result(hObject,handles);

function plot_ssr(hObject,handles)

ax = handles.rejAxes_2; %lower axes
axes(ax)
cla reset % the option "reset" clear both the main and the secondary axis
hold on

plot((1:length(handles.MaxMed)),handles.MaxMed,'bo-')
% plot((1:length(handles.MaxMed)),handles.MaxMed(2:end),'b-')
if ~handles.legal(handles.currentTr)
    %ritrovare a quale iterazione è stato eliminato il trial
    i=1;
    while handles.TrOk(handles.currentTr,i)
        i=i+1;
    end %se all'iterazione i è risultato escluso, era il MaxMed dell'iterazione i-1
    plot(i-1,handles.MaxMed(i-1),'rs','MarkerFaceColor','r')
end

xlabel('iteration')
ylabel('max ssr')

islegal(hObject,handles)

zebraguidata(hObject,handles);



function compute_manual(hObject,handles)

handles.segmOn = 1:max(handles.segm);

handles.legal = ones(1,handles.nTr); %logical array; handles.legal(n)=1 if the n-th trial is accepted

leg = sum(handles.legal); %leg = number of legal trials

%update the panel on the GUI
handles.n_legal_txt.String = num2str(leg);
handles.n_rej_txt.String = num2str(handles.nTr - leg);

zebraguidata(hObject,handles);

plot_threshold(hObject,handles);

function compute_ht(hObject,handles)

handles.threshold = str2double(handles.set_thresh_edit.String);

handles.segmOn=[];
if handles.trialsON_ck.Value
    handles.segmOn = 1:str2double(handles.trialsON_edit.String);
else
    handles.segmOn = 1:max(handles.segm);
end

handles.legal = ones(1,handles.nTr); %logical array; handles.legal(n)=1 if the n-th trial is accepted

for i=1:handles.nTr
    for j=1:max(handles.segmOn)
    %     if max(abs(handles.newLFP(i,:)-mean(handles.newLFP(i,:)))) > thr  %check if the maximum deviation from the mean of the trial exceeds the threshold.
        if max(abs(handles.segmLFP(i,j,:))) > handles.threshold %normal threshold
            handles.legal(i)=0;
        end
    end
end

% for i=1:handles.nTr
%     
%     %     if max(abs(handles.newLFP(i,:)-mean(handles.newLFP(i,:)))) > thr  %check if the maximum deviation from the mean of the trial exceeds the threshold.
%         if max(abs(handles.segmLFP(i,handles.segmOn,:))) < handles.threshold %normal threshold
%             handles.legal(i)=1;
%         end
%     
% end

zebraguidata(hObject,handles);

plot_threshold(hObject,handles);

function plot_threshold(hObject,handles)
ax = handles.rejAxes_2; %lower axes
axes(ax)
cla reset % the option "reset" clear both the main and the secondary axis
hold on

if ~handles.BPcheckA(6) 
    y(:,:)=handles.newLFP(handles.currentTr,:);
else
    y(:,:)=handles.zebraMainData.LFP(handles.currentCh,handles.currentTr,:);
end
x=(0:handles.trLen-1)*handles.zebraMainData.sp;
plot(x,y);

%this plot function is called both in 'ht' mode and in 'manual' mode, but
%the threshold lines must be plotted only in 'ht' mode.
if strcmp( handles.mode, 'ht' ) 
    yL=ylim;
    if handles.threshold > yL(2)-10
        yL(2) = handles.threshold + 10;
    end
    
    if -handles.threshold < yL(1)+10
        yL(1) = -handles.threshold - 10;
    end
    xt=xlim;
    plot([xt(1) handles.segmLen*handles.segmOn(end)*handles.zebraMainData.sp],[handles.threshold handles.threshold],'r')
    plot([xt(1) handles.segmLen*handles.segmOn(end)*handles.zebraMainData.sp],[-handles.threshold -handles.threshold],'r')
    
    ylim(yL); %the reason of this mess is that i don't want the threshold lines to be plotted on the edges of the axes. 
end
ylabel('uV')
xlabel('')

islegal(hObject,handles) %update the status of the current trial on the GUI

zebraguidata(hObject,handles)

function islegal(hObject,handles)
if handles.legal(handles.currentTr)
    handles.trial_status_txt.String='ACCEPTED';
    handles.trial_status_txt.ForegroundColor='g';
    handles.discard_bt.String='Reject trial';
else
    handles.trial_status_txt.String='REJECTED';
    handles.trial_status_txt.ForegroundColor='r';
    handles.discard_bt.String='Accept trial';
end
leg = sum(handles.legal);
handles.n_legal_txt.String = num2str(leg);
handles.n_rej_txt.String = num2str(handles.nTr - leg);
%zebraguidata(hObject,handles); non dovrebbe essercene bisogno
    
% --- Executes on button press in discard_bt.
function discard_bt_Callback(hObject, eventdata, handles)
% hObject    handle to discard_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.legal(handles.currentTr) = ~handles.legal(handles.currentTr);
zebraguidata(hObject,handles);

islegal(hObject,handles)

plot_result(hObject,handles)

% handles=zebraguidata(hObject);
% zebraguidata(hObject,handles);

function plot_result(hObject,handles)
ax = handles.rejAxes_histo;
axes(ax)
cla
hold on

title('Cleaned trace average')

handles.meanLFP=[];
% handles.meanLFP = mean(handles.newLFP(logical(handles.legal),:),1);
exportBP = find(handles.BPcheckE);
if exportBP == 6
    if sum(handles.legal)==1
        handles.meanLFP = mean(squeeze(handles.zebraMainData.LFP(handles.currentCh,logical(handles.legal),:)),2)';
    else
        handles.meanLFP = mean(squeeze(handles.zebraMainData.LFP(handles.currentCh,logical(handles.legal),:)),1);
    end
else
    if sum(handles.legal)==1
        handles.meanLFP = mean(squeeze(handles.zebraMainData.bandPassed_LFP(exportBP,handles.currentCh,logical(handles.legal),:)),2)';
    else
        handles.meanLFP = mean(squeeze(handles.zebraMainData.bandPassed_LFP(exportBP,handles.currentCh,logical(handles.legal),:)),1);
    end
end

x = (1:length(handles.meanLFP))*handles.zebraMainData.sp;
plot(x, handles.meanLFP, 'LineWidth', 0.5)
ylim auto
yL=ylim;
plot(x(1:handles.segmLen*handles.segmOn(end)), handles.meanLFP(1:handles.segmLen*handles.segmOn(end)),'b', 'LineWidth', 2)

% 
% handles.stdLFP = std (handles.newLFP(logical(handles.legal),:),1);
% xx = [x fliplr(x)];
% yy = [handles.meanLFP+handles.stdLFP fliplr(handles.meanLFP-handles.stdLFP)];
% poli = patch (xx,yy,'b');
% poli.EdgeAlpha = 0;
% poli.FaceAlpha = 0.1;
ylim(yL);
zebraguidata(hObject,handles);
ylabel('uV')
xlabel('')
set(handles.plot1_limMax, 'String', num2str(yL(2)));
set(handles.plot1_limMin, 'String', num2str(yL(1)));
% if handles.minMax_rb.Value
%     handles.plotMode='mm';
%     computeMinMax(hObject,handles)
% end
% if handles.nrgDistr_rb.Value
%     handles.plotMode='e';
%     computeNrgDistr(hObject,handles)
% end

function segm_edit_Callback(hObject, eventdata, handles)
% hObject    handle to segm_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segm_edit as text
%        str2double(get(hObject,'String')) returns contents of segm_edit as a double
handles.segm=str2double(handles.segm_edit.String);
zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function segm_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segm_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bins_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bins_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bins_edit as text
%        str2double(get(hObject,'String')) returns contents of bins_edit as a double
handles.bins=str2double(handles.bins_edit.String);
zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function bins_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bins_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function computeMinMax(hObject,handles)
%currentTr=handles.zebraMainData.currentTrial;

handles.maxV=max(handles.zebraMainData.LFP(handles.currentCh,:,:),[],3);
handles.minV=min(handles.zebraMainData.LFP(handles.currentCh,:,:),[],3);
handles.maxDelta=handles.maxV-handles.minV;

%histogram(log10(maxDelta1),20)
% histogram(maxDelta,20)
fromV=[];
toV=[];
plotHist(hObject,handles,handles.maxDelta,handles.bins) %fromE and toE will be implemented later

% handles.plotMode='mm'; %CAMBIA!!
% params.yLim=[];
% params.xLim=[];
plot_mm(hObject,handles);

zebraguidata(hObject,handles);

function plot_mm(hObject,handles)

% axes(handles.rejAxes_2)
% To clear the active side, use cla. To clear both sides of the axes and remove the right y-axis, use cla reset.

cla reset
hold on

%define x and y
x=handles.minV;
y=handles.maxV;
Plot2=plot(x,y);
Plot2.LineStyle='none';
Plot2.Marker='o';

xlabel('minV (uV)')
ylabel('maxV (uV)')

%to be implemented maybe later
if ~isempty(handles.xLim)
    xlim(handles.xLim)
end
if ~isempty(handles.yLim)
    ylim(handles.yLim)
end

% --- Executes on button press in SynchZebra_bt.
function SynchZebra_bt_Callback(hObject, eventdata, handles)
% hObject    handle to SynchZebra_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
synchroFunction(hObject,handles)

function synchroFunction(hObject,handles)

handles.zebraMain=findobj('Tag', 'ZebraMainFig');
if ~isempty(handles.zebraMain)
    handles.zebraMainData=zebraguidata(handles.zebraMain);
end

if isfield(handles.zebraMainData, 'dtaComplete')
    handles.trLen=handles.zebraMainData.dtaLen;
    handles.trLen_txt.String=['Trial lenght = ' num2str(handles.trLen) ' points'];
    
    handles.nTr=handles.zebraMainData.nTrials;

    if (handles.zebraMainData.nCh>1)
        handles.ch_slider.Max=handles.zebraMainData.nCh;
        handles.ch_slider.Min=1;
        handles.ch_slider.SliderStep=[1/(handles.zebraMainData.nCh-1) 1/(handles.zebraMainData.nCh-1)];
        handles.ch_slider.Value=1;
        handles.ch_slider.Visible='on';
        handles.ch_txt.Visible='on';    
        handles.ch_txt.String=['Ch ' num2str(handles.zebraMainData.currentCh)];
    else
        set (handles.ch_slider,'Visible','off');
        set (handles.ch_txt,'Visible','off');
    end

    if (handles.nTr>1)
        handles.tr_slider.Max=handles.nTr;
        handles.tr_slider.Min=1;
        handles.tr_slider.SliderStep=[1/(handles.nTr-1) 1/(handles.nTr-1)];
        handles.tr_slider.Value=1;
        handles.tr_slider.Visible='on';
        handles.trial_txt.Visible='on';    
        handles.trial_txt.String=['Trial ' num2str(handles.zebraMainData.currentTrial)];
    else
        set (handles.ch_slider,'Visible','off');
        set (handles.ch_txt,'Visible','off');
    end


    handles.currentCh=handles.zebraMainData.currentCh;
    handles.currentTr=handles.zebraMainData.currentTrial;
    
    tmpCk=handles.zebraMainData.BPcomputed;
    if ~tmpCk(1)
        handles.LP_rb.Enable='off';
        handles.LPe_rb.Enable='off';
    else
        handles.LP_rb.Enable='on';
        handles.LPe_rb.Enable='on';
    end
    if ~tmpCk(2)
        handles.BP1_rb.Enable='off';
        handles.BP1e_rb.Enable='off';
    else
        handles.BP1_rb.Enable='on';
        handles.BP1e_rb.Enable='on';
    end
    if ~tmpCk(3)
        handles.BP2_rb.Enable='off';
        handles.BP2e_rb.Enable='off';
    else
        handles.BP2_rb.Enable='on';
        handles.BP2e_rb.Enable='on';
    end
    if ~tmpCk(4)
        handles.BP3_rb.Enable='off';
        handles.BP3e_rb.Enable='off';
    else
        handles.BP3_rb.Enable='on';
        handles.BP3e_rb.Enable='on';
    end
    if ~tmpCk(5)
        handles.HP_rb.Enable='off';
        handles.HPe_rb.Enable='off';
    else
        handles.HP_rb.Enable='on';
        handles.HPe_rb.Enable='on';
    end
    
end
zebraguidata(hObject,handles);


% --- Executes on slider movement.
function ch_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.currentCh=handles.ch_slider.Value;
handles.ch_txt.String=['Ch ' num2str(handles.currentCh)];
zebraguidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function ch_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function computeNrgDistr(hObject,handles)

%calculate the energy of each segment -> segmNrg is a matrix in the form of
%trial*segment
handles.segmNrg = rms(handles.segmLFP,3);
%handles.segmNrg = rms(handles.newLFP,3);
if handles.trialsON_ck.Value
    handles.segmOn = 1:str2double(handles.trialsON_edit.String);
%     handles.segmNrg = handles.segmNrg(:,handles.segmOn);
%     handles.currentHist = handles.segmNrg;
else
%     handles.currentHist = handles.segmNrg;
    dim=size(handles.segmNrg);
    handles.segmOn = 1:dim(2);
end
handles.segmNrg = handles.segmNrg(:,handles.segmOn);
handles.currentHist = handles.segmNrg;
%plotHist alla fine


segmCentres=(handles.segmLen+1)/2:handles.segmLen:(handles.segmLen+1)/2+handles.segmLen*(handles.segm-1);
if handles.segm>1
    segmSepar=handles.segmLen+0.5:handles.segmLen:(handles.segmLen+0.5)+handles.segmLen*(handles.segm-2);
    handles.segmSeparT=(segmSepar-1)*handles.zebraMainData.sp;
end
handles.timeNrg=(segmCentres-1)*handles.zebraMainData.sp;

%gab prova
handles.timeNrg=handles.timeNrg(handles.segmOn);
% handles.plotMode='e'; %CAMBIA!

%plot_e alla fine
plotHist(hObject,handles,handles.segmNrg,handles.bins)
handles=zebraguidata(hObject);

plot_e(hObject,handles)%qui dovrò far plottare LFP filtrato o meno


function plotHist (hObject,handles,x,bins)
axes(handles.rejAxes_histo)
cla
flagLog=handles.histoLog_ck.Value;
if flagLog
    %histogram(log10(x),bins)
    x=10*log10(x);
else
    %histogram(x,bins)
end
histogram(x,bins)
% handles.rejAxes_histo.UserData=x;
ylabel('Absolute frequency')
if flagLog
    xlabel('power (dB)')
else
    xlabel('power')
end

ylim auto
title('power distribution')

handles.fromHist_edit.String = num2str( floor(min(min(x))) );
handles.toHist_edit.String = num2str( ceil(max(max(x))) );
selection_bt_Callback(hObject,[],handles)


% --- Executes on button press in histoLog_ck.
function histoLog_ck_Callback(hObject, eventdata, handles)
% hObject    handle to histoLog_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% x=handles.rejAxes_histo.UserData;

if isfield(handles,'currentHist')
    x=handles.currentHist;
    
    zebraguidata(hObject,handles);

    plotHist(hObject,handles,x,handles.bins)
    plot_e (hObject,handles)
end
%zebraguidata(hObject,handles);


function plot_e (hObject,handles) %qui dovrò far plottare LFP filtrato o meno, a seconda del fiter_ck

axes(handles.rejAxes_2)
% To clear the active side, use cla. To clear both sides of the axes and remove the right y-axis, use cla reset.
cla reset
hold on
% yyaxis left

if ~handles.BPcheckA(6)
    y(:,:)=handles.newLFP(handles.currentTr,:);
else
    y(:,:)=handles.zebraMainData.LFP(handles.currentCh,handles.currentTr,:);
end
x=(0:handles.trLen-1)*handles.zebraMainData.sp;

Plot2=plot(x,y);
%         plot(x,abs(y));    %per controllo: è molto simile alla rms ovviamente 
ylabel('uV')

for i=1:length(handles.segmSeparT)
    plot([handles.segmSeparT(i) handles.segmSeparT(i)],ylim, '--k')
end

yyaxis right
nrg=handles.segmNrg(handles.currentTr,:);
if handles.histoLog_ck.Value
    plot(handles.timeNrg,10*log10(nrg),'o-')
    ylabel('RMS (dB)')
else
    plot(handles.timeNrg,nrg,'o-')
    ylabel('RMS')
end


xlim([0,x(end)])

islegal(hObject,handles)

zebraguidata(hObject,handles)




% --- Executes on slider movement.
function tr_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.currentTr=handles.tr_slider.Value;
handles.trial_txt.String=['Trial ' num2str(handles.currentTr)];
%compPlot_btn_Callback(hObject,[],handles)

switch handles.mode
    case 'e'
        plot_e(hObject,handles)
    
    case 'ht'
        plot_threshold(hObject,handles)
        %plot_result(hObject,handles)
        
    case 'manual'
        plot_threshold(hObject,handles)
        plot_result(hObject,handles)
        
    case 'ssr'
        plot_ssr(hObject,handles)
end

%zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function tr_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function fromHist_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fromHist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function fromHist_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fromHist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function toHist_edit_Callback(hObject, eventdata, handles)
% hObject    handle to toHist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function toHist_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toHist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selection_bt.
function selection_bt_Callback(hObject, eventdata, handles)
% hObject    handle to selection_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flagLog=handles.histoLog_ck.Value;
fromX = round(str2double(handles.fromHist_edit.String),3);
toX = round(str2double(handles.toHist_edit.String),3);
% checkMin=[];
% checkMax=[];
% for i=1:handles.zebraMainData.nTr
%     checkMin=handles.zebraMainData.LFP(handles.currentCh,i,:)>fromX;
%     checkMax=handles.zebraMainData.LFP(handles.currentCh,i,:)<toX;
% end
% checkTr=checkMin.*checkMax;

%find trials to be excluded, the ones in checkR vector.
if flagLog
    [checkR checkC]=find(10*log10(handles.segmNrg)<fromX | 10*log10(handles.segmNrg)>toX);
else
    [checkR checkC]=find(handles.segmNrg<fromX | handles.segmNrg>toX);
end
checkR=unique(checkR);

handles.legal=ones(1,handles.nTr);
handles.legal(checkR)=0;

zebraguidata(hObject,handles);

islegal(hObject,handles);
% if handles.filter_ck.Value
%     temp(:,:)=handles.newLFP;
% else
%     temp(:,:)=handles.zebraMainData.LFP(handles.currentCh,:,:);
% end
% temp(checkR,:)=[];
% dimTemp=size(temp);
% ciao=1;
% fid=fopen([handles.zebraMainData.dtaComplete(1:end-4) '_clean.asc'],'w+');
% fprintf(fid,'% s',handles.zebraMainData.expID);
% fprintf(fid,'%s','    ');
% for i=1:dimTemp(1)
%     fprintf(fid,'% d',2048);
%     fprintf(fid,'%s','    ');
% end
% for i=1:handles.trLen
%     tempDim=size(temp);
%     for j=1:tempDim(1)
%             fprintf(fid,'% f',temp(j,i));
%             fprintf(fid,'%s','    ');
%     end
% end
% fclose(fid);



% --- Executes on button press in save_bt.
function save_bt_Callback(hObject, eventdata, handles)
% hObject    handle to save_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if handles.filter_ck.Value
%     temp(:,:)=handles.newLFP;
% else
%     temp(:,:)=handles.zebraMainData.LFP(handles.currentCh,:,:);
% end
exportBP = find(handles.BPcheckE);
if exportBP == 6
    temp(:,:)=handles.zebraMainData.LFP(handles.currentCh,:,:);
else
    temp(:,:)=handles.zebraMainData.bandPassed_LFP(exportBP,handles.currentCh,:,:);
end

excluded=~handles.legal;
temp(excluded,:)=[];
dimTemp = size(temp);
if (dimTemp(2)==1 && dimTemp(1)>dimTemp(2)) %i have a single trial; in this case it appears as a column vector
    temp = temp';
    dimTemp = size(temp);
end
fid=fopen([handles.zebraMainData.dtaComplete(1:end-4) '_cleanWholeTrace.asc'],'w+');

%the following procedure works only for MeyerVEPs
fprintf(fid,'% s',handles.zebraMainData.expID);
fprintf(fid,'%s','    ');
for i=1:dimTemp(1)
    fprintf(fid,'% d',2048);
    fprintf(fid,'%s','    ');
end
for i=1:handles.trLen
    %tempDim=size(temp);
    %for j=1:tempDim(1)
    for j=1:dimTemp(1)
            fprintf(fid,'% f',temp(j,i));
            fprintf(fid,'%s','    ');
    end
end
fclose(fid);


% --- Executes on button press in trialsON_ck.
function trialsON_ck_Callback(hObject, eventdata, handles)
% hObject    handle to trialsON_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trialsON_ck



function trialsON_edit_Callback(hObject, eventdata, handles)
% hObject    handle to trialsON_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trialsON_edit as text
%        str2double(get(hObject,'String')) returns contents of trialsON_edit as a double


% --- Executes during object creation, after setting all properties.
function trialsON_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialsON_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_ck_Callback(hObject, eventdata, handles)
% 
% if hObject.Value
%     handles.noFilt_rb.Enable='off';
% else
%     handles.noFilt_rb.Enable='on';
% end
% zebraguidata(hObject,handles);


function BPa_bg_SelectionChangedFcn(hObject, eventdata, handles)


%handles.BPcheckA stores the selected BP; LP index is 1, BPn index is n+1,
%HP index is 5, noFilt index is 6

if handles.LP_rb.Value
    handles.BPcheckA(1)=1;
else
    handles.BPcheckA(1)=0;
end

if handles.BP1_rb.Value
    handles.BPcheckA(2)=1;
else
    handles.BPcheckA(2)=0;
end

if handles.BP2_rb.Value
    handles.BPcheckA(3)=1;
else
    handles.BPcheckA(3)=0;
end

if handles.BP3_rb.Value
    handles.BPcheckA(4)=1;
else
    handles.BPcheckA(4)=0;
end

if handles.HP_rb.Value
    handles.BPcheckA(5)=1;
else
    handles.BPcheckA(5)=0;
end

if ~handles.noFilt_rb.Value
%     handles.filter_ck.Enable='on';
    handles.BPcheckA(6)=0;
else
%     handles.filter_ck.Enable='off';
    handles.BPcheckA(6)=1;
end
zebraguidata(hObject,handles);

% --- Executes on button press in filter_ck.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to filter_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of filter_ck


function mode_rb_SelectionChangedFcn(hObject, eventdata, handles)

%rejAxes2

if handles.minMax_rb.Value
    handles.mode='mm';
end
if handles.nrgDistr_rb.Value
    handles.mode='e';
    handles.discard_bt.Enable='off';
    handles.selection_bt.Enable='on';
end
if handles.ht_rb.Value
    handles.mode='ht';
    handles.discard_bt.Enable='off';
    handles.selection_bt.Enable='off';
end
if handles.manual_rb.Value
    handles.mode='manual';
    handles.discard_bt.Enable='on';
    handles.selection_bt.Enable='off';
end
if handles.ssr_rb.Value
    handles.mode='ssr';
    handles.discard_bt.Enable='off';
    handles.selection_bt.Enable='off';
    handles.trialsON_ck.Value=0;
    %handles.trialsON_ck.Enable='off';
    handles.segm_edit.String='1';
    segm_edit_Callback(hObject,[],handles);
    handles=zebraguidata(hObject);
    handles.gif_bt.Enable='on';
else
    handles.gif_bt.Enable='off';
    %handles.trialsON_ck.Enable='on';
end

zebraguidata(hObject,handles);



function set_thresh_edit_Callback(hObject, eventdata, handles)
1;



function set_thresh_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in result_bt.
function result_bt_Callback(hObject, eventdata, handles)
% hObject    handle to result_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function plot1_limMax_Callback(hObject, eventdata, handles)

ax = handles.rejAxes_histo;
axes(ax)

yL=ylim;
yL(2) = str2double(hObject.String);
ylim(yL);


function plot1_limMax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot1_limMin_Callback(hObject, eventdata, handles)

ax = handles.rejAxes_histo;
axes(ax)

yL=ylim;
yL(1) = str2double(hObject.String);
ylim(yL);


function plot1_limMin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function saveMeanMat_bt_Callback(hObject, eventdata, handles)

tmp=handles.meanLFP;
save([handles.zebraMainData.dtaComplete(1:end-4) '_Mean.mat'], 'tmp');



function iterationsFrac_edit_Callback(hObject, eventdata, handles)

handles.iterFraction = str2double(hObject.String);
zebraguidata(hObject,handles);


function iterationsFrac_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function iterationsN_edit_Callback(hObject, eventdata, handles)

handles.iterN = str2double(hObject.String);
zebraguidata(hObject,handles);


function iterationsN_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in BPe_bg.
function BPe_bg_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in BPe_bg 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.BPcheckA stores the selected BP; LP index is 1, BPn index is n+1,
%HP index is 5, noFilt index is 6

if handles.LPe_rb.Value
    handles.BPcheckE(1)=1;
else
    handles.BPcheckE(1)=0;
end

if handles.BP1e_rb.Value
    handles.BPcheckE(2)=1;
else
    handles.BPcheckE(2)=0;
end

if handles.BP2e_rb.Value
    handles.BPcheckE(3)=1;
else
    handles.BPcheckE(3)=0;
end

if handles.BP3e_rb.Value
    handles.BPcheckE(4)=1;
else
    handles.BPcheckE(4)=0;
end

if handles.HPe_rb.Value
    handles.BPcheckE(5)=1;
else
    handles.BPcheckE(5)=0;
end

if ~handles.noFilte_rb.Value
%     handles.filter_ck.Enable='on';
    handles.BPcheckE(6)=0;
else
%     handles.filter_ck.Enable='off';
    handles.BPcheckE(6)=1;
end
zebraguidata(hObject,handles);


% --- Executes on button press in gif_bt.
function gif_bt_Callback(hObject, eventdata, handles)
% hObject    handle to gif_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%progressively plot all the means after rejection using medians; then save a gif
h = figure;
hold on

filename = [handles.zebraMainData.dtaComplete(1:end-4) '_cleaning.gif'];

exportBP = find(handles.BPcheckE);

dimMed=size(handles.med);
for i=1:dimMed(1)
    if exportBP == 6
        tmp = mean(squeeze(handles.zebraMainData.LFP(handles.currentCh,logical(handles.TrOk(:,i)),:)),1);
    else
        tmp = mean(squeeze(handles.zebraMainData.bandPassed_LFP(exportBP,handles.currentCh,logical(handles.TrOk(:,i)),:)),1);
    end
%     tmp = mean(handles.newLFP(handles.TrOk(:,i),:));
    plot((0:2047)*handles.zebraMainData.sp,tmp)
    ylim([-15,10])
    title(['progress, i=' num2str(i)])
    drawnow
    
    pause(0.2)
    
    
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end

if exportBP == 6
    tmp = mean(squeeze(handles.zebraMainData.LFP(handles.currentCh,logical(handles.TrOk(:,i)),:)),1);
else
    tmp = mean(squeeze(handles.zebraMainData.bandPassed_LFP(exportBP,handles.currentCh,logical(handles.TrOk(:,i)),:)),1);
end

plot((0:2047)*handles.zebraMainData.sp,tmp,'k','lineWidth',2)

frame = getframe(h);
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',10); 
