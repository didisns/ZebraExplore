function varargout = multiChannelAnalyzer(varargin)
% MULTICHANNELANALYZER MATLAB code for multiChannelAnalyzer.fig
%      MULTICHANNELANALYZER, by itself, creates a new MULTICHANNELANALYZER or raises the existing
%      singleton*.
%
%      H = MULTICHANNELANALYZER returns the handle to a new MULTICHANNELANALYZER or the handle to
%      the existing singleton*.
%
%      MULTICHANNELANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTICHANNELANALYZER.M with the given input arguments.
%
%      MULTICHANNELANALYZER('Property','Value',...) creates a new MULTICHANNELANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multiChannelAnalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multiChannelAnalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multiChannelAnalyzer

% Last Modified by GUIDE v2.5 31-May-2018 12:58:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multiChannelAnalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @multiChannelAnalyzer_OutputFcn, ...
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


% --- Executes just before multiChannelAnalyzer is made visible.
function multiChannelAnalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multiChannelAnalyzer (see VARARGIN)

% Choose default command line output for multiChannelAnalyzer
handles.output = hObject;

% establish the link with the calling process
handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

% Update handles structure
zebraguidata(hObject, handles);

% UIWAIT makes multiChannelAnalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multiChannelAnalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function sliderChA_Callback(hObject, eventdata, handles)

function sliderChA_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderChB_Callback(hObject, eventdata, handles)

function sliderChB_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function fromTxt_Callback(hObject, eventdata, handles)

function fromTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function toTxt_Callback(hObject, eventdata, handles)

function toTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function wndTxt_Callback(hObject, eventdata, handles)

function wndTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ovlTxt_Callback(hObject, eventdata, handles)

function ovlTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function XcorrBtn_Callback(hObject, eventdata, handles)
% Compute the Xcorrelation between the selected channels.

% refresh parentHandles
handles.parentHandles = zebraguidata(handles.parentObject);

% Initialize data
lag = [];
handles.corr = [];
handles.lagTime = [];
zebraguidata(hObject,handles);

ch(1) = get(handles.sliderChA,'Value');
ch(2) = get(handles.sliderChB,'Value');
%trial =
fromTime = str2double(get(handles.fromTxt,'String'));
toTime = str2double(get(handles.toTxt,'String'));
sp = handles.parentHandles.sp;
interval(1) = int32(fromTime/sp +1);
if interval(1) < 1, interval(1) = 1;
end
interval(2) = int32(toTime/sp +1);
if (interval(2) > handles.parentHandles.dtaLen), interval(2) = handles.parentHandles.dtaLen;
end
corrWind = str2double(get(handles.corrWindow,'String'));
corrWind = int32(corrWind/handles.parentHandles.sp);
[LFP1, LFP2] = obtainData(handles,ch,interval);

[handles.corr, lag] = xcorr(LFP1,LFP2,corrWind);
handles.lagTime = lag * handles.parentHandles.sp;
zebraguidata(hObject, handles)
plotXcorr  (hObject, handles)
plotLFP(hObject, handles)

function plotXcorr (hObject, handles)
axes (handles.xCorrSpectra)
l(1) = str2double(get(handles.lagMin,'String'));
l(2) = str2double(get(handles.lagMax,'String'));

% Find out what is the index interval associated to this frequency
% interval. Parse the handles.f array and extract the indexes (ind) included int he interval.
lag = handles.lagTime;
lag(lag<l(1) | lag>l(2)) = 0;
ind = find(lag);
indLeft = ind(1);
indRight = ind(end);

plot(handles.lagTime(indLeft:indRight),handles.corr(indLeft:indRight))
axis ([l(1) l(2) -inf inf]);
set(gca,'fontsize',8)

function plotLFP (hObject, handles)

% plot the LFP signal under analysis in the LFP axes. The plot is limited
% to the analysis interval defined by the from, to fields.

ch(1) = get(handles.sliderChA,'Value');
ch(2) = get(handles.sliderChB,'Value');

fromTime = str2double(get(handles.fromTxt,'String'));
toTime = str2double(get(handles.toTxt,'String'));
sp = handles.parentHandles.sp;
interval(1) = int32(fromTime/sp +1);
if interval(1) < 1, interval(1) = 1;
    fromTime = 0;
end
interval(2) = int32(toTime/sp +1);
if (interval(2) > handles.parentHandles.dtaLen), interval(2) = handles.parentHandles.dtaLen;
    toTime = sp*(handles.parentHandles.dtaLen-1);
end
[LFP1, LFP2] = obtainData(handles,ch,interval);
time = fromTime:handles.parentHandles.sp:toTime;
axes(handles.LFP)
plot(time,LFP1,time,LFP2)
axis ([fromTime toTime -inf inf]);
set(gca,'fontsize',8)

function [LFP1 LFP2]=obtainData(handles, ch, interval)
switch get(handles.dataSource,'selectedObject');
    case handles.workLFP
        LFP1 = handles.parentHandles.workLFP(ch(1),1,interval(1):interval(2)); 
        LFP2 = handles.parentHandles.workLFP(ch(2),1,interval(1):interval(2)); 
        LFP1 = permute(LFP1,[3 2 1]);
        LFP2 = permute(LFP2,[3 2 1]);
    case handles.deleakedLFP
        LFP1 = handles.parentHandles.deleakedLFP(ch(1),1,interval(1):interval(2)); 
        LFP2 = handles.parentHandles.deleakedLFP(ch(2),1,interval(1):interval(2));         
        LFP1 = permute(LFP1,[4 3 2 1]);
        LFP2 = permute(LFP2,[4 3 2 1]);
    case handles.HP
        LFP1 = handles.parentHandles.bandPassed_LFP(5,ch(1),1,interval(1):interval(2));
        LFP2 = handles.parentHandles.bandPassed_LFP(5,ch(2),1,interval(1):interval(2));
        LFP1 = permute(LFP1,[4 3 2 1]);
        LFP2 = permute(LFP2,[4 3 2 1]);
    case handles.BP1
        LFP1 = handles.parentHandles.bandPassed_LFP(2,ch(1),1,interval(1):interval(2));
        LFP2 = handles.parentHandles.bandPassed_LFP(2,ch(2),1,interval(1):interval(2));
        LFP1 = permute(LFP1,[4 3 2 1]);
        LFP2 = permute(LFP2,[4 3 2 1]);
    case handles.BP2
        LFP1 = handles.parentHandles.bandPassed_LFP(3,ch(1),1,interval(1):interval(2));
        LFP2 = handles.parentHandles.bandPassed_LFP(3,ch(2),1,interval(1):interval(2));
        LFP1 = permute(LFP1,[4 3 2 1]);
        LFP2 = permute(LFP2,[4 3 2 1]);
    case handles.BP3
        LFP1 = handles.parentHandles.bandPassed_LFP(4,ch(1),1,interval(1):interval(2));
        LFP2 = handles.parentHandles.bandPassed_LFP(4,ch(2),1,interval(1):interval(2));
        LFP1 = permute(LFP1,[4 3 2 1]);
        LFP2 = permute(LFP2,[4 3 2 1]);
end        
        
function lagMin_Callback(hObject, eventdata, handles)
plotXcorr (hObject, handles)

function lagMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lagMax_Callback(hObject, eventdata, handles)
plotXcorr (hObject, handles)

function lagMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function corrWindow_Callback(hObject, eventdata, handles)

function corrWindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function shortTimeXcorr_Callback(hObject, eventdata, handles)

% compute the short time cross correlation


ch(1) = get(handles.sliderChA,'Value');
ch(2) = get(handles.sliderChB,'Value');

wnd = str2double(get(handles.wndTxt,'String'));
step = str2double(get(handles.ovlTxt,'String'));
% now compute the number of windows
sp = handles.parentHandles.sp;
wndI = int32(wnd/sp);
total = handles.parentHandles.dtaLen;
wndPnt = int32(wnd/sp);
stepPnt = int32(step/sp);
midWndBegin = ceil(wndPnt/2);
midWndEnd = total - midWndBegin;
midPoints = midWndBegin:stepPnt:midWndEnd;
nwnd = length (midPoints);
% preallocate the output
tmpCorr = zeros(2*wnd+1,nwnd);
% now compute the CC for each window
% pre allocate the STCC temporary array
tmpCorr = zeros(2*wndI+1,nwnd);
for i=1:nwnd
    time = wnd/2 + (i-1)*step;
    interval(1) = int32(midPoints(i)-midWndBegin);
    interval(2) = int32(midPoints(i)+midWndBegin);
    if interval(1) == 0, interval(1) = 1;
    end
    [LFP1, LFP2] = obtainData(handles,ch,interval);
    [handles.corr, lag] = xcorr(LFP1,LFP2,wndI);
    tmpCorr(:,i) = handles.corr(:);
    tmpTime(i) = time;
end    
handles.lagTime = lag * handles.parentHandles.sp;
handles.STxCorr = tmpCorr;
axes(handles.shortTimeXcorr)
imagesc(handles.STxCorr)
zebraguidata(hObject, handles)
 

function crossSpectrum_Callback(hObject, eventdata, handles)
% Compute the cross spectrum between the selected channels.

% refresh parentHandles
handles.parentHandles = zebraguidata(handles.parentObject);

% Initialize data
lag = [];
handles.corr = [];
handles.lagTime = [];
zebraguidata(hObject,handles);

ch(1) = get(handles.sliderChA,'Value');
ch(2) = get(handles.sliderChB,'Value');
%trial =
fromTime = str2double(get(handles.fromTxt,'String'));
toTime = str2double(get(handles.toTxt,'String'));
sp = handles.parentHandles.sp;
interval(1) = int32(fromTime/sp +1);
if interval(1) < 1, interval(1) = 1;
end
interval(2) = int32(toTime/sp +1);
if (interval(2) > handles.parentHandles.dtaLen), interval(2) = handles.parentHandles.dtaLen;
end
corrWind = str2double(get(handles.corrWindow,'String'));
corrWind = int32(corrWind/handles.parentHandles.sp);
% I am leving the possibility of computing the Xspectum on a band passed
% data. This might be a good idea to make small peaks more visible in
% presence of a massive peak, such as the one caused by slow wave activity.

[LFP1, LFP2] = obtainData(handles,ch,interval);

% prepare the parameters for Chronux
params.tapers = [5 9];
params.Fs = 1/handles.parentHandles.sp;
params.tpass = [0 100];
params.err = [];
[handles.Coher,phi,handles.xS12,S1,S2,handles.f] = coherencyc (LFP1,LFP2,params);

zebraguidata(hObject, handles)
plotXspectrum  (hObject, handles)
plotLFP(hObject, handles)

function plotXspectrum (hObject, handles)
axes (handles.crossSpectra)
l(1) = str2double(get(handles.freqMin,'String'));
l(2) = str2double(get(handles.freqMax,'String'));

% Find out what is the index interval associated to this frequency
% interval. Parse the handles.f array and extract the indexes (ind) included int he interval.
freq = handles.f;
freq(freq<l(1) | freq>l(2)) = 0;
ind = find(freq);
indLeft = ind(1);
indRight = ind(end);

%mx = max(handles.xS12(indLeft:indRight));
plot(handles.f(indLeft:indRight),handles.xS12(indLeft:indRight),'LineWidth',1)
axis ([l(1) l(2) -inf inf]);
set(gca,'fontsize',8)
set(gca, 'XScale', 'log')

function shortTimeXspextrum_Callback(hObject, eventdata, handles)

function freqMin_Callback(hObject, eventdata, handles)
plotXspectrum (hObject, handles)

function freqMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freqMax_Callback(hObject, eventdata, handles)
plotXspectrum (hObject, handles)

function freqMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function holdOnxCorr_Callback(hObject, eventdata, handles)
if get(hObject, 'Value'),
    % checkbox is on
    hold all (handles.xCorrSpectra)
else
    hold off (handles.xCorrSpectra)
end    
cla(handles.xCorrSpectra)
plotXcorr (hObject, handles)

function holdXspectra_Callback(hObject, eventdata, handles)
if get(hObject, 'Value'),
    % checkbox is on
    hold all (handles.crossSpectra)
else
    hold off (handles.crossSpectra)
end
cla(handles.crossSpectra)
plotXspectrum (hObject, handles)
