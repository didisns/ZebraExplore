function varargout = importTRCfilesOptions(varargin)
% IMPORTTRCFILESOPTIONS MATLAB code for importTRCfilesOptions.fig
%      IMPORTTRCFILESOPTIONS, by itself, creates a new IMPORTTRCFILESOPTIONS or raises the existing
%      singleton*.
%
%      H = IMPORTTRCFILESOPTIONS returns the handle to a new IMPORTTRCFILESOPTIONS or the handle to
%      the existing singleton*.
%
%      IMPORTTRCFILESOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPORTTRCFILESOPTIONS.M with the given input arguments.
%
%      IMPORTTRCFILESOPTIONS('Property','Value',...) creates a new IMPORTTRCFILESOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before importTRCfilesOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to importTRCfilesOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help importTRCfilesOptions

% Last Modified by GUIDE v2.5 03-Sep-2017 12:17:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @importTRCfilesOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @importTRCfilesOptions_OutputFcn, ...
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


% --- Executes just before importTRCfilesOptions is made visible.
function importTRCfilesOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to importTRCfilesOptions (see VARARGIN)

% Choose default command line output for importTRCfilesOptions
handles.output = hObject;
handles.dirIn = varargin{1};              
handles.fileIn = varargin{2};              
handles.completePath = varargin{3};              
% Update handles structure
zebraguidata(hObject, handles);
set(handles.dirTxt,'String',['Dir: ' handles.dirIn]);
set(handles.fileTxt,'String',['File name: ' handles.fileIn]);

% UIWAIT makes importTRCfilesOptions wait for user response (see UIRESUME)
uiwait(handles.ImportTRCselectOptions);

% --- Outputs from this function are returned to the command line.
function varargout = importTRCfilesOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.dataOut;

%varargout = num2cell(handles.dataOut.srate);
varargout = struct2cell(handles.dataOut);
close(hObject);     % close the modal dialog figure.

function triggerCk_Callback(hObject, eventdata, handles)

function channelsIn_Callback(hObject, eventdata, handles)

function channelsIn_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function triggerChStr_Callback(hObject, eventdata, handles)

function triggerChStr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ImportTRC_Callback(hObject, eventdata, handles)

%   Detailed explanation goes here

%   >> [DATAOUT]=readtrc(SETTINGS)
%
%
% INPUT
%   SETTINGS is a struct holding the parameters for reading the .TRC file
%   SETTINGS has the following fields:
%
%       SETTINGS.filename :                 Name of file to be imported
%       SETTINGS.loadevents.state :         'yes' for loading event triggers
%                                           'no' for not
%       SETTINGS.loadevents.type :          'marker' for event triggers inserted 
%                                           on 'MKR channel
%                                           'eegchan' for triggers inserted on 
%                                           EEG channels
%                                           'none' or
%                                           'both'
%       SETTINGS.loadevents.dig_ch1:        number of name of eegchan marker channel 1
%       SETTINGS.loadevents.dig_ch1_label:  label to give events on eegchan marker channel 1
%       SETTINGS.loadevents.dig_ch2:        number of name of eegchan marker channel 2
%       SETTINGS.loadevents.dig_ch2_label:  label to give events on eegchan marker channel 2
%       SETTINGS.chan_adjust_status:        1 for adjusting amp of channels 0 for not
%       SETTINGS.chans                      channels to load, [ ] for all
%       (default)
%       SETTINGS.chan_adjust                channels to adjust
%
%   Alternant method: enter 11 input arguments each corresponding to one of the above
%   fields.  If only one arg is used,  it must be the struct above.  If more, there
%   must be 11 inputs in the order above i.e. OUT=readtrc(filename,eventstate....etc).
%

%dirIn = 'C:\Users\Gix\Documents\Work\Science\Shank3\Meyer\Exp July 2017 onward\2017 07 04\FEDERICO\';
%fileIn = 'EEG_256.TRC';
fileInComplete = [handles.dirIn handles.fileIn];

TRCparam.filename = fileInComplete;
triggerFlag = get(handles.triggerCk,'Value');
if triggerFlag,
    TRCparam.loadevents.state = 'yes';          % 'yes' for loading event triggers
    TRCparam.loadevents.type = 'eegchan';       % 'marker' for event triggers on MKR channel  
                                                % 'eegchan' for triggers inserted on EEG channels
else
    TRCparam.loadevents.state = 'no';           % 'yes' for loading event triggers
end

TRCparam.loadevents.dig_ch1 = get(handles.triggerChStr, 'String');          % channel number where eegchan marker are saved
TRCparam.loadevents.dig_ch1_label = 'stim';                                % label to give events on eegchan marker channel 1
TRCparam.loadevents.dig_ch2 = '';            % number of name of eegchan marker channel 2
TRCparam.loadevents.dig_ch2_label = [];      % label to give events on eegchan marker channel 2
TRCparam.chan_adjust_status = 0;             % 1 for adjusting amp of channels 0 for not

TRCparam.chans = get(handles.channelsIn,'String'); 
refChannel = get(handles.referenceStr,'String');
if (isnan(str2double(refChannel))),
    % no reference channel selected. Do nothing but remember later on
else
    % add the reference channel as first element in the input string list
    TRCparam.chans = [refChannel ' ' TRCparam.chans];
end    

%TRCparam.chans = ['11 19 20'];             % channels to load, [ ] for all
TRCparam.chan_adjust = [];                  % channels to adjust

%[DATAOUT] = readtrc(TRCparam);
[handles.dataOut] = readtrcGMR(TRCparam);

% OUTPUT
%   DATAOUT with same fields as EEGLAB EEG file structure. (see eeg_checkset.m)
%   The following fields are use: 
%       DATAOUT.data 
%       DATAOUT.filename
%       DATAOUT.filepath
%       DATAOUT.srate
%       DATAOUT.setname
%       DATAOUT.pnts
%       DATAOUT.nbchan
%       DATAOUT.trials
%       DATAOUT.xmin
%       DATAOUT.ref
%       TRC.event(j).latency    index of event j

set(handles.subjectTxt,'String',['Subject name: ' handles.dataOut.subject]);
set(handles.dateTxt,'String',['Date: ' handles.dataOut.comments]);

set(handles.sfTxt,'String',['Sampling frequency (Hz): ' num2str(handles.dataOut.srate)]);
handles.dataLength = handles.dataOut.pnts/handles.dataOut.srate;
set(handles.dataDurationTxt,'String',['Data length (s): ' num2str(handles.dataLength)]);
handles.nEvents = length(handles.dataOut.event);
set(handles.eventsTxt,'String',['Events: ' num2str(handles.nEvents)]);
nDataIn = size(handles.dataOut.data);

% subtract the reference channel if defined.
if ~(isnan(str2double(refChannel))),
    % Subtract the reference channel (first element of the data) from the
    % other channels
    for i=2:nDataIn(1)
        handles.dataOut.data (i,:) = handles.dataOut.data (i,:) - handles.dataOut.data (1,:);
    end
end    


% x = 0:(dataOut.pnts-1);
% x = x/handles.dataOut.srate;
% figure
% plot (x, handles.dataOut.data);
zebraguidata(hObject, handles);


% --- Executes on button press in endBtn.
function endBtn_Callback(hObject, eventdata, handles)
% hObject    handle to endBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume (gcbf);
%close(gcbf);    % gcbf is the hndle of the figure that contains the command



function referenceStr_Callback(hObject, eventdata, handles)
% hObject    handle to referenceStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of referenceStr as text
%        str2double(get(hObject,'String')) returns contents of referenceStr as a double


% --- Executes during object creation, after setting all properties.
function referenceStr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to referenceStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in subtractRefCk.
function subtractRefCk_Callback(hObject, eventdata, handles)
% hObject    handle to subtractRefCk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subtractRefCk
