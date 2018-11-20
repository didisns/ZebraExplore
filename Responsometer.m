function varargout = Responsometer(varargin)
% RESPONSOMETER MATLAB code for Responsometer.fig
%      RESPONSOMETER, by itself, creates a new RESPONSOMETER or raises the existing
%      singleton*.
%
%      H = RESPONSOMETER returns the handle to a new RESPONSOMETER or the handle to
%      the existing singleton*.
%
%      RESPONSOMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESPONSOMETER.M with the given input arguments.
%
%      RESPONSOMETER('Property','Value',...) creates a new RESPONSOMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Responsometer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Responsometer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Responsometer

% Last Modified by GUIDE v2.5 03-Oct-2018 11:40:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Responsometer_OpeningFcn, ...
                   'gui_OutputFcn',  @Responsometer_OutputFcn, ...
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


% --- Executes just before Responsometer is made visible.
function Responsometer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Responsometer (see VARARGIN)

% Choose default command line output for Responsometer
handles.output = hObject;

handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

handles.baselines = {}; %baselines will be a cell array;
                        %row1 = "semgments" calculated by computeMeanOnTrials.m
                        %row2 = folder name
                        %one column for each imported "baseline" file
                        
handles.stimuli = {};   %stimuli will be a cell array;
                        %row1 = "semgments" calculated by computeMeanOnTrials.m
                        %row2 = folder name
                        %row3 = "templates" calculated by computeMeanOnTrials.m
                        %one column for each imported "baseline" file


% Update handles structure
zebraguidata(hObject, handles);

% UIWAIT makes Responsometer wait for user response (see UIRESUME)
% uiwait(handles.Responsometer);


% --- Outputs from this function are returned to the command line.
function varargout = Responsometer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in eraseAll_bt.
function eraseAll_bt_Callback(hObject, eventdata, handles)
% hObject    handle to eraseAll_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Calc_bt.
function Calc_bt_Callback(hObject, eventdata, handles)
% hObject    handle to Calc_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function criterion_txt_Callback(hObject, eventdata, handles)
% hObject    handle to criterion_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of criterion_txt as text
%        str2double(get(hObject,'String')) returns contents of criterion_txt as a double


% --- Executes during object creation, after setting all properties.
function criterion_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to criterion_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stimulus_ppp.
function stimulus_ppp_Callback(hObject, eventdata, handles)
% hObject    handle to stimulus_ppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stimulus_ppp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stimulus_ppp


% --- Executes during object creation, after setting all properties.
function stimulus_ppp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimulus_ppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stimulusImport_bt.
function stimulusImport_bt_Callback(hObject, eventdata, handles)
% hObject    handle to stimulusImport_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.stimuli)
    handles.stimulus_sldr.Max = handles.parentHandles.nCh;
    handles.stimulus_sldr.Value = 1;
    handles.stimulus_sldr.Min = 1;
    handles.stimulus_sldr.SliderStep = [round(1/(handles.stimulus_sldr.Max-1)), round(1/(handles.stimulus_sldr.Max-1))];  
    handles.stimulusChannel_txt.String = ['Ch 1/' num2str(handles.stimulus_sldr.Max)];
else
    if handles.stimulus_sldr.Max ~= handles.parentHandles.nCh
        error('stimuli with different number of channels')
    end
end
handles.stimuli{1,end+1} = handles.parentHandles.segment;
handles.stimuli{3,end+1} = handles.parentHandles.template;

[~, handles.stimuli{2,end}, ~] = fileparts(handles.parentHandles.dir_in(1:end-1));

handles.stimulus_ppp.String = handles.baselines(2,:)';
handles.stimulus_ppp.Value = length(handles.stimulus_ppp.String);

zebraguidata(hObject,handles)

% --- Executes on button press in stimulusErase_bt.
function stimulusErase_bt_Callback(hObject, eventdata, handles)
% hObject    handle to stimulusErase_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to baselineErase_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
toBeDel = handles.stimulus_ppp.Value;
handles.stimuli(:,toBeDel) = [];
if handles.stimulus_ppp.Value == size(handles.stimuli,2)+1
    handles.stimulus_ppp.Value = handles.stimulus_ppp.Value-1;
end
handles.stimulus_ppp.String = handles.stimuli(2,:)';
zebraguidata(hObject,handles)

% --- Executes on selection change in baselineList_lbx.
function baselineList_lbx_Callback(hObject, eventdata, handles)
% hObject    handle to baselineList_lbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns baselineList_lbx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from baselineList_lbx


% --- Executes during object creation, after setting all properties.
function baselineList_lbx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baselineList_lbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in baselineImport_bt.
function baselineImport_bt_Callback(hObject, eventdata, handles)
% hObject    handle to baselineImport_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isempty(handles.baselines)
    handles.baseline_sldr.Max = handles.parentHandles.nCh;
    handles.baseline_sldr.Value = 1;
    handles.baseline_sldr.Min = 1;
    handles.baseline_sldr.SliderStep = [round(1/(handles.baseline_sldr.Max-1)), round(1/(handles.baseline_sldr.Max-1))];  
    handles.baselineChannel_txt.String = ['Ch 1/' num2str(handles.baseline_sldr.Max)];
else
    if handles.baseline_sldr.Max ~= handles.parentHandles.nCh
        error('baselines with different number of channels')
    end
end
handles.baselines{1,end+1} = handles.parentHandles.segment; %baselines is a matrix with format: channels * trials * time * multistimulus ;  

[~, handles.baselines{2,end}, ~] = fileparts(handles.parentHandles.dir_in(1:end-1));

handles.baselineList_lbx.String = handles.baselines(2,:)';
handles.baselineList_lbx.Value = length(handles.baselineList_lbx.String);

zebraguidata(hObject,handles)

% --- Executes on button press in baselineErase_bt.
function baselineErase_bt_Callback(hObject, eventdata, handles)
% hObject    handle to baselineErase_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
toBeDel = handles.baselineList_lbx.Value;
handles.baselines(:,toBeDel) = [];
if handles.baselineList_lbx.Value == size(handles.baselines,2)+1
    handles.baselineList_lbx.Value = handles.baselineList_lbx.Value-1;
end
handles.baselineList_lbx.String = handles.baselines(2,:)';
zebraguidata(hObject,handles)

% --- Executes on slider movement.
function stimulus_sldr_Callback(hObject, eventdata, handles)
% hObject    handle to stimulus_sldr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function stimulus_sldr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimulus_sldr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function baseline_sldr_Callback(hObject, eventdata, handles)
% hObject    handle to baseline_sldr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.baseline_sldr.Value = round(handles.baseline_sldr.Value);
handles.baselineChannel_txt.String = ['Ch ' num2str(handles.baseline_sldr.Value) '/' num2str(handles.baseline_sldr.Max)];


% --- Executes during object creation, after setting all properties.
function baseline_sldr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baseline_sldr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
