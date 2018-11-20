function varargout = ZebraSays(varargin)
% ZEBRASAYS MATLAB code for ZebraSays.fig
%      ZEBRASAYS, by itself, creates a new ZEBRASAYS or raises the existing
%      singleton*.
%
%      H = ZEBRASAYS returns the handle to a new ZEBRASAYS or the handle to
%      the existing singleton*.
%
%      ZEBRASAYS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZEBRASAYS.M with the given input arguments.
%
%      ZEBRASAYS('Property','Value',...) creates a new ZEBRASAYS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ZebraSays_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ZebraSays_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ZebraSays

% Last Modified by GUIDE v2.5 24-May-2017 16:41:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ZebraSays_OpeningFcn, ...
                   'gui_OutputFcn',  @ZebraSays_OutputFcn, ...
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


% --- Executes just before ZebraSays is made visible.
function ZebraSays_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ZebraSays (see VARARGIN)

% Choose default command line output for ZebraSays
handles.output = hObject;

% recover the handle to ZebraExplore
handles.zebraMain = findobj('Tag','ZebraMainFig');
if ~isempty(handles.zebraMain)
    % get handles and other user-defined data associated to Gui1
    handles.zebraMainData = zebraguidata(handles.zebraMain);
 
    % test of access of ZebraMainData
    % x = getappdata(handles.zebraMain,'lp')
    % maybe you want to set the text in Gui2 with that from Gui1
    % set(handles.text1,'String',get(handles.zebraMainData.edit1,'String'));
    % maybe you want to get some data that was saved to the Gui1 app
    %x = getappdata(h,'x');

end
 
% Update handles structure
zebraguidata(hObject, handles);

% UIWAIT makes ZebraSays wait for user response (see UIRESUME)
% uiwait(handles.ZebraSays);


% --- Outputs from this function are returned to the command line.
function varargout = ZebraSays_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in clearOutable.
function clearOutable_Callback(hObject, eventdata, handles)

% This function clear all elements from the table

handles = zebraguidata(hObject);    % this to refresh the local values of handles
handles.zebraMainData = zebraguidata(handles.zebraMain); % load the local copy of the parent's handles

% x = getappdata(handles.zebraMain,'lp')
% x = handles.zebraMainData.sp
handles.zebraMainData.cnt = 0;      % Zero the number of data in the table (handles.cnt,11)
handles.zebraMainData.dtaSummary = {'' '' '' 0 0 0 0 0 0 0 ''};
set(handles.outable,'Data',handles.zebraMainData.dtaSummary);
zebraguidata(handles.zebraMain, handles.zebraMainData);
zebraguidata(hObject, handles);

% --- Executes when entered data in editable cell(s) in outEventsSummary.
function outEventsSummary_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to outEventsSummary (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
