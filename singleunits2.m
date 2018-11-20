function varargout = singleunits2(varargin)
% SINGLEUNITS2 MATLAB code for singleunits2.fig
%      SINGLEUNITS2, by itself, creates a new SINGLEUNITS2 or raises the existing
%      singleton*.
%
%      H = SINGLEUNITS2 returns the handle to a new SINGLEUNITS2 or the handle to
%      the existing singleton*.
%
%      SINGLEUNITS2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLEUNITS2.M with the given input arguments.
%
%      SINGLEUNITS2('Property','Value',...) creates a new SINGLEUNITS2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before singleunits2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to singleunits2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help singleunits2

% Last Modified by GUIDE v2.5 16-Nov-2018 16:31:04

% Comment by Didi: Notes on how to use the program.
% 1) Single Units extraction was done according to Quiroga et al. In order
% to stick to the recommendations in the paper, use high pass filtered data
% of 300 Hz set in the ZebraExplore main program.

% 2) Set the channel and sampling period yourself: the single unit program
% uses the channel and sampling period as set by the user here, NOT the
% ones set in the main ZebraExplore figure

% 3) This program uses results from the EventDetection program. To use
% this, you must have computed in EventDetection first, and then in the
% ZebraExplore main figure press Update Handles EventDetection. Note that
% all information used here from EventDetection is the start and end time
% of up and down states. Depending on your goal, you can therefore select
% on which channel you want to identify up and down states by setting the
% channel in EventDetection. It is recommended use the channel with the clearest up states to
% detect up and down state times in EventDetection, and then keep these times
% for both channel calculations in the single units program

% 4) The order of execution is important: always click calculate in the SU
% box on top first, only then you can plot and/or calculate in up and down
% states

% 5) To save the timing of positive and negative single units, as well as
% of up and downs states, click the save times button, and make sure the
% lpt number is set correctly is you used a sum (more than 1) of recordings
% at the ephys setup

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @singleunits2_OpeningFcn, ...
                   'gui_OutputFcn',  @singleunits2_OutputFcn, ...
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


% --- Executes just before singleunits2 is made visible.
function singleunits2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to singleunits2 (see VARARGIN)

% Choose default command line output for singleunits2
handles.output = hObject;
% recover the handle to ZebraExplore
handles.zebraMain = findobj('Tag','ZebraMainFig');
if ~isempty(handles.zebraMain)
    % get handles and other user-defined data associated to Gui1
    handles.zebraMainData = zebraguidata(handles.zebraMain);
 
    % test of access of ZebraMainData
    % x = getappdata(handles.zebraMain,'lp')
    % maybe you want to set the text in Gui2 with that from Gui1
    % set(handles.text1,'String',get(handles.zebraMainData.threshold,'String'));
    % maybe you want to get some data that was saved to the Gui1 app
    %x = getappdata(h,'x');

end
%handles.zebraMainData.evDetection = eventDetection(hObject,handles);
handles.zebraMainData.evDetData = zebraguidata (handles.zebraMainData.evDetection);
% Update handles structure
axes(handles.axes10)
title('')
axes(handles.axes11)
title('')
zebraguidata(hObject, handles);

% UIWAIT makes singleunits2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function findpeaks(hObject, eventdata, handles)
% function to identify number and timing of positively and negatively
% deflecting peaks

% let's start by creating vectors contaning the time points and the high
% pass filtered data
% first obtain the current handles
handles = zebraguidata(hObject);
channelnumber = str2double(get(handles.chnumber, 'String'));
LFP_hp_temp = handles.zebraMainData.bandPassed_LFP(5,channelnumber,:,:);    
% note: LFP has as size of: 1 1 1 900000 so use only the last column, as
% done in the column below
handles.LFP_hp = LFP_hp_temp(:,:)'; % we now have the high pass filtered data in a single column
% Now find the threshold, as proposed by Quiroga et al.
std = median(abs(handles.LFP_hp)/0.6745);
thr = 4*std;    
   
% identify the peaks higher than threshold
handles.x = 1:length(handles.LFP_hp);
len = length(handles.LFP_hp);
% setting the to be calculated values/matrixes to find positive peaks
handles.nmax = 0; 
handles.maxAt = []; 

% we can now look at episodes above threshold
for i = 1:len % loop through the ratio matrix
    if i == 1 && handles.LFP_hp(i) >= thr % for the first time point, if it is higher than threshold it is the first peak
        handles.nmax = handles.nmax + 1; %increase the peak number by one
        handles.maxAt(handles.nmax) = handles.x(i); % use this time as the start time of the current peak
    elseif i > 1 && handles.LFP_hp(i) >= thr && handles.LFP_hp(i-1) < thr % we test that it is the first time point higher than threshold
        handles.nmax = handles.nmax + 1; % then we count it as a new peak
        handles.maxAt(handles.nmax) = handles.x(i); % and label its start time
    end
end

%we now do the same for downward peaks
% setting the to be calculated values/matrixes to find negative peaks
handles.nmin = 0; 
handles.minAt = []; 

for i = 1:len
    if i == 1 && handles.LFP_hp(i) <= -thr 
        handles.nmin = handles.nmin + 1; 
        handles.minAt(handles.nmin) = handles.x(i); 
    elseif i > 1 && handles.LFP_hp(i) <= -thr && handles.LFP_hp(i-1) > -thr 
        handles.nmin = handles.nmin + 1; 
        handles.minAt(handles.nmin) = handles.x(i);
    end
end
% update the handle structure
zebraguidata(hObject, handles);
    
function plotLFP(hObject, eventdata, handles)
% function to plot the LFP

% first obtain the current handles
handles = zebraguidata(hObject);

% Obtain  the LFP vector: 
channelnumber = str2double(get(handles.chnumber, 'String'));
LFPreal_temp = handles.zebraMainData.LFP(channelnumber,1,:);
LFPreal = LFPreal_temp(:,:)';

% Find the Y axes:
%min and max values of the LFP
LFPreal_min = min(LFPreal); LFPreal_max = max(LFPreal);
if LFPreal_min < 0
    LFP_min_temp = ceil(10*abs(LFPreal_min));
    y_min = -LFP_min_temp/10;
else LFP_min_temp = ceil(10*LFPreal_min)
    y_min = LFP_min_temp/10;
end       
if LFPreal_max < 0
    LFP_max_temp = ceil(10*abs(LFPreal_max));
    y_max = -LFP_max_temp/10;
else LFP_max_temp = ceil(10*LFPreal_max);
    y_max = LFP_max_temp/10;
end
    
% Now find the X axes:
% convert x axis to time
acqF = str2double(get(handles.edit_acqF, 'string'));
handles.time = handles.x./acqF;
if isempty(get(handles.edit_tstart_plot, 'String')) % if the user did not enter the start time, i.e. the edit text box is empty
    handles.x_min_time = handles.time(1); % the start time in seconds
    handles.x_min_index = handles.x(1); % index of the start time
else handles.x_min_time = str2double(get(handles.edit_tstart_plot, 'String')); % if the user entered the start time, find the entered number
    if handles.x_min_time == 0 
        handles.x_min_index = handles.x(1); % if it is 0, we use as index the first one of axes
    else
        handles.x_min_index = handles.x_min_time*acqF; % otherwhise we can simply convert the time to index
    end
end
% the same for the end time
if isempty(get(handles.edit_tend_plot, 'String'));
    handles.x_max_time = handles.time(end);
    handles.x_max_index = handles.x(end);
else handles.x_max_time = str2double(get(handles.edit_tend_plot, 'String'));
    handles.x_max_index = handles.x_max_time*acqF;
end

% if end time is smaller than start time
if handles.x_max_time <= handles.x_min_time
    error('end time has to be larger than start time')
end
    
% Now plot the LFP: 
axes(handles.axes10);    
xlim([handles.x_min_time handles.x_max_time]); ylim([y_min y_max]);    
plot(handles.time(handles.x_min_index:handles.x_max_index), LFPreal(handles.x_min_index:handles.x_max_index), 'k');
% end by updating handles structure
zebraguidata(hObject, handles);

function plotSUmax(hObject, eventdata, handles)
% function to plot the positive 'peaks' (single units)
handles = zebraguidata(hObject);  
axes(handles.axes11);
% use the same x limits as created for the LFP
xlim([handles.x_min_time handles.x_max_time]);    
y_su_max = ismember(handles.x, handles.maxAt); % creates a vector of the same length of x, with values 0 or 1 if the same value 
% is a member of maxAt, in other words if there is a spike.
plot(handles.time(handles.x_min_index:handles.x_max_index), y_su_max(handles.x_min_index:handles.x_max_index), 'r');
% end by updating handles structure
zebraguidata(hObject, handles);

function plotSUmin(hObject, eventdata, handles)
% function to plot the negative 'peaks' (single units)
handles = zebraguidata(hObject);  
axes(handles.axes11);
xlim([handles.x_min_time handles.x_max_time]);    
y_su_min = ismember(handles.x, handles.minAt);
plot(handles.time(handles.x_min_index:handles.x_max_index), y_su_min(handles.x_min_index:handles.x_max_index), 'b');
% end by updating handles structure
zebraguidata(hObject, handles);

function plotSUall(hObject, eventdata, handles)
% function to plot all peaks, positive ones in red, negative in blue
handles = zebraguidata(hObject);
axes(handles.axes11);
xlim([handles.x_min_time handles.x_max_time]); 
y_su_min = ismember(handles.x, handles.minAt);
y_su_max = ismember(handles.x, handles.maxAt);
plot(handles.time(handles.x_min_index:handles.x_max_index), y_su_min(handles.x_min_index:handles.x_max_index), 'b',...
    handles.time(handles.x_min_index:handles.x_max_index), y_su_max(handles.x_min_index:handles.x_max_index), 'r');
% end by updating handles structure
zebraguidata(hObject, handles);

function Extract_SWA(hObject, eventdata, handles)
% function to extract the indexes of up and down states

% first obtained the start and end points of the up states obtained with
% EventDetection. Choose in EventDetection from which channel you want to
% identify up states.

% the fromT toT values in the upstates list of evDetection is longer than
% NUS and NDS (this happens during fusing/removal up states)..these last
% fromT and toT values contains the last up state time and so if I would
% include these I would count several single units several times. To
% prevent this, I first find the correct value of number of up and down
% states and then extract only fromT and toT of this length
nus = handles.zebraMainData.evDetData.NUS;
nds = handles.zebraMainData.evDetData.NDS;
us_start = (handles.zebraMainData.evDetData.upstates.fromT(1:nus))';
us_end = (handles.zebraMainData.evDetData.upstates.toT(1:nus))';
ds_start = (handles.zebraMainData.evDetData.downstates.fromT(1:nds))';
ds_end = (handles.zebraMainData.evDetData.downstates.toT(1:nds))';
% test if each start point also has an end point, i.e. if the data are
% correct
if length(us_start) ~= length(us_end) || length(ds_start) ~= length(ds_end)
    error('input data from EvDet do not have correct lengths')
end

% Then transform the time points to index point
acqF = str2double(get(handles.edit_acqF, 'string'));
us_startI = us_start.*acqF;
us_endI = us_end.*acqF;
ds_startI = ds_start.*acqF;
ds_endI = ds_end*acqF;

% add to the handles structure
handles.us_startI = us_startI;
handles.us_endI = us_endI;
handles.ds_startI = ds_startI;
handles.ds_endI = ds_endI;
handles.nus = nus;
handles.nds = nds;

% end by updating the handles structure
zebraguidata(hObject, handles);

function Extract_SUpos_SWA (hObject, eventdata, handles)
% find the number of positive single units located within up states and down states
% first set the needed values
nus = handles.nus;
len = length(handles.maxAt);
su_in_us = 0; % this is going to count the number of single units in upstates
% now we will first loop through every up state
for us = 1:nus
    idx1 = handles.us_startI(us); % and for each up state find the start..
    idx2 = handles.us_endI(us); % and end point
    % then we can loop through the indexes of su
    for sui = 1:len
        if handles.maxAt(sui) > idx1 && handles.maxAt(sui) < idx2
            su_in_us = su_in_us + 1;
        end 
        % add an if statement to stop when you have reached higher than the
        % current up state, just to save time
        if handles.maxAt(sui) > idx2
            break
        end
    end
end

% Now let's do the same for the down states
nds = handles.nds;
su_in_ds = 0;
for ds = 1:nds
    idx1 = handles.ds_startI(ds); 
    idx2 = handles.ds_endI(ds);    
    for sui = 1:len
        if handles.maxAt(sui) > idx1 && handles.maxAt(sui) < idx2
            su_in_ds = su_in_ds + 1;
        end 
        if handles.maxAt(sui) > idx2
            break
        end
    end    
end

% test if correct
%if su_in_us + su_in_ds > handles.nmax
%    error ('single units in upstates + down is larger than total single units: something wrong in the code')
%end

% add to the handles structure and update handles structure
handles.pos_su_in_us = su_in_us;
handles.pos_su_in_ds = su_in_ds;
zebraguidata(hObject, handles);

function Extract_SUneg_SWA (hObject, eventdata, handles)
% find the number of negative single units located within up states and down states
% works the same as the function for the positive single units
nus = handles.nus;
len = length(handles.minAt);
su_in_us = 0;
for us = 1:nus
    idx1 = handles.us_startI(us); 
    idx2 = handles.us_endI(us);    
    for sui = 1:len
        if handles.minAt(sui) > idx1 && handles.minAt(sui) < idx2
            su_in_us = su_in_us + 1;
        end
        if handles.minAt(sui) > idx2
            break
        end
    end    
end

nds = handles.nds;
su_in_ds = 0;
for ds = 1:nds
    idx1 = handles.ds_startI(ds); 
    idx2 = handles.ds_endI(ds);    
    for sui = 1:len
        if handles.minAt(sui) > idx1 && handles.minAt(sui) < idx2
            su_in_ds = su_in_ds + 1;
        end
        if handles.minAt(sui) > idx2
            break
        end
    end    
end

%if su_in_us + su_in_ds > handles.nmin
%    error ('single units in upstates + down is larger than total single units: something wrong in the code')
%end

handles.neg_su_in_us = su_in_us;
handles.neg_su_in_ds = su_in_ds;
zebraguidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = singleunits2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function chnumber_Callback(hObject, eventdata, handles)
% hObject    handle to chnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chnumber as text
%        str2double(get(hObject,'String')) returns contents of chnumber as a double


% --- Executes during object creation, after setting all properties.
function chnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculate_SU.
function calculate_SU_Callback(hObject, eventdata, handles)
% hObject    handle to calculate_SU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% function will execute the findpeaks function and display counted positive
% and negative peaks
if get(hObject, 'Value')
    findpeaks(hObject,handles);
    % update handle structure now that the function has been executed and
    % new variables have been added to the handles variable
    handles = zebraguidata(hObject);    
    % Show the total number of upward going and downward going peaks 
    set(handles.edit_nmax, 'string', num2str(handles.nmax));
    set(handles.edit_nmin, 'string', num2str(handles.nmin));   
 
else
    % if checkbox is unclicked, empty the textboxes
    set(handles.edit_nmax, 'string', ''), set(handles.edit_nmin, 'string', '');   
end
handles = zebraguidata(hObject); 
% Hint: get(hObject,'Value') returns toggle state of calculate_SU



function edit_nmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nmax as text
%        str2double(get(hObject,'String')) returns contents of edit_nmax as a double


% --- Executes during object creation, after setting all properties.
function edit_nmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nmin as text
%        str2double(get(hObject,'String')) returns contents of edit_nmin as a double


% --- Executes during object creation, after setting all properties.
function edit_nmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plot.
function checkbox_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% plot all the single units
if get(hObject, 'Value')
    plotLFP(hObject, eventdata, handles)
    handles = zebraguidata(hObject);
    plotSUall(hObject, eventdata, handles)        
else axes(handles.axes10); cla;
    axes(handles.axes11); cla;
end
handles = zebraguidata(hObject); 
% Hint: get(hObject,'Value') returns toggle state of checkbox_plot


% --- Executes on button press in checkbox_plotLFPSUmax.
function checkbox_plotLFPSUmax_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plotLFPSUmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% plot the positive single units only
if get(hObject, 'Value')
    plotLFP(hObject, eventdata, handles)
    handles = zebraguidata(hObject);
    plotSUmax(hObject, eventdata, handles)        
else axes(handles.axes10); cla;
    axes(handles.axes11); cla;
end
handles = zebraguidata(hObject); 
% Hint: get(hObject,'Value') returns toggle state of checkbox_plotLFPSUmax


% --- Executes on button press in checkbox_plotLFPSUmin.
function checkbox_plotLFPSUmin_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plotLFPSUmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% plot the negative single units only
if get(hObject, 'Value')
    plotLFP(hObject, eventdata, handles)
    handles = zebraguidata(hObject);
    plotSUmin(hObject, eventdata, handles)    
else axes(handles.axes10); cla;
    axes(handles.axes11); cla;
end
handles = zebraguidata(hObject);
% Hint: get(hObject,'Value') returns toggle state of checkbox_plotLFPSUmin



function edit_tstart_plot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tstart_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tstart_plot as text
%        str2double(get(hObject,'String')) returns contents of edit_tstart_plot as a double


% --- Executes during object creation, after setting all properties.
function edit_tstart_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tstart_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_scroll_left.
function pushbutton_scroll_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_scroll_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%function to scroll though the plots of LFP and single units
if get(hObject, 'Value')
    currentx_min = str2double(get(handles.edit_tstart_plot, 'String'));
    currentx_max = str2double(get(handles.edit_tend_plot, 'String'));
    % find the current distance between the start and end time points
    diff = currentx_max-currentx_min;
    % use it to calculate the new timepoints
    newx_min = currentx_min-diff;
    newx_max = currentx_max-diff;
    if newx_min >= 0 && newx_max <= handles.time(end) %test that the new time points fall within the time range of the entire recording
        % if yes, you can undate the tstart and tend edit boxes
        set(handles.edit_tstart_plot, 'string', num2str(newx_min)); 
        set(handles.edit_tend_plot, 'string', num2str(newx_max));
        %and then re-execute the plotting
        handles = zebraguidata(hObject);
        plotLFP(hObject, eventdata, handles);
        handles = zebraguidata(hObject);
        % find out which single unit plotting function to use
        if get(handles.checkbox_plot, 'Value')
            plotSUall(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmin, 'Value')
            plotSUmin(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmax, 'Value')
            plotSUmax(hObject, eventdata, handles);
        end
    else axes(handles.axes10); cla; % if the time range is outside of the recording range, empty the plots
        axes(handles.axes11); cla;
    end
    % reset the pushbutton value to 0, so that when you click again it
    % changes to 1 and executes the code to scroll again, so you don't have
    % to unclick and click again
    set(hObject,'value',0);
end

% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%function to reset the time to the entire range of the LFP recording
if get(hObject, 'Value')
    % set the edit boxes of the start time and end time to empty again
    set(handles.edit_tstart_plot, 'string', '');
    set(handles.edit_tend_plot, 'string', '');
    % then replot everything, the plotting functions will automatically
    % determine the entire range
    plotLFP(hObject, eventdata, handles);
    handles = zebraguidata(hObject);
    % find out which single unit plotting function to use
    if get(handles.checkbox_plot, 'Value')
        plotSUall(hObject, eventdata, handles);
    elseif get(handles.checkbox_plotLFPSUmin, 'Value')
        plotSUmin(hObject, eventdata, handles);
    elseif get(handles.checkbox_plotLFPSUmax, 'Value')
        plotSUmax(hObject, eventdata, handles);
    end
 % and again reset the value of the pushbutton
 set(hObject,'value',0); 
end

% --- Executes on button press in pushbutton_scroll_right.
function pushbutton_scroll_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_scroll_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% function to scroll right through the plots, the same as the left
% scrolling function above
if get(hObject, 'Value')
    currentx_min = str2double(get(handles.edit_tstart_plot, 'String'));
    currentx_max = str2double(get(handles.edit_tend_plot, 'String'));
    diff = currentx_max-currentx_min;    
    newx_min = currentx_min+diff;
    newx_max = currentx_max+diff;
    if newx_min >= 0 && newx_max <= handles.time(end)
        set(handles.edit_tstart_plot, 'string', num2str(newx_min));
        set(handles.edit_tend_plot, 'string', num2str(newx_max));
        handles = zebraguidata(hObject);
        plotLFP(hObject, eventdata, handles);
        handles = zebraguidata(hObject);
        if get(handles.checkbox_plot, 'Value')
            plotSUall(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmin, 'Value')
            plotSUmin(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmax, 'Value')
            plotSUmax(hObject, eventdata, handles);
        end
    else axes(handles.axes10); cla;
        axes(handles.axes11); cla;
    end
    set(hObject,'value',0);
end


function edit_tend_plot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tend_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tend_plot as text
%        str2double(get(hObject,'String')) returns contents of edit_tend_plot as a double


% --- Executes during object creation, after setting all properties.
function edit_tend_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tend_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on edit_tstart_plot and none of its controls.
function edit_tstart_plot_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_tstart_plot (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% key press function: if enter is typed within the start time edit text
% box, the plotting functions are re-executed
k = get(gcf, 'CurrentKey');
if strcmp(k, 'return') % if the latest pressed key is enter, also called return  
        plotLFP(hObject, eventdata, handles);
        handles = zebraguidata(hObject);
        if get(handles.checkbox_plot, 'Value')
            plotSUall(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmin, 'Value')
            plotSUmin(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmax, 'Value')
            plotSUmax(hObject, eventdata, handles);
        end 
end


% --- Executes on key press with focus on edit_tend_plot and none of its controls.
function edit_tend_plot_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_tend_plot (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% same key press function for the end time
k = get(gcf, 'CurrentKey');
if strcmp(k, 'return');   
        plotLFP(hObject, eventdata, handles);
        handles = zebraguidata(hObject);
        if get(handles.checkbox_plot, 'Value')
            plotSUall(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmin, 'Value')
            plotSUmin(hObject, eventdata, handles);
        elseif get(handles.checkbox_plotLFPSUmax, 'Value')
            plotSUmax(hObject, eventdata, handles);
        end 
end



function edit_acqF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_acqF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_acqF as text
%        str2double(get(hObject,'String')) returns contents of edit_acqF as a double


% --- Executes during object creation, after setting all properties.
function edit_acqF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_acqF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_calculateSUUPDOWN.
function checkbox_calculateSUUPDOWN_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_calculateSUUPDOWN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value')
    Extract_SWA (hObject, eventdata, handles);
    handles = zebraguidata(hObject);
    Extract_SUpos_SWA (hObject, eventdata, handles);
    handles = zebraguidata(hObject);
    Extract_SUneg_SWA (hObject, eventdata, handles);
    handles = zebraguidata(hObject);
    % creation of a list containing all output data
    % first column: entire trace
    % second column: in up states
    % third column: in down states
    list = {1 1 1};
    % first row: number of positive peaks
    list(1 ,1) = num2cell(handles.nmax);
    list(1, 2) = num2cell(handles.pos_su_in_us);
    list(1, 3) = num2cell(handles.pos_su_in_ds);
    % second row: fraction of postive peaks
    list(2, 1) = num2cell((handles.nmax/handles.nmax)*100);
    list(2, 2) = num2cell((handles.pos_su_in_us/handles.nmax)*100);
    list(2, 3) = num2cell((handles.pos_su_in_ds/handles.nmax)*100);
    % third row: number of negative peaks
    list(3, 1) = num2cell(handles.nmin);
    list(3, 2) = num2cell(handles.neg_su_in_us);
    list(3, 3) = num2cell(handles.neg_su_in_ds);
    % fourth row: fraction of negative peaks
    list(4, 1) = num2cell((handles.nmin/handles.nmin)*100);
    list(4, 2) = num2cell((handles.neg_su_in_us/handles.nmin)*100);
    list(4, 3) = num2cell((handles.neg_su_in_ds/handles.nmin)*100);
    % fifth row: number of total peaks
    total = handles.nmin+handles.nmax;
    total_in_us = handles.neg_su_in_us+handles.pos_su_in_us;
    total_in_ds = handles.neg_su_in_ds+handles.pos_su_in_ds;
    list(5 ,1) = num2cell(total);
    list(5, 2) = num2cell(total_in_us);
    list(5, 3) = num2cell(total_in_ds);
    % sixth row: fraction of total peaks
    list(6, 1) = num2cell((total/total)*100);
    list(6, 2) = num2cell((total_in_us/total)*100);
    list(6, 3) = num2cell((total_in_ds/total)*100);
    
    % finally set this list as output data in the table
    set(handles.uitable_output,'Data',list);
else emptylist = {};
    set(handles.uitable_output, 'Data', emptylist);
end

% Hint: get(hObject,'Value') returns toggle state of checkbox_calculateSUUPDOWN


% --- Executes on button press in checkbox_savetimes.
function checkbox_savetimes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_savetimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value')
    acqF = str2double(get(handles.edit_acqF, 'string'));
    nus = handles.nus;
    nds = handles.nds;
    %create the several columns first: first timing positive single units
    %states, then timing negative single units, timing up states start and end and timing
    %down states start and end    
    A = handles.maxAt./acqF;
    B = handles.minAt./acqF;
    C = handles.zebraMainData.evDetData.upstates.fromT(1:nus);
    D = handles.zebraMainData.evDetData.upstates.toT(1:nus);
    E = handles.zebraMainData.evDetData.downstates.fromT(1:nds);    
    F = handles.zebraMainData.evDetData.downstates.toT(1:nds); 
    
    % then create the to be printed struct    
    output = struct('TimeposSU', A, 'TimenegSU', B, 'UPstart', C, 'UPend', D, 'DOWNstart', E, 'DOWNend', F);    
   
    
    %then create the text file. First find the correct name. I added the
    %possibility to add ltp number so that when several sums are recorded
    %at the ephys setup, which appear in the same folder, each of these
    %sums the data is separately stored    
    name = get(handles.edit_ltpnumber, 'string');
    channel = get(handles.chnumber, 'String');
    fileout = [handles.zebraMainData.dir_in 'SUanalysis' name 'Ch' channel '.dat'];     
    lengthstruct = length(fieldnames(output));
    fields = fieldnames(output);
    % find the field with the longest length
    ls = [handles.nmax handles.nmin handles.nus handles.nds];
    maximum = max(ls);
    
    % open and start printing
    fid = fopen(fileout, 'w');
    fprintf(fid,'%s\t', fields {:});
    fprintf(fid, '\n');
    for i =  1:maximum 
        for ff = 1:lengthstruct 
            if i > length(output.(fields{ff}))
                fprintf(fid, '\t');
            else fprintf(fid,'%g\t', output.(fields{ff})(i));            
            end 
        end
        fprintf(fid,'\n');
    end
    fid = fclose(fid);   

end    
    
% Hint: get(hObject,'Value') returns toggle state of checkbox_savetimes



function edit_ltpnumber_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ltpnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ltpnumber as text
%        str2double(get(hObject,'String')) returns contents of edit_ltpnumber as a double


% --- Executes during object creation, after setting all properties.
function edit_ltpnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ltpnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
