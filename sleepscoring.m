function varargout = sleepscoring(varargin)
% SLEEPSCORING MATLAB code for sleepscoring.fig
%      SLEEPSCORING, by itself, creates a new SLEEPSCORING or raises the existing
%      singleton*.
%
%      H = SLEEPSCORING returns the handle to a new SLEEPSCORING or the handle to
%      the existing singleton*.
%
%      SLEEPSCORING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLEEPSCORING.M with the given input arguments.
%
%      SLEEPSCORING('Property','Value',...) creates a new SLEEPSCORING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sleepscoring_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sleepscoring_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sleepscoring

% Last Modified by GUIDE v2.5 02-Nov-2018 17:26:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sleepscoring_OpeningFcn, ...
                   'gui_OutputFcn',  @sleepscoring_OutputFcn, ...
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


% --- Executes just before sleepscoring is made visible.
function sleepscoring_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sleepscoring (see VARARGIN)

% Choose default command line output for sleepscoring
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
%handles.numbersleep = 0;
%handles.timesleep = [];
%handles.numberhightheta = 0;
%handles.timehightheta = [];
% Update handles structure
zebraguidata(hObject, handles);
% UIWAIT makes sleepscoring wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sleepscoring_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function find_NREM_episodes(hObject, eventdata, handles)
% let's start by creating a vector contaning time points in the left column
% and theta over delta band power in the right column
handles = zebraguidata(hObject);
channel = handles.zebraMainData.currentCh;
timepoints_temp = handles.zebraMainData.BPpower(1,2,channel,1,:);
timepoints = timepoints_temp(:,:)';
theta_temp = handles.zebraMainData.BPpower(2,3,channel,1,:);
theta = theta_temp(:,:)';
delta_temp = handles.zebraMainData.BPpower(2,2,channel,1,:);
delta = delta_temp(:,:)';
handles.ratio = delta./theta;
len = length(timepoints);
handles.numbersleep = 0; % settin the to be calculated values/matrixes
handles.timesleep = []; % settin the to be calculated values/matrixes
handles.endtimesleep = []; % settin the to be calculated values/matrixes
% we can now look at episodes below threshold of 0.75 that count as sleep
for i = 1:len % loop through the ratio matrix
    if i == 1 && handles.ratio(i) >= 1.25 % for the first time point, if it is lower than 0.75 it is the first sleep episodes
        handles.numbersleep = handles.numbersleep + 1; %increase the sleep episodes number by one
        handles.timesleep(handles.numbersleep) = timepoints(i); % use this time as the begin time of the current sleep episodes
    elseif i > 1 && handles.ratio(i) >= 1.25 && handles.ratio(i-1) < 1.25 % we test that it is the first time point lower than threshold
        handles.numbersleep = handles.numbersleep + 1; % then we count it as a new sleep episodes
        handles.timesleep(handles.numbersleep) = timepoints(i); % and label the start time for the current episode
    elseif i > 1 && handles.ratio(i) < 1.25 && handles.ratio(i-1) >= 1.25 % if it reaches over threshhold again
        handles.endtimesleep(handles.numbersleep) = timepoints(i); % we label it as the endpoint of the current episode
    elseif i == len && handles.ratio(i) >= 1.25 % if we reached the end of the trace and we are still in an episode
        handles.endtimesleep(handles.numbersleep) = timepoints(i); % we label this time point as the end of the episode
    end
end
% we can now calculate length of episodes and the end of each episode
if length(handles.endtimesleep) < length(handles.timesleep)
    handles.endtimesleep(handles.numbersleep) = timepoints(end);
end
handles.sleeplength = handles.endtimesleep-handles.timesleep;
handles.averagesleeplength = mean(handles.sleeplength);

% we now do the same for high theta episodes (higher than 1.25)
handles.numberhightheta = 0;
handles.timehightheta = [];
handles.endtimehightheta = [];
for ii = 1:len
    if ii == 1 && handles.ratio(ii)<= 0.75
        handles.numberhightheta = handles.numberhightheta + 1;
        handles.timehightheta(handles.numberhightheta) = timepoints(ii);
    elseif ii > 1 && handles.ratio(ii) <= 0.75 && handles.ratio(ii-1)> 0.75
        handles.numberhightheta = handles.numberhightheta + 1;
        handles.timehightheta(handles.numberhightheta) = timepoints(ii);
    elseif ii > 1 && handles.ratio(ii) > 0.75 && handles.ratio(ii-1) <= 0.75
        handles.endtimehightheta(handles.numberhightheta) = timepoints(ii);
    elseif ii == len && handles.ratio(ii) <= 0.75
        handles.endtimehightheta(handles.numberhightheta) = timepoints(ii);
    end
end
if length(handles.endtimehightheta) < length(handles.timehightheta)
    handles.endtimehightheta(handles.numberhightheta) = timepoints(end);
end
handles.highthetalength = handles.endtimehightheta-handles.timehightheta;
handles.averagehighthetalength = mean(handles.highthetalength);


zebraguidata(hObject, handles);

function findpower(hObject, eventdata, handles)
% Will find the average delta and theta within sleep and high theta episodes
handles = zebraguidata(hObject);
% Then I have to obtain the delta and theta powers again

channel = handles.zebraMainData.currentCh;
timepoints_temp = handles.zebraMainData.BPpower(1,2,channel,1,:);
timepoints = timepoints_temp(:,:)';
stepsize = timepoints(2)-timepoints(1); % needed to find the number of time steps within a sleep episode
theta_temp = handles.zebraMainData.BPpower(2,3,channel,1,:);
theta = theta_temp(:,:)';
delta_temp = handles.zebraMainData.BPpower(2,2,channel,1,:);
delta = delta_temp(:,:)';
len_sleep_vector = handles.numbersleep;
len_theta_vector = handles.numberhightheta;
deltainsleep = zeros(len_sleep_vector, 1); % set empty vector
thetainsleep = zeros(len_sleep_vector, 1); % set empty vector
for i = 1:len_sleep_vector % for each sleep episode
    d = 0; % we set the total deltapower to 0
    t = 0; % same for theta power
    l = handles.sleeplength(i); % calculate the time length of the episode
    if l == 0
        break
    end
    steps = l/stepsize; % in order to find the number of time steps within the episode
    for tt = handles.timesleep(i):stepsize:handles.endtimesleep(i) % then loop through from begin to endtime of the episode        
        ttt = find(abs(timepoints-tt)<0.001); % find the index of the current timepoint        
        d = delta(ttt) + d; % and add the delta of each time point to the total delta
        t = theta(ttt) + t; % same for theta        
    end
    deltainsleep(i) = d/steps; % we then find average delta by dividing total delta of episode by number of steps in episode 
    thetainsleep(i) = t/steps; % same for theta
end  
% now we can find the average delta power of all the episodes
handles.averagedeltasleep = mean(deltainsleep);
handles.averagethetasleep = mean(thetainsleep);
% so we do the same for the high theta episodes
deltaintheta = zeros(len_theta_vector,1);
thetaintheta = zeros(len_theta_vector, 1);
for ii = 1:len_theta_vector
    d2 = 0;
    t2 = 0;
    l2 = handles.highthetalength(ii);
    if l2 == 0
        break
    end
    steps2 = l2/stepsize;
    for c = handles.timehightheta(ii):stepsize:handles.endtimehightheta(ii)
        cc = find(abs(timepoints-c)<0.001);
        d2 = d2 + delta(cc);
        t2 = t2 + theta(cc);        
    end
    deltaintheta(ii) = d2/steps2;
    thetaintheta(ii) = t2/steps2;
end
handles.averagedeltaintheta = mean(deltaintheta);
handles.averagethetaintheta = mean(thetaintheta);
% end by updating handles structure
zebraguidata(hObject, handles);

% --- Executes on button press in checkbox_calculate.
function checkbox_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value')    
    find_NREM_episodes(hObject,handles);
    handles = zebraguidata(hObject);
    set(handles.edit_nsleep, 'string', num2str(handles.numbersleep));
    set(handles.edit_nhightheta, 'string', num2str(handles.numberhightheta));
    set(handles.edit_avlengthsleep, 'string', num2str(handles.averagesleeplength));
    set(handles.edit_avlengthhightheta, 'string', num2str(handles.averagehighthetalength));
    
    findpower(hObject, handles);
    handles = zebraguidata(hObject);
    set(handles.edit_avdeltapowersleep, 'string', num2str(handles.averagedeltasleep));
    set(handles.edit_avthetapowersleep, 'string', num2str(handles.averagethetasleep));
    set(handles.edit_avdeltapowerhightheta, 'string', num2str(handles.averagedeltaintheta));
    set(handles.edit_avthetapowerhightheta, 'string', num2str(handles.averagethetaintheta));
    
    % plot a histogram ratio distribution
    axes(handles.axes_histo);
    mratio = (ceil(max(handles.ratio)/0.04))*0.04;
    edges_ratio = 0.1:0.04:mratio;
    handles.h_ratio = histogram(handles.ratio, edges_ratio);
    handles.h_ratio.FaceColor = 'g';
    xlabel('delta/beta-sigma ratio'), ylabel('frequency')
    
    % plot a histogram length sleep distribution
    axes(handles.axes_histosleepduration) 
    m = (ceil(max(handles.sleeplength)/0.8))*0.8;
    edges = 0:0.8:m;
    h_sleepduration = histogram(handles.sleeplength, edges);
    h_sleepduration.FaceColor = 'r';
    xlabel('sleep length (s)'), ylabel('frequency');
    
    %plot a cumulative distribution of sleep length
    axes(handles.axes_cdf);
    %nr = normalize(handles.sleeplength, 'range');    
    handles.c = cdfplot(handles.sleeplength);    
    xlabel('sleep length (s)'), ylabel('cumulative distribution');  

    % update handles structure
    zebraguidata(hObject, handles);
    
else
    set(handles.edit_nsleep, 'string', ' ');
    set(handles.edit_nhightheta, 'string', ' ');
    set(handles.edit_avlengthsleep, 'string', ' ');
    set(handles.edit_avlengthhightheta, 'string', ' ');
    set(handles.edit_avdeltapowersleep, 'string', ' ');
    set(handles.edit_avthetapowersleep, 'string', ' ');
    set(handles.edit_avdeltapowerhightheta, 'string', ' ');
    set(handles.edit_avthetapowerhightheta, 'string', ' '); 
    axes(handles.axes_histo); 
    cla;
    axes(handles.axes_histosleepduration);
    cla;
    axes(handles.axes_cdf);
    cla;
end


function edit_nsleep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nsleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nsleep as text
%        str2double(get(hObject,'String')) returns contents of edit_nsleep as a double


% --- Executes during object creation, after setting all properties.
function edit_nsleep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nsleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nhightheta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nhightheta as text
%        str2double(get(hObject,'String')) returns contents of edit_nhightheta as a double

% --- Executes during object creation, after setting all properties.
function edit_nhightheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_avlengthsleep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avlengthsleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avlengthsleep as text
%        str2double(get(hObject,'String')) returns contents of edit_avlengthsleep as a double


% --- Executes during object creation, after setting all properties.
function edit_avlengthsleep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avlengthsleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_avlengthhightheta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avlengthhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avlengthhightheta as text
%        str2double(get(hObject,'String')) returns contents of edit_avlengthhightheta as a double


% --- Executes during object creation, after setting all properties.
function edit_avlengthhightheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avlengthhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_avdeltapowersleep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avdeltapowersleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avdeltapowersleep as text
%        str2double(get(hObject,'String')) returns contents of edit_avdeltapowersleep as a double


% --- Executes during object creation, after setting all properties.
function edit_avdeltapowersleep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avdeltapowersleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_avdeltapowerhightheta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avdeltapowerhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avdeltapowerhightheta as text
%        str2double(get(hObject,'String')) returns contents of edit_avdeltapowerhightheta as a double


% --- Executes during object creation, after setting all properties.
function edit_avdeltapowerhightheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avdeltapowerhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_avthetapowersleep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avthetapowersleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avthetapowersleep as text
%        str2double(get(hObject,'String')) returns contents of edit_avthetapowersleep as a double


% --- Executes during object creation, after setting all properties.
function edit_avthetapowersleep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avthetapowersleep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_avthetapowerhightheta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_avthetapowerhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_avthetapowerhightheta as text
%        str2double(get(hObject,'String')) returns contents of edit_avthetapowerhightheta as a double


% --- Executes during object creation, after setting all properties.
function edit_avthetapowerhightheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_avthetapowerhightheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_saveplots.
function checkbox_saveplots_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_saveplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value')    
    %save the y data of the cdf
    fileOutY = [handles.zebraMainData.dir_in 'cdf-nnm-y.dat'];        
    fid = fopen(fileOutY,'w');
    for ii = 1:size(handles.c.YData)
        fprintf(fid,'%g\t',handles.c.YData(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    %save the x data of the cdf
    fileOutX = [handles.zebraMainData.dir_in 'cdf-nnm-x.dat'];        
    fid = fopen(fileOutX,'w');
    for ii = 1:size(handles.c.XData)
        fprintf(fid,'%g\t',handles.c.XData(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    % save the ratio histogram Y data
    fileOutRY = [handles.zebraMainData.dir_in 'ratio-histo-y.dat'];        
    fid = fopen(fileOutRY,'w');
    for ii = 1:size(handles.h_ratio.BinEdges);
        fprintf(fid,'%g\t',handles.h_ratio.BinEdges(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    % save the ratio histogram X Data (should be the same for all plots
    % though)
    fileOutRX = [handles.zebraMainData.dir_in 'ratio-histo-x.dat'];        
    fid = fopen(fileOutRX,'w');
    for ii = 1:size(handles.h_ratio.Values);
        fprintf(fid,'%g\t',handles.h_ratio.Values(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
       
end
    
% Hint: get(hObject,'Value') returns toggle state of checkbox_saveplots
