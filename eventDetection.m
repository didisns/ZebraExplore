function varargout = eventDetection(varargin)
% EVENTDETECTION MATLAB code for eventDetection.fig
%      EVENTDETECTION, by itself, creates a new EVENTDETECTION or raises the existing
%      singleton*.
%
%      H = EVENTDETECTION returns the handle to a new EVENTDETECTION or the handle to
%      the existing singleton*.
%
%      EVENTDETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENTDETECTION.M with the given input arguments.
%
%      EVENTDETECTION('Property','Value',...) creates a new EVENTDETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eventDetection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eventDetection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eventDetection

% Last Modified by GUIDE v2.5 09-Apr-2018 18:20:16

% Modificaton November 2018 by Didi: addition of display medianV, SD, and
% pw of the other hemisphere in order to be able to find amplitude
% of up states in the mosaic hemisphere that are hardly present.
% NOTE: this code is written assuming there are two hemispheres, i.e. two
% channels, if there are more, nothing will happen to the program and it
% will calculate only from the selected channel as usual

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eventDetection_OpeningFcn, ...
                   'gui_OutputFcn',  @eventDetection_OutputFcn, ...
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
end

function eventDetection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eventDetection (see VARARGIN)

% Choose default command line output for eventDetection
handles.output = hObject;

% establish the link with the calling process
handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

% local structures and variables
%handles.USlist = struct('N',1,'tfrom',1,'tto',1,'delta',1);
%handles.USlist = {1 1 1 1 1 1};

%handles.upstates = struct('fromI',1,'fromT',1,'toI',1,'toT',1,'deltaT',1,'medianV',1,'SD',1,'pw',1,...
%    'medianV_hem2', 1, 'SD_hem2', 1, 'pw_hem2', 1, 'enable',1);
%handles.downstates = struct('fromI',1,'fromT',1,'toI',1,'toT',1,'deltaT',1,'medianV',1,'SD',1,'pw',1,...
%    'medianV_hem2', 1, 'SD_hem2', 1, 'pw_hem2', 1, 'enable',1);

% Update handles structure
zebraguidata(hObject, handles);
end

function varargout = eventDetection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function baselineThr_Callback(hObject, eventdata, handles)
end

function baselineThr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function minDuration_Callback(hObject, eventdata, handles)
end

function minDuration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function maxInterruption_Callback(hObject, eventdata, handles)
end

function maxInterruption_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function eventThr_Callback(hObject, eventdata, handles)
end

function eventThr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function extractEvents_Callback(hObject, eventdata, handles)
% refers to the 'compute' button
handles.parentHandles = zebraguidata(handles.parentObject);    % this to refresh the local values of parentHandles

handles.upstates = struct('fromI',1,'fromT',1,'toI',1,'toT',1,'deltaT',1,'medianV',1,'SD',1,'pw',1,...
    'medianV_hem2', 1, 'SD_hem2', 1, 'pw_hem2', 1, 'enable',1);
handles.downstates = struct('fromI',1,'fromT',1,'toI',1,'toT',1,'deltaT',1,'medianV',1,'SD',1,'pw',1,...
    'medianV_hem2', 1, 'SD_hem2', 1, 'pw_hem2', 1, 'enable',1);

handles.SPG = [];
handles.LFP = [];
handles.LFP_hem2 = [];
handles.timeSPG = [];
handles.upstates.fromI = [];
handles.upstates.toI = [];
handles.downstates.fromI = [];
handles.downstates.toI = [];
handles.USlist = [];
handles.USlist = {1 1 1 1 1 1 1 1 1 1};

% set the threshold for event detection based on the used method
if get(handles.autoThr,'Value') % if autosetting state thresshold is clicked, automatically determine threshold
    lowThr = handles.autoblThr; % this value is calculated later in analyzeData function based on power distribution!
    highThr = handles.autoevThr; % this value is calculated later in analyzeData function based on power distribution!
else
    lowThr = str2double(get(handles.baselineThr,'String')); %otherwise, obtain user defined threshold values
    highThr = str2double(get(handles.eventThr,'String'));
end

% Identify the data set to be analyzed
dataSource = 9;
if get(handles.detectionBP1,'Value'), 
    if handles.parentHandles.BPcomputed(2), dataSource = 2; %BPcomputed is a vector containing a 0 if the 
       % specific band pass has not been analyzed, and a 1 if it has been analyzed
    end
end
if get(handles.detectionBP2,'Value'),
    if handles.parentHandles.BPcomputed(3), dataSource = 3;
    end
end
if get(handles.detectionBP3,'Value'), 
    if handles.parentHandles.BPcomputed(4), dataSource = 4;
    end
end
if get(handles.detectionHP,'Value'), 
    if handles.parentHandles.BPcomputed(5), dataSource = 5;
    end
end
if get(handles.detectionSPG,'Value'), dataSource = 0;
end

if dataSource < 9
    ct = handles.parentHandles.currentTrial;
    cCh = handles.parentHandles.currentCh;
    len = handles.parentHandles.dtaLen;
    spp = handles.parentHandles.sp;
    logFlag = get(handles.powerLogCk,'Value'); % if 'log of envelope power' is checked
    flagLeak = get (handles.parentHandles.displayLeakCk,'Value');
    handles.tmin = str2double(get(handles.leftTime,'String')); % time values entered by user in the bottom left box
    handles.tmax = str2double(get(handles.rightTime,'String'));

    handles.LFP(1:len) = handles.parentHandles.workLFP (cCh,ct,1:len);
    if handles.parentHandles.nCh == 2
        if cCh == 1
            Ch_hem2 = 2;
        elseif cCh == 2 
            Ch_hem2 = 1;
        end
        handles.LFP_hem2(1:len) = handles.parentHandles.workLFP (Ch_hem2, ct, 1:len);
    end
        

    if dataSource == 0          % length of the spectral power file
        sl = handles.parentHandles.spgl;
        handles.timeSPG = handles.parentHandles.spgt;
    else
        sl = handles.parentHandles.BPpowerL(dataSource); % length of the vector power of band passed data 
        handles.timeSPG(1:sl) = handles.parentHandles.BPpower(1,dataSource,cCh,ct,1:sl); % spectrogram data itself for current BP setting
    end

    if (handles.parentHandles.deleakedFlag && flagLeak) % If we want the deleaked version of the spg
        if logFlag, handles.SPG (1:sl) = log10(handles.parentHandles.spgPlotDeleaked(cCh,ct,1:sl)); %calculate the log of spg if clicked
        else
            handles.SPG (1:sl) = handles.parentHandles.spgPlotDeleaked(cCh,ct,1:sl);
        end
    else    
        if logFlag, % if not the deleaked version, but log is clicked:
            switch dataSource
                case 0      % SPG
                    handles.SPG (1:sl) = log10(handles.parentHandles.spgPlot(cCh,ct,1:sl));
                case {2, 3, 4, 5}      % BP1, BP2, BP3, HP 
                    handles.SPG (1:sl) = log10(handles.parentHandles.BPpower(2,dataSource,cCh,ct,1:sl));
            end
        else % finally, if not the the deleaked version, and log is not clicked:
            handles.SPG (1:sl) = handles.parentHandles.spgPlot(cCh,ct,1:sl);
            switch dataSource
                case 0      % SPG
                    handles.SPG (1:sl) = handles.parentHandles.spgPlot(cCh,ct,1:sl);
                case {2, 3, 4, 5}      % BP1, BP2, BP3, HP 
                    handles.SPG (1:sl) = handles.parentHandles.BPpower(2,dataSource,cCh,ct,1:sl);
            end
        end
    end

    downStatesBolean = handles.SPG < lowThr; %logical array of SPG below threshold
    downStatesI = find(handles.SPG< lowThr); % indexes of SPG where it is below threshold (i.e. is downstate)
    upStatesBolean = handles.SPG > highThr; % same for up states (above threshold)
    upStatesI = find(handles.SPG>highThr);
    limboBolean = ~(downStatesBolean | upStatesBolean); % I don't understand this, it seems to create an empty logical array 

    % now the vector handles.spgt (downstatesBolean) contains the timing of the
    % upstates. This has to be converted to the indexes of the workLFP file

    lus = length (upStatesI);
    lds = length (downStatesI);

    % process upstates first
    if lus>0 
        firstI = upStatesI (1);
        firstTime = handles.timeSPG (firstI); 
        handles.upstates.fromT (1) = firstTime;
        handles.upstates.fromI (1) = 1+firstTime/spp; % spp is the sampling period, this converts the time to the LFP index
        nus = 1; % keeps track of the number of up states
        for i=2:lus
            nextI = upStatesI(i); % loop through all the indexes with higher than threshold value
            if nextI==firstI+1
                % continue the old US
                firstI=nextI;
                if i == lus
                    % we have reached the end of the track. close the US
                    endTime = handles.timeSPG (firstI);
                    handles.upstates.toT (nus) = endTime;
                    handles.upstates.toI (nus) = int32(1+endTime/spp);
                    handles.upstates.deltaT (nus) = endTime-handles.upstates.fromT (nus);
                    handles.upstates.enable (nus) = 1;
                end
            else
                firstTime = handles.timeSPG (firstI);
                % beginning of a new US. Close the old US
                handles.upstates.toT (nus) = firstTime;
                handles.upstates.toI (nus) = int32(1+firstTime/spp);
                handles.upstates.deltaT (nus) = firstTime-handles.upstates.fromT (nus);
                handles.upstates.enable (nus) = 1;
                nus = nus + 1;
                % open new US
                handles.upstates.fromT (nus) = handles.timeSPG (nextI);
                handles.upstates.fromI (nus) =  int32(1+handles.timeSPG (nextI)/spp);
                firstI = nextI;
            end    
        end
    end

    nus = length (handles.upstates.toI);      % just to make sure...

    % process downstates
    if nus>0
        firstI = downStatesI (1);
        firstTime = handles.timeSPG (firstI); 
        handles.downstates.fromT (1) = firstTime;
        handles.downstates.fromI (1) = 1+firstTime/spp;
        nds = 1;

        for i=2:lds
            nextI = downStatesI(i);
            if nextI==firstI+1
                % continue the old DS
                firstI=nextI;
                 if i == lds
                      % we have reached the end of the track. close the DS
                      endTime = handles.timeSPG (firstI);
                      handles.downstates.toT (nds) = endTime;
                      handles.downstates.toI (nds) = int32(1+endTime/spp);
                      handles.downstates.deltaT (nds) = endTime-handles.downstates.fromT (nds);
                      handles.downstates.enable (nds) = 1;
                 end
            else
                firstTime = handles.timeSPG (firstI);
                % beginning of a new DS. Close the old DS
                handles.downstates.toT (nds) = firstTime;
                handles.downstates.toI (nds) = int32(1+firstTime/spp);
                handles.downstates.deltaT (nds) = firstTime-handles.downstates.fromT (nds);
                handles.downstates.enable (nds) = 1;
                nds = nds + 1;
                % open new DS
                handles.downstates.fromT (nds) = handles.timeSPG (nextI);
                handles.downstates.fromI (nds) =  int32(1+handles.timeSPG (nextI)/spp);
                firstI = nextI;
            end
        end
    end    
    nds = length(handles.downstates.toI);

    zebraguidata(hObject, handles);

    % defragment DS and US across short state gap. The max size of the filled
    % gap is defined in the GUI
    % These are the rules: a brief interruption shorter than maxInterruption is
    % filled in by assigning to the neighboring US that are fused in a longer
    % one. Isolated US briefer than minDuration are attributed to the limbo.

    maxGap = str2double(get(handles.maxInterruption,'String'));
    minUS = str2double(get(handles.minDuration,'String'));

    if get(handles.removeCk,'Value') %if the remove state interruptions button is clicked
        % process US first
        usi = 2 
        while usi <= nus % there must be at least two up states to do this operation, loop through from 2 till nus
            % perform the fusions first since brief US might get fused together
            % to form a longer and legal US.
            interval = [handles.upstates.fromT(usi) handles.upstates.toT(usi-1)];
            distance = interval(1) - interval(2);
            if (distance<=maxGap)
                % first fuse the two upstates
                handles.upstates.toT(usi-1) = handles.upstates.toT(usi);
                handles.upstates.toI(usi-1) = handles.upstates.toI(usi);
                handles.upstates.deltaT(usi-1) = handles.upstates.toT(usi-1) - handles.upstates.fromT(usi-1);
                % second remove the US pointed to by usi
                for k = usi:nus-1
                    % shift all USs
                    handles.upstates.toT(k) = handles.upstates.toT(k+1);
                    handles.upstates.toI(k) = handles.upstates.toI(k+1); 
                    handles.upstates.fromT(k) = handles.upstates.fromT(k+1);
                    handles.upstates.fromI(k) = handles.upstates.fromI(k+1); 
                    handles.upstates.deltaT(k) = handles.upstates.deltaT (k+1);
                end
                nus = nus - 1;
            else
                usi = usi+1;
            end    
        end
        % remove brief states 
        for usi=1:nus
            if handles.upstates.deltaT(usi) <= minUS, handles.upstates.enable(usi) = 0; 
            end 
        end
        % remove the disabled states from the list
        cnt = nus;
        for usi=nus:-1:1
            if handles.upstates.enable(usi)
                % do nothing
            else
                cnt = cnt - 1;
                % shift down if the US is not the last one of the track
                if usi < nus
                    handles.upstates.fromI (cnt:-1:usi) = handles.upstates.fromI (cnt+1:-1:usi+1);
                    handles.upstates.toI (cnt:-1:usi) = handles.upstates.toI (cnt+1:-1:usi+1);
                    handles.upstates.fromT (cnt:-1:usi) = handles.upstates.fromT (cnt+1:-1:usi+1);
                    handles.upstates.toT (cnt:-1:usi) = handles.upstates.toT (cnt+1:-1:usi+1);
                    handles.upstates.deltaT (cnt:-1:usi) = handles.upstates.deltaT (cnt+1:-1:usi+1);                
                    handles.upstates.enable (cnt:-1:usi) = handles.upstates.enable (cnt+1:-1:usi+1);
                end
            end 
        end
        nus = cnt;

        % process DS
        for usi=2:nds
            interval = [handles.downstates.fromT(usi) handles.downstates.toT(usi-1)];
            distance = interval(1) - interval(2);
            if (distance<=maxGap)
                % fuse the contiguous downstates
                handles.downstates.toT(usi-1) = handles.downstates.toT(usi);
                handles.downstates.toI(usi-1) = handles.downstates.toI(usi);
                handles.downstates.deltaT(usi-1) = handles.downstates.toT(usi-1) - handles.downstates.fromT(usi-1);
                handles.downstates.enable(usi) = 0;    
            end    
        end

        for usi=1:nds
            if handles.downstates.deltaT(usi) <= minUS, handles.downstates.enable(usi) = 0; 
            end
        end

        cnt = nds;
        for usi=nds:-1:1
            if handles.downstates.enable(usi)
                % do nothing
            else
                cnt = cnt - 1;
                % shift down
                if usi < nds
                    handles.downstates.fromI (cnt:-1:usi) = handles.downstates.fromI (cnt+1:-1:usi+1);
                    handles.downstates.toI (cnt:-1:usi) = handles.downstates.toI (cnt+1:-1:usi+1);
                    handles.downstates.fromT (cnt:-1:usi) = handles.downstates.fromT (cnt+1:-1:usi+1);
                    handles.downstates.toT (cnt:-1:usi) = handles.downstates.toT (cnt+1:-1:usi+1);
                    handles.downstates.deltaT (cnt:-1:usi) = handles.downstates.deltaT (cnt+1:-1:usi+1);                
                    handles.downstates.enable (cnt:-1:usi) = handles.downstates.enable (cnt+1:-1:usi+1);                
                end
            end
        end
        nds = cnt;
    end
    
    
    % refine the edges of the state detection. Experimental and not great
    % not to be used as of March 22, 2017
    zebraguidata(hObject, handles);
    refineEdges (hObject, handles);
    handles = zebraguidata(hObject); 

    % At this stage all US and DS have been extracted and we can compute the
    % relative metrics.
    % Compute median, PW and SD of each US and DS.
    for usi=1:nds   % compute DS first since you need this to compute US size
       i1 = handles.downstates.fromI(usi);
       i2 = handles.downstates.toI(usi);
       handles.downstates.medianV(usi) = median(handles.LFP(i1:i2));       
       handles.downstates.SD(usi) = std(handles.LFP(i1:i2)); 
       handles.downstates.pw(usi) = rms(handles.LFP(i1:i2));
       if handles.parentHandles.nCh == 2
           handles.downstates.medianV_hem2(usi) = median(handles.LFP_hem2(i1:i2));
           handles.downstates.SD_hem2(usi) = std(handles.LFP_hem2(i1:i2)); 
           handles.downstates.pw_hem2(usi) = rms(handles.LFP_hem2(i1:i2));
       end
    end

    DSsearch = 1;
    dsi1 = handles.downstates.fromI(DSsearch);
    dsi2 = handles.downstates.toI(DSsearch);
    for usi=1:nus
       i1 = handles.upstates.fromI(usi);
       i2 = handles.upstates.toI(usi);
       handles.upstates.medianV(usi) = median(handles.LFP(i1:i2));
       handles.upstates.SD(usi) = std(handles.LFP(i1:i2));
       handles.upstates.pw(usi) = rms(handles.LFP(i1:i2)); 
       if handles.parentHandles.nCh == 2
           handles.upstates.medianV_hem2(usi) = median(handles.LFP_hem2(i1:i2));
           handles.upstates.SD_hem2(usi) = std(handles.LFP_hem2(i1:i2));
           handles.upstates.pw_hem2(usi) = rms(handles.LFP_hem2(i1:i2));
       end
       % now search for the closest DS
       while dsi2 < i1 && DSsearch <= nds % while the end of the downstate is smaller than the start of the current up
           % state and smaller than the number of down states
           DSsearch = DSsearch + 1; % it will increase by one: we found the down state following the current up state
           if DSsearch <= nds
               dsi2 = handles.downstates.toI(DSsearch);
           end
       end
       % now DSsearch points to the DS immediatey after unless the data ends
       % with a US. A second exception is when the data starts with an US
       if DSsearch == 1     % data begins with US
           baseline = handles.downstates.medianV(1);
           if handles.parentHandles.nCh == 2
               baseline_hem2 = handles.downstates.medianV_hem2(1);
           end
       else
           if DSsearch > nds    % data ends up with a US
               baseline = handles.downstates.medianV(nds);
               if handles.parentHandles.nCh == 2
                   baseline_hem2 = handles.downstates.medianV_hem2(nds);
               end
           else
               % OK, this is a middle of the road US!
               baseline = (handles.downstates.medianV(DSsearch)+handles.downstates.medianV(DSsearch-1))/2;
               if handles.parentHandles.nCh == 2
                   baseline_hem2 = (handles.downstates.medianV_hem2(DSsearch)+handles.downstates.medianV_hem2(DSsearch-1))/2;
               end
           end
       end 
       % this is added by Didi, in order to access the value later, so that
       % I can add a downward going filter
       newmedianV(usi) = handles.upstates.medianV(usi) - baseline;
       if handles.parentHandles.nCh == 2
            newmedianV_hem2(usi) = handles.upstates.medianV_hem2(usi) - baseline_hem2;
       end
    end
    
    % This part is added by Didi in November, 2018. Meant to remove
    % upstates from the list that have a positive medianV (i.e. that are
    % just artefacts)
    for us = 1:nus    
       i1 = handles.upstates.fromI(us);
       i2 = handles.upstates.toI(us);       
        if newmedianV(us) > 0
            handles.upstates.enable(us) = 0;
        end
    end 
    
    % now remove the disabled up states from the list as was done above
    cnt = nus;
    for usi=nus:-1:1
        if handles.upstates.enable(usi)
        % do nothing
        else
            cnt = cnt - 1;
            % shift down if the US is not the last one of the track
            if usi < nus
                handles.upstates.fromI (cnt:-1:usi) = handles.upstates.fromI (cnt+1:-1:usi+1);
                handles.upstates.toI (cnt:-1:usi) = handles.upstates.toI (cnt+1:-1:usi+1);
                handles.upstates.fromT (cnt:-1:usi) = handles.upstates.fromT (cnt+1:-1:usi+1);
                handles.upstates.toT (cnt:-1:usi) = handles.upstates.toT (cnt+1:-1:usi+1);
                handles.upstates.deltaT (cnt:-1:usi) = handles.upstates.deltaT (cnt+1:-1:usi+1);                
                handles.upstates.enable (cnt:-1:usi) = handles.upstates.enable (cnt+1:-1:usi+1);
                handles.upstates.medianV(cnt:-1:usi) = handles.upstates.medianV(cnt+1:-1:usi+1);
                newmedianV(cnt:-1:usi) = newmedianV(cnt+1:-1:usi+1);
                handles.upstates.SD(cnt:-1:usi) = handles.upstates.SD(cnt+1:-1:usi+1);
                handles.upstates.pw(cnt:-1:usi) = handles.upstates.pw(cnt+1:-1:usi+1);       
                if handles.parentHandles.nCh == 2
                    handles.upstates.medianV_hem2(cnt:-1:usi) = handles.upstates.medianV_hem2(cnt+1:-1:usi+1);
                    newmedianV_hem2(cnt:-1:usi) = newmedianV_hem2(cnt+1:-1:usi+1);
                    handles.upstates.SD_hem2(cnt:-1:usi) = handles.upstates.SD_hem2(cnt+1:-1:usi+1);
                    handles.upstates.pw_hem2(cnt:-1:usi) = handles.upstates.pw_hem2(cnt+1:-1:usi+1);
                end
            end
        end 
        end
        nus = cnt;

    handles.NUS = nus;
    handles.NDS = nds;    

    % now create the list:
    
    for usi = 1:nus
       handles.USlist (usi,5) = num2cell(newmedianV(usi));   
       handles.USlist (usi,6) = num2cell(handles.upstates.SD(usi));   
       handles.USlist (usi,7) = num2cell(handles.upstates.pw(usi));
       if handles.parentHandles.nCh == 2
           handles.USlist (usi,8) = num2cell(newmedianV_hem2(usi));   
           handles.USlist (usi,9) = num2cell(handles.upstates.SD_hem2(usi));   
           handles.USlist (usi,10) = num2cell(handles.upstates.pw_hem2(usi));
       end
    end
    
    % Creation of the list of USs 
    %handles.USlist (1,1) = num2cell(1);
    %handles.USlist (1,2) = num2cell(handles.upstates.fromT(1));
    handles.USlist (1:nus,1) = num2cell(1:nus); % the first cell contains the up state number
    handles.USlist (1:nus,2) = num2cell(handles.upstates.fromT(1:nus)); % the second cell contains the start 
    % time of the corresponding up state
    handles.USlist (1:nus,3) = num2cell(handles.upstates.toT(1:nus)); % then the end time
    handles.USlist (1:nus,4) = num2cell(handles.upstates.toT(1:nus) - handles.upstates.fromT(1:nus)); % and the length of the up state

    % Add the up state info to the table in the GUI
    set(handles.UStable,'Data',handles.USlist);
    
    zebraguidata(hObject, handles);

    % plot the state IDs
    plotStates (hObject, handles)

    % create the joined files
    handles.US = [];
    handles.DS = [];

    ifrom = handles.upstates.fromI(1);
    ito = handles.upstates.toI(1);
    handles.US = handles.LFP(ifrom:ito);

    for usi=2:nus
        ifrom = handles.upstates.fromI(usi);
        ito = handles.upstates.toI(usi);
        if get(handles.joinOffsetChk,'Value')
            offset = handles.LFP(ifrom)-handles.US(end); 
            handles.US = [handles.US handles.LFP(ifrom:ito)-offset];
        end    
        if get(handles.joinCk,'Value'), handles.US = [handles.US tempData(ifrom:ito)];
        end
    end

    for usi=1:nds
        ifrom = handles.downstates.fromI(usi);
        ito = handles.downstates.toI(usi);
        if get(handles.joinOffsetChk,'Value') %Gab
            if ~isempty(handles.DS) %Gab
                offset = handles.LFP(ifrom)-handles.DS(end); %Gab
                handles.DS = [handles.DS handles.LFP(ifrom:ito)-offset];    %Gab
            else %Gab
                handles.DS = [handles.DS handles.LFP(ifrom:ito)]; %Gab
            end %Gab
        else %Gab
        handles.DS = [handles.DS handles.LFP(ifrom:ito)];
        end
    end    

    % and finally lets compute the power spectra!
    handles.params = struct('tapers',[],'Fs',1,'fpass',[]);

    flagChronux = get(handles.parentHandles.ChronuxCk,'Value');
    if flagChronux
        params.tapers = [5 9] ;%[5 9];
        params.Fs = handles.parentHandles.acqF;
        params.fpass = [0 1000];
    %    zebraguidata(hObject, handles);
        [handles.pssUS,handles.fUS] = mtspectrumc(handles.US, params);
        handles.fUS = handles.fUS';
        [handles.pssDS,handles.fDS] = mtspectrumc(handles.DS, params);
        handles.fDS = handles.fDS';
    else
        [handles.pssUS,handles.fUS] = pwelch(handles.US,[],[],[],handles.parentHandles.acqF);
        [handles.pssDS,handles.fDS] = pwelch(handles.DS,[],[],[],handles.parentHandles.acqF);
    end
    %lUS = size(handles.fUS);
    %lDS = size(handles.fDS);

    axes(handles.powerSpectra);
    hold on;
    lim = axis;
    cla

    %flagN = get(handles.spectreData.normSpectraCk,'Value');
    % plot the mean traces
    plot(handles.fUS,10*log10(handles.pssUS.'),'red');            % spectra in decibels
    plot(handles.fDS,10*log10(handles.pssDS.'),'green');

    lim(1) = str2double(get(handles.minFreq,'String'));
    lim(2) = str2double(get(handles.maxFreq,'String'));
    lim(3) = str2double(get(handles.minPw,'String'));
    lim(4) = str2double(get(handles.maxPw,'String'));
    axis (lim);

    set(handles.powerSpectra,'XScale','log');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');

    zebraguidata(hObject, handles);
end
end

function analyzeData_Callback(hObject, eventdata, handles)
% This function computes the distribution of the spectral power in the
% 'envelope' band width. The distribution is fitted by a double gaussian
% and the limits for the event detection are computed

handles.parentHandles = zebraguidata(handles.parentObject);    % this to refresh the local values of parentHandles

handles.Hcenters = [];
handles.cdist = [];
handles.cdistN = [];
handles.Hnelements = [];
handles.normNelements = [];

% first, extract the data to be analyzed
dataSource = 9;
if get(handles.detectionSPG,'Value')
    dataSource = 0;
end
if get(handles.detectionBP1,'Value')
    if handles.parentHandles.BPcomputed(2), dataSource = 2;
    end
end
if get(handles.detectionBP2,'Value')
    if handles.parentHandles.BPcomputed(3), dataSource = 3;
    end
end
if get(handles.detectionBP3,'Value')
    if handles.parentHandles.BPcomputed(4), dataSource = 4;
    end
end
if get(handles.detectionHP,'Value')
    if handles.parentHandles.BPcomputed(5), dataSource = 5;
    end
end

if dataSource < 9
    ct = handles.parentHandles.currentTrial;
    cCh = handles.parentHandles.currentCh;
    if dataSource == 0
        sl = handles.parentHandles.spgl;
    else
        sl = handles.parentHandles.BPpowerL(dataSource);
    end
    logFlag = get(handles.powerLogCk,'Value');
    flagLeak = get (handles.parentHandles.displayLeakCk,'Value');
    leftH = str2double(get(handles.fromH, 'String'));
    rightH = str2double(get(handles.toH, 'String'));
    bins = str2double(get(handles.binNumber, 'String'));

    if (handles.parentHandles.deleakedFlag && flagLeak)
        % deleaked only for the SPG. 
        if logFlag
            temp (1:sl) = log10(handles.parentHandles.spgPlotDeleaked(cCh,ct,1:sl));
        else
            temp (1:sl) = handles.parentHandles.spgPlotDeleaked(cCh,ct,1:sl);
        end
    else
        if logFlag
            switch dataSource
                case 0      % SPG
                    temp (1:sl) = log10(handles.parentHandles.spgPlot(cCh,ct,1:sl));
                case {2, 3, 4, 5}      % BP1, BP2, BP3, HP 
                    temp (1:sl) = log10(handles.parentHandles.BPpower(2,dataSource,cCh,ct,1:sl));
            end
        else
            switch dataSource
                case 0      % SPG
                    temp (1:sl) = handles.parentHandles.spgPlot(cCh,ct,1:sl);
                case {2, 3, 4, 5}      % BP1, BP2, BP3, HP 
                    temp (1:sl) = handles.parentHandles.BPpower(2,dataSource,cCh,ct,1:sl);
            end
        end
    end

    % check about auto scaling or not
    if leftH==0, [nelements,centers] = hist(temp,bins);
    else
       xValues = leftH:(rightH-leftH)/bins:rightH;
       [nelements,centers] = hist(temp,xValues);
    end

    % compute the 1) normalized histo, 2) cumulative distribution, 3)
    % normalized cumulative distribution
    totalN = sum (nelements);
    handles.normNelements = nelements/totalN;
    for j=1:bins+1
        handles.cdist (j) = sum(nelements(1:j));
    end
    handles.cdistN = handles.cdist / totalN;

    handles.Hcenters = centers;
    handles.Hnelements = nelements;

    % plot the histogram
    axes (handles.histoPlot);
    cla
    bar(centers,nelements);
    axis ([leftH rightH 0 inf], 'auto y');
    handles.meanPW(cCh,ct) = mean(temp);
    handles.stdPW(cCh,ct) = std(temp);
    nsamples = length(centers);
    nmeasures = sum(nelements);

    tll = handles.meanPW(cCh,ct)-handles.stdPW(cCh,ct)/3;
    set(handles.baselineThr,'String',num2str(tll));
    tll = handles.meanPW(cCh,ct)+handles.stdPW(cCh,ct)/3;
    set(handles.eventThr,'String',num2str(tll));

    %f = fit(centers.',nelements.','gauss2');

    % fit the histogram with a double gaussian. 10.0266 = 2*2*sqr(2*pi)
    param(2)= handles.meanPW(cCh,ct)-handles.stdPW(cCh,ct); % mean first component
    param(3)= handles.stdPW(cCh,ct)/2;                      % SD first component
    amplitude(1)= nmeasures/(10.0266*param(3));                 % amplitude first component
    param(5)= handles.meanPW(cCh,ct)+handles.stdPW(cCh,ct); % mean second component
    param(6)= handles.stdPW(cCh,ct)/2;                      % SD second component
    amplitude(2)= nmeasures/(10.0266*param(6));                 % amplitude first component
    gaussFnc = [];

    % July 10, 2017. First, fit with fixed means and STDs to obtain a
    % preliminary estimates of the gaussian amplitudes.
    ampliOut = fminsearch(@gaussFit1,amplitude);
    param(1) = ampliOut(1);
    param(4) = ampliOut(2);
    % next, recompute everything.
    paramOut = fminsearch(@gaussFit,param);

    % paramOut contains the 6 parameters that define the best fit.
    st = ['Gauss 1; mean:' num2str(paramOut(2),'%2.2f') ' ±SD: ' num2str(paramOut(3),'%2.2f')];
    set (handles.gauss1param, 'String', st);
    st = ['Gauss 2; mean:' num2str(paramOut(5),'%2.2f') ' ±SD: ' num2str(paramOut(6),'%2.2f')];
    set (handles.gauss2param, 'String', st);

    hold on
    plot (centers,gaussFnc,'red','LineWidth',2);
    yl = ylim;
    line([paramOut(2)+paramOut(3) paramOut(2)+paramOut(3)],[yl(1) yl(2)]);
    line([paramOut(5)-paramOut(6) paramOut(5)-paramOut(6)],[yl(1) yl(2)]);

    % now compute the state thresholds
    xinterc = fzero(@gaussDiff,(paramOut(2)+paramOut(5))/2);

    handles.autoblThr = xinterc;
    handles.autoevThr = xinterc;
end

zebraguidata(gcbo, handles);

    function res = gaussFit(a)
        % gmodel = a(1)*exp(-(centers-a(2))^2/2*a(3).^2) + a(4)*exp(-(centers-a(5))^2/2*a(6).^2)
        % nelements is the array that contains the histogram bin counts: this
        % is the function to fit
        gaussFnc = a(1)*exp(-(centers-a(2)).^2/(2*(a(3)^2))) + a(4)*exp(-(centers-a(5)).^2/(2*(a(6)^2)));
        res = sum((nelements-gaussFnc).^2);
    end
    function res = gaussFit1(a)
        % gmodel = a(1)*exp(-(centers-a(2))^2/2*a(3).^2) + a(4)*exp(-(centers-a(5))^2/2*a(6).^2)
        % nelements is the array that contains the histogram bin counts: this
        % is the function to fit
        gaussFnc = a(1)*exp(-(centers-param(2)).^2/(2*(param(3)^2))) + a(2)*exp(-(centers-param(5)).^2/(2*(param(6)^2)));
        res = sum((nelements-gaussFnc).^2);
    end

    function y = gaussDiff(x)
        y = paramOut(1)*exp(-(x-paramOut(2)).^2/(2*(paramOut(3)^2)));
        y = y - paramOut(4)*exp(-(x-paramOut(5)).^2/(2*(paramOut(6)^2)));
    end
end

function binNumber_Callback(hObject, eventdata, handles)
end

function binNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function leftTime_Callback(hObject, eventdata, handles)
end

function leftTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function rightTime_Callback(hObject, eventdata, handles)
end

function rightTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function replot_Callback(hObject, eventdata, handles)
plotStates (hObject, handles)
end

function erodeCk_Callback(hObject, eventdata, handles)
end

function refineEdges(hObject, handles, tempData)

erodeFlag = get (handles.erodeCk,'Value');
if erodeFlag
    % now refine the edges of the states
    for usi=1:handles.NUS
        % check the edges of the US: if they are farther away from the median
        % than 2*SD the event is eroded to 2*SD. Check first left edge.
        edgeI = handles.upstates.fromI(usi);
        value = tempData(edgeI);
        while ~(abs(value-handles.upstates.medianV(usi))<0.5*handles.upstates.SD(usi)) 
            % erode to the right
            if edgeI==handles.parentHandles.dtaLen, break
            end
            edgeI = edgeI+1;
            handles.upstates.fromI(usi) = edgeI;
            value = tempData(edgeI);
            % check not to have passed the end of the US! 
        end
        edgeI = handles.upstates.toI(usi);
        value = tempData(edgeI);
        while ~(abs(value-handles.upstates.medianV)<2*handles.upstates.SD) 
            if edgeI>1
                % erode to the left
                edgeI = edgeI-1;
                handles.upstates.toI(usi) = edgeI;
                value = tempData(edgeI);
                % check not to have passed the beginning of the US!
            else
                break
            end    
        end    
        % now update the upstate definition
        handles.upstates.toT(usi) = (handles.upstates.toI(usi)-1)*spp;
        handles.upstates.fromT(usi) = (handles.upstates.fromI(usi)-1)*spp;
        handles.upstates.deltaT(usi) = handles.upstates.toT(usi)-handles.upstates.fromT (usi);
        
        handles.upstates.medianV(usi) = median(tempData(handles.upstates.fromI(usi):handles.upstates.toI(usi))); 
        handles.upstates.SD(usi) = std(tempData(handles.upstates.fromI(usi):handles.upstates.toI(usi))); 
    end
end

zebraguidata(hObject, handles);
end

function minPw_Callback(hObject, eventdata, handles)
end

function minPw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function maxPw_Callback(hObject, eventdata, handles)
end

function maxPw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function minFreq_Callback(hObject, eventdata, handles)
end

function minFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function maxFreq_Callback(hObject, eventdata, handles)
end

function maxFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function savePWS_Callback(hObject, eventdata, handles)

% save as two different files the PS computed from the defragmented UP and
% Downstate cumulative files.

handles.parentHandles = zebraguidata(handles.parentObject);    % this to refresh the local values of parentHandles

cTr = handles.parentHandles.currentTrial;
cCh = handles.parentHandles.currentCh;

fileOutUS = [handles.parentHandles.dir_in 'USpws Ch' num2str(cCh) ' Trial' num2str(cTr) '.dat'];
fileOutDS = [handles.parentHandles.dir_in 'DSpws Ch' num2str(cCh) ' Trial' num2str(cTr) '.dat'];
usl = length (handles.fUS);
tmpOut1 = [];
tmpOut2 = [];
% US and DS have to be saved as separate files since they have different x
% axis due to the variable length of the consolidated data.
tmpOut1 (1,:) = handles.fUS;        
tmpOut1 (2,:) = handles.pssUS;        
tmpOut1 (3,:) = 10*log10(handles.pssUS);        

tmpOut2 (1,:) = handles.fDS;        
tmpOut2 (2,:) = handles.pssDS;        
tmpOut2 (3,:) = 10*log10(handles.pssDS);        

% create the format descriptor
fmtSt=['%12.10f %12.10f %12.10f\n'];       % 3 columns and start new line
fid = fopen(fileOutUS,'w');
fprintf (fid,fmtSt,tmpOut1);
fclose(fid);
fid = fopen(fileOutDS,'w');
fprintf (fid,fmtSt,tmpOut2);
fclose(fid);

end

function removeCk_Callback(hObject, eventdata, handles)
end

function autoThr_Callback(hObject, eventdata, handles)
end

function saveSegmentedStates_Callback(hObject, eventdata, handles)
%GAB: saves 2 separate files for up & down states
currCh=num2str(handles.parentHandles.currentCh);
currTr=num2str(handles.parentHandles.currentTrial);
fmt='%8.6f \n';
fid=fopen([handles.parentHandles.dir_in 'DS_ch' currCh '_tr' currTr '.dat'],'w');
fprintf(fid,fmt,handles.DS);
fid=fopen([handles.parentHandles.dir_in 'US_ch' currCh '_tr' currTr '.dat'],'w');
fprintf(fid,fmt,handles.US);
fclose(fid);
end

function plotStates (hObject, handles)

len = handles.parentHandles.dtaLen;
spp = handles.parentHandles.sp;
handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));

axes(handles.plotStates);
cla
hold on
time = 0:spp:spp*(len-1);
plot (time,handles.LFP,'black');

for usi=1:handles.NUS
   i1 = handles.upstates.fromI(usi);
   i2 = handles.upstates.toI(usi);
   plot(time(i1:i2),handles.LFP(i1:i2),'red'); 
end

for usi=1:handles.NDS
   i1 = handles.downstates.fromI(usi);
   i2 = handles.downstates.toI(usi);
   plot(time(i1:i2),handles.LFP(i1:i2),'green'); 
end
axis ([handles.tmin handles.tmax -inf inf])

axes(handles.statesTrack);
set(gca, 'XTick', []);
cla
hold on
plot (handles.timeSPG,handles.SPG);
axis ([handles.tmin handles.tmax -inf inf])
hold off

end

function traslateLeft_Callback(hObject, eventdata, handles)
handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));
delta = handles  .tmax - handles.tmin;

set(handles.leftTime,'String',handles.tmin-delta);
set(handles.rightTime,'String',handles.tmin);
zebraguidata(hObject, handles);

plotStates (hObject, handles)
end

function traslateRight_Callback(hObject, eventdata, handles)
handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));
delta = handles.tmax - handles.tmin;

set(handles.leftTime,'String',handles.tmax);
set(handles.rightTime,'String',handles.tmax+delta);
zebraguidata(hObject, handles);

plotStates (hObject, handles)
end

function saveHCD_Callback(hObject, eventdata, handles)
% Save the distributions of the integrated spectral power
% March 22, 2017.

tmpOut1 = [];

handles.parentHandles = zebraguidata(handles.parentObject);    % this to refresh the local values of parentHandles
cTr = handles.parentHandles.currentTrial;
cCh = handles.parentHandles.currentCh;

fileOut = [handles.parentHandles.dir_in 'PWSdistribution Ch' num2str(cCh) ' Trial' num2str(cTr) '.dat'];
tmpOut1 (1,:) = handles.Hcenters;        
tmpOut1 (2,:) = handles.Hnelements;    
tmpOut1 (3,:) = handles.normNelements;
tmpOut1 (4,:) = handles.cdist;    
tmpOut1 (5,:) = handles.cdistN;

% create the format descriptor
fmtSt=['%12.10f %12.10f %12.10f %12.10f %12.10f\n'];  % 5 columns and start new line
fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut1);
fclose(fid);

end

function toH_Callback(hObject, eventdata, handles)
end

function toH_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function fromH_Callback(hObject, eventdata, handles)
end
   
function fromH_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function detectionBP1_Callback(hObject, eventdata, handles)
end

function detectionBP2_Callback(hObject, eventdata, handles)
end

function detectionBP3_Callback(hObject, eventdata, handles)
end

function detectionHP_Callback(hObject, eventdata, handles)
end

function detectionSPG_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in SaveDownstSpec.
function SaveDownstSpec_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDownstSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveDownstSpectra(handles)
end
