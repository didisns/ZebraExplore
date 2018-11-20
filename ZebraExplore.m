function varargout = ZebraExplore(varargin)

% Zebra Explore: November 2016
% Visualisation of long records of EF data for the evaluation and classification of
% spontaneous asynchronous activity.

% ZEBRAEXPLORE MATLAB code for ZebraExplore.fig
%
%      H = ZEBRAEXPLORE returns the handle to a new ZEBRAEXPLORE or the handle to
%      the existing singleton*.
%
%      ZEBRAEXPLORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZEBRAEXPLORE.M with the given input arguments.
%
%      ZEBRAEXPLORE('Property','Value',...) creates a new ZEBRAEXPLORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ZebraExplore_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ZebraExplore_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ZebraExplore

% Last Modified by GUIDE v2.5 20-Nov-2018 11:26:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ZebraExplore_OpeningFcn, ...
                   'gui_OutputFcn',  @ZebraExplore_OutputFcn, ...
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


% --- Executes just before ZebraExplore is made visible.
function ZebraExplore_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ZebraExplore (see VARARGIN)

% Choose default command line output for ZebraExplore
handles.output = hObject;

% Initialize variable and parameters and stuff them in the handles
% structure

handles.sp = str2double(handles.spEdit.String);            % sampling period in s
handles.acqF = 1/handles.sp;    % sampling frequency
handles.insetW = 5;             % width of the inset window (s)
handles.lp = 0.1;               % width of the low pass filter window % MIGHT BE UNUSED
handles.dir_in = '';
handles.gainAmp = str2double(handles.gain.String);         % gain setting of the amplifier
handles.nCh = str2double(handles.chEdit.String);                % number of channels
handles.nTrials = str2double(handles.trialsEdit.String);            % number of trials
handles.currentCh = handles.chSlider.Value;          % current channel

handles.filterOrder = str2double(handles.filtOrder.String);        % butterworth filter order
handles.notch = [str2double(handles.notchFrom.String) str2double(handles.notchTo.String)];        % notch filter for 50 Hz removal
handles.frequencyBand = [0 str2double(handles.lpTo.String);...
    str2double(handles.alphaFrom.String) str2double(handles.alphaTo.String);...
    str2double(handles.betaFrom.String) str2double(handles.betaTo.String); ...
    str2double(handles.gammaFrom.String) str2double(handles.gammaTo.String);...
    str2double(handles.hpFrom.String) 0.5/handles.sp;...
    0 str2double(handles.leakageLP.String)];
handles.notchOrder = 3;         % order 6 for the notch filter % MIGHT BE UNUSED
handles.xnow = -999;

handles.fileInNum = 0;      % keep track of the number of analyzed files for the output table
handles.eventNum = 0;       % keep track of the number of analyzed events
handles.cnt = 0;            % counter for full trace summary table
handles.EPcnt = 0;          % counter for the summary table of evoched responses

% the following structure is required for passing the proper parameters to
% the Chronux routines.
handles.params = struct('tapers',[1 2],'Fs',0.001,'fpass',[0 10],'pad',0,'err',[]);

% handles.dtaSummary = struct('f1','file name','f2',1,'f3',2,'f4',3,'f5',4,'f6','Comments');

% Definition of a cell array to store the output table
% The fields are consistent with te GUIDE definition of outable provided in
% the ZebraSays figure.
% dtaSummary holds the complete file name time of file last modification
% (it corresponds to the original creation time) and the RMS powers in the 5 bands 

handles.dtaSummary = {'dir' 'file' 'time' 1 2 3 4 5 6 7 'notes'};

handles.eventSummary = {'dir' 'file' 1 2 3 4 5 6 7 8 9 10 11 12 13 'news'}; 
% Cell array for the detected events
% Field description: 
% 1: event number in the given file
% 2: time at the centre of the event
% 3: duration of the event
% 4-6: RMS power of the baseline in the three bands
% 7-9: RMS power of the event
% 10-12: power normalised to the baseline
% 
% 14: news and views (editable text)

handles.respSummary = {'file' 1 2 3 4 5 6 7 8 9 10 'news'}; 
% Cell array for trial responses: field description
% 1: file name
% 2: channel number
% 3: trial number
% 4: repeat number
% 5: Peak response
% 6: time to peak
% 7: Mean BP1 baseline power (Db)
% 8: Mean BP1 response power (Db)
% 9: normalized power (8-7)
% 10: BP1 power multiplied by the scaled template
% 11: comments

handles.dtaPWsummary = {'dir' 'file' 'time' 4 5 6 7 8 9 10 11 12 13};
handles.dtaPWlocalSummary = {'dir' 'file' 'time' 4 5 6 7 8 9 10 11 12 13};
% Cell array for power distribution: field description
% 1: file name
% 2: channel number
% 3: time of data file creation
% 4: channel
% 5: trial
% 6: band pass index (1-4 corresponding to BP1-3 and HP
% 7: Shannon entropy 
% 8: main mode: mean log of power 
% 9: main mode: std
% 10: main mode: relative area
% 11: secondary mode: mean of log of power
% 12: secondary mode: std 
% 13: secondary mode: relative area

% handles to the external figures
handles.outWnd = ZebraSays;      % handle to the output window
handles.zebraSaysData = zebraguidata(handles.outWnd);    % and its controls

handles.evDetection = eventDetection(hObject,handles);
handles.evDetData = zebraguidata (handles.evDetection);

%GAB
handles.autoMode=0;
handles.storedLFP = [];

% Update handles structure
zebraguidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ZebraExplore_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

function select_data_Callback(hObject, eventdata, handles)

% Open new data file according to the selected data type

% Check if any accessory window is selected and, if that is the case,
% closes it.

if get(handles.mainVisualizer, 'Value') == 0,
    set (handles.mainVisualizer, 'Value',1);
    clearMainWindow (hObject, handles, 'mainVisualizer')
end

flagN = get(handles.normalDataRB,'Value');
if (flagN), fmtSt = '*.dat';
end
flagL = get(handles.lauraEEG,'Value');
if (flagL), fmtSt = '*.eeg';
end
flagMV = get(handles.meyerVepRB,'Value');
if (flagMV), fmtSt = '*.asc';
end
flagSSP = get(handles.meyerSspRB,'Value');
if (flagSSP), fmtSt = '*.asc';
end
flagAxP = get(handles.axopatchCk,'Value');
if (flagAxP), fmtSt = '*.abf';
end
flagTRC = get(handles.TRCimportRB,'Value');
if (flagTRC), fmtSt = '*.trc';
end
flagMEA = get(handles.meaRB,'Value');
if (flagMEA), fmtSt = '*.mat';
end


if handles.autoMode    %:::GABRIELE::: 08/12/2017, for AutoOpening
    selection=handles.selection;
    path_in=[handles.path_in '\'];
    zebraguidata(hObject,handles);
    handles.fileTime=[]; % needed by OpenFile
    openFile(hObject,handles,selection,path_in)
    handles=zebraguidata(hObject);
else
    [selection,path_in,filterI] = uigetfile ([handles.dir_in fmtSt]);
    fileInfo = {'name' 'date' 'bytes' 'isdir' 'datenum'};
    zebraguidata(hObject, handles);
    if selection ~= 0
        fileInfo = dir([path_in selection]);
        timeNum = datestr(fileInfo.datenum);
        handles.fileTime = timeNum (13:20);
        zebraguidata(hObject, handles);
        openFile(hObject,handles,selection,path_in)
        handles = zebraguidata(hObject);
    end
end
zebraguidata(hObject, handles);

function openDir_Callback(hObject, eventdata, handles)
flagN = get(handles.normalDataRB,'Value');
if (flagN), fmtSt = '*.dat';
end
dirArray = {'name' 'date' 'bytes' 'isdir' 'datenum'};
fileInfo = {'name' 'date' 'bytes' 'isdir' 'datenum'};
selection = uigetdir('C:\Users\Gix\Documents\Work\Science\')
if selection ~= 0
    fileName = '\ltp001.dat';
    dirArray = dir([selection '\Dir_*']);
    ndir = length(dirArray);
    for i=1:ndir
        path_in = [selection '\' dirArray(i).name];
        fileInfo = dir([path_in fileName]);
        timeNum = datestr(fileInfo.datenum);
        handles.fileTime = timeNum (13:20);
        zebraguidata(hObject, handles);
        openFile(hObject,handles,fileName,path_in)
        handles = zebraguidata(hObject);
    end    
end
zebraguidata(hObject, handles);

function openFile (hObject, handles, selection, path_in)

% this function open and process the file specified in [path_in selection]

%%
handles.file_in = selection;
handles.dir_in = path_in;
handles.dtaComplete=[path_in selection];
set(handles.plotSlate,'title',handles.dtaComplete);

% Output to console of status and timing
tic
disp([num2str(toc) ' Begin to load data file'])
% initialize all arrays
LFPin = [];
handles.LFP = [];
handles.bandPassed_LFP = [];
handles.workLFP = [];
handles.spg = [];
handles.spgDeleaked = [];
handles.spgPlot = [];
handles.spgPlotDeleaked = [];
handles.meanLFP = [];
handles.meanBP = [];
handles.meanSpg = [];
handles.meanSpgPlot = [];
handles.meanPWS = [];
handles.EIrat = [];
handles.meanEIrat = [];
handles.EPamplitude = [];
handles.EIresponse = [];
handles.EPoffset = [];
handles.templ = [];

% flags to indicate the frequency domain computations that have been
% performed on the data set
handles.BPcomputed(1:6) = 0;
% 1 = LP; 2-4 = BP1-3; 5 = HP; 6 = LP for leakage

handles.EPcnt = 0;
handles.xnow = -100;
exCh = str2double(get(handles.excludeCh,'String'));
sCh = str2double(get(handles.synchCh,'String'));
handles.deleakedFlag = 0;       % this flag is set to 1 if SPGs have been deleaked.
zebraguidata(hObject, handles);      % refresh the content of handle

% Now the processing pipelines differentiate depending on the input mode.
flagN = get(handles.normalDataRB,'Value');
flagL = get(handles.lauraEEG,'Value');
flagMV = get(handles.meyerVepRB,'Value');
flagSSP = get(handles.meyerSspRB,'Value');
flagAxP = get(handles.axopatchCk,'Value');
flagTRC = get(handles.TRCimportRB,'Value');
flagMEA = get(handles.meaRB,'Value');
flagEx = get(handles.excludeChFlag,'Value');
flagExSynch =  get(handles.internalSynchCk,'Value');
flagSkipFirstTrial = get(handles.skip1trialCk,'Value');

% This function opens a specified file either in ASCII or in binary format.
% the variable nCh holds the number of channels present in the data file
% the resulting LFP data will be a matrix of nSweeps columns.
%%

%% TRC data import
if flagTRC
    % open the dialogue window for TRC data import
    % The TRCimport must be opened by passing the file path and it returns the structure TRC that contina
    % the data, all the associated info and the events.
    TRCoutput = cell([40 1]);
    [TRCoutput{1:40}] = importTRCfilesOptions(handles.dir_in,handles.file_in,handles.dtaComplete);
    disp([num2str(toc) ' End of TRC import'])
    %TRCimportData = zebraguidata (TRCimportDialogue);
    % let's fill the data structure with the imported data
    set (handles.microVck,'Value',1);   % select microV unit
%         setGain (hObject, handles, 1000);
%         handles = zebraguidata(hObject);
    handles.expID = TRCoutput{1};
    nCh = TRCoutput{9};
    handles.SPGlowDB = -80;                             % set limits for spectrogram plot
    set(handles.lowerDB,'String',handles.SPGlowDB);
    handles.SPGhighDB = 0;
    set(handles.upperDB,'String',handles.SPGhighDB);
    
    setCh (hObject, handles, nCh);
    handles = zebraguidata(hObject);
    handles.dtaLen = TRCoutput{11};
    acqF = TRCoutput{12};
    setSP (hObject, handles, acqF);
    handles = zebraguidata(hObject);
    handles.nTrials = length(TRCoutput{26}) - 1;     % if no events nTrials = -1 !
    LFPin = zeros(handles.nCh, handles.dtaLen);      % preallocation of LFPin
    LFPin = TRCoutput{16};
    if (handles.nTrials>1)
        set (handles.skip1trialCk,'value',0);       % no need to skip the first trial
        stimTrigger (1:handles.nTrials) = 0;        % preallocation of trigger
        for j=1:handles.nTrials
            stimTrigger (j) = TRCoutput{26}(j+1).latency; 
            % the first trigger is skipped because it is associated to the
            % beginning of the equiluminant grey and not to a chekerboard stimulus.
        end
        %Compute the inter-trial distance.
        trialInterdistance (1:handles.nTrials-1) = stimTrigger (2:handles.nTrials)-stimTrigger (1:handles.nTrials-1);
        % Lets find the minimal inter trial distance
        minLen = min(trialInterdistance);      % this will be the new handles.dtaLen
        % Let's extract the trials from the available data channels.
        disp([num2str(toc) ' Extract trials from the data channels.'])
        for (i=1:handles.nCh)
            for (j=1:handles.nTrials)
                temp(i,j,1:minLen) = LFPin (i,stimTrigger(j):stimTrigger(j)+minLen-1);
            end    
        end
        handles.LFP = [];
        handles.LFP = temp;
        handles.dtaLen = minLen;
    else
        handles.nTrials = 1;
        handles.LFP(1:handles.nCh,1,1:handles.dtaLen) = LFPin(1:handles.nCh,1:handles.dtaLen);
    end
else
    % create an handle for the file opened in read mode
    fid=fopen(handles.dtaComplete,'r');
end
%%

%% Import NI or Axiopatch data
if (flagN || flagL || flagAxP)
    % This is the default mode. Data from the NI or Axon board. nCh channels and only 1
    % trial. The trigger channel is included. The number of channels must
    % be defined at priori only for the NI boards
    set (handles.mVck,'Value',1);
    if (flagAxP)
        % in this modality data are generated by the Axopatch program and
        % are saved as binary data. sp and Ch number obtained from data header
        disp([num2str(toc) ' Begin loading of Axopatch data'])
        [LFPin,si,h] = abfload(handles.dtaComplete);
        % si    scalar           the sampling interval in micro s
        % h     struct           information on file (selected header parameters)
        % number of channels
        LFPin = LFPin';    % abfloat return the channels as columns. Traspose to set them as rows
        setGain (hObject, handles, 1000);
        handles = zebraguidata(hObject);
        setSP (hObject, handles, 1/si*1000000);
        handles = zebraguidata(hObject);
        setCh (hObject, handles, h.nADCNumChannels);
        handles = zebraguidata(hObject);
        handles.dtaLen = h.sweepLengthInPts;
    else    
        if flagL    % set some parameters specific to Laura's mode.
            setSP (hObject, handles, 512);
            handles = zebraguidata(hObject);
            setGain (hObject, handles, 5000);
            handles = zebraguidata(hObject);
            setCh (hObject, handles, 2);
            handles = zebraguidata(hObject);
            set(handles.synchCh, 'String', '1');
            set(handles.SPGwindow,'String','16');
            set(handles.SPGoverlap,'String','2');
        else
            set(handles.SPGwindow,'String','256');
            set(handles.SPGoverlap,'String','16');
            handles.nCh = str2double(get(handles.chEdit,'String'));
        end
        if flagL
            disp([num2str(toc) ' Begin loading of IN EEG data.'])
            % this is required because these data use a ',' as separator of the
            % decimal part. How stupid is this???
            c = textscan(fid,'%s');
            l = size(c{1});
            handles.dtaLen = l(1)/2-1;  % remove first row
            LFPin = zeros(handles.nCh,handles.dtaLen);
            for i=1:handles.dtaLen
                % replace the comma with the decimal point
                for j=1:handles.nCh
                    runningIndex = 2 + handles.nCh * (i-1) + j;
                    LFPin (j,i) = str2double(strrep(c{1}{runningIndex},',','.'));
                end 
            end
        else
            disp([num2str(toc) ' Begin loading of National Instruments EEG data.'])
            % Create the string describing the data format
            fmtSt='';
            for ii=1:handles.nCh
               fmtSt=[fmtSt '%f'] ;
            end
            LFPin = fscanf(fid,fmtSt,[handles.nCh,inf]);
            l = size(LFPin);
            handles.dtaLen = l(2);
        end
    end
    %% the three submode of DEFAULT mode continue from here
    handles.nTrials = 1;    % in this open mode there is only one trial unless
    % the trials are extracted later on by using the trigger channel 
    handles.SPGlowDB = -80;
    set(handles.lowerDB,'String',handles.SPGlowDB);
    handles.SPGhighDB = 0;
    set(handles.upperDB,'String',handles.SPGhighDB);
    if flagEx && ~flagExSynch    % if ext synch the channel is excluded only from the output
        if (handles.nCh<=1),handles.nCh=2;   % This situation occurs if the original data has 2 channels 
        end                                  % and one forgets to restore the ch number
        % remove one channel
        cnt = 0;
        for i=1:handles.nCh
            if (i~= exCh), cnt = cnt+1;
                handles.LFP (cnt,1,:) = LFPin(i,:);
            end
        end    
        handles.nCh = handles.nCh-1;
    else
        handles.LFP (1:handles.nCh,1,:) = LFPin(1:handles.nCh,:);
    end
    if flagExSynch 
       disp([num2str(toc) ' Extract single trials from a multitrial data set.']) 
       maxTrigger = str2double(get(handles.triggerMaxLength,'String')); 
       indexTrials = [];
       trialInterdistance = [];
       temp = [];
       % first, analyze the synch channel to extract the timing of each trial
       % binarize the trigger channel
       for i=1:handles.dtaLen
           if (handles.LFP(sCh,1,i)<1), handles.LFP (sCh,1,i) = 0;
           else
               handles.LFP (sCh,1,i) = 10;
           end
       end
       [sortedSynch,sortedInd] = sort(handles.LFP (sCh,1,:),'descend');
       % the indexes should be clusterized and the <t> of each cluster will
       % give the timing of each stimulus. What is left is easy.
       % Parse the sortedInd vector. A trigger event is represented by a cluster of
       % neighboring indexes. 
       j = 1;   % index running on the sorted arrays
       flg = get(handles.skip1trialCk,'Value');
       trialCounter = 0;
       sortI = sortedInd(j);
       oldSortI = 1;
       oldJ = 0;
       while (sortedSynch(j)>1)   
           % The code handles two different trigger signals
           % The <yy 2018 implementation of Escher generates the trigger through
           % the DTR line of the RS232 serial port with levels on= +10V off= -10V 
           % Later build wil use a digital output producing the ordinary TTL signal
           
           % Legality of the trigger event
           % A trigger event is legal if it meets 3 conditions:
           % 1) it is longer than 30 sampling periods
           % 2) it is at least 100 sp far away from the previous trigger event
           % 3) it is shorter than triggerMaxLength (still in SP units)
           
%           if (sum(sortedSynch(j:j+2))>29 && sortI>(oldSortI+100))
           if (sortI>(oldSortI+100) || trialCounter==0)
               % if this condition is true the new event is at least 100 sp away
               % from the previous event, or it is the first event
               % now we must measure its length
               len = 1;
               runningJ = j;
               while sortedInd(runningJ+1) == sortedInd(runningJ)+1
                   % we are counting the SPs in the event
                   len = len + 1;
                   runningJ = runningJ+1;
               end
               % now test the event for legality
               if (len>2 && len<maxTrigger)
                   % event is legal
                   trialCounter = trialCounter + 1;
                   indexTrials (trialCounter) = sortI;
                   if (trialCounter>1)
                       if flg          % if true ignore the first trial
                           if (trialCounter>2)
                               % ignore the first trial during the computation
                               % of the trial duration
                               trialInterdistance(trialCounter-2) = sortI-indexTrials (trialCounter-1);
                           end
                       else
                           trialInterdistance(trialCounter-1) = sortI-indexTrials (trialCounter-1);
                       end
                   end
                   oldSortI = sortI;
               end
               j = runningJ;
           end
           j = j+1;
           sortI = sortedInd(j);
       end
       % Lets find the minimal inter trial distance
       minLen = min(trialInterdistance); % this will be the new handles.dtaLen
       % Let's extract the trials from the available data channels.
       trTimeOffset = str2double(handles.trTimeOffset_txt.String); %offset in seconds GAB 30/10/2018
       trTimeOffset_i = round(trTimeOffset/handles.sp); 
       for i=1:handles.nCh
            skipped = 0;
            for j=1:trialCounter
                if (indexTrials(j)+minLen-1)>handles.dtaLen
                    % bad synch: the data finishes before of the end of the
                    % last sweep
                    trialCounter = trialCounter - 1;
                else
                    if (indexTrials(j)+trTimeOffset_i) > 0
                        temp(i,j-skipped,1:minLen) = handles.LFP (i,1,indexTrials(j)+trTimeOffset_i:indexTrials(j)+minLen-1+trTimeOffset_i); %GAB 30/10/2018 added +trTimeOffset_i
                    else
                        skipped = skipped + 1;
                    end
                end
            end    
       end
       handles.LFP = [];
       handles.LFP = temp;
       handles.nTrials = trialCounter-skipped;
       handles.expID = '';
       handles.dtaLen = minLen;
    end
    % Close the data file
    fclose(fid);
end
%%

%% Import Meyer ascii VEP
if (flagMV)
    disp([num2str(toc) ' Begin loading of Meyer ASCII data.'])
    setSP (hObject, handles, 2048/str2double(get(handles.SweepTimeEdit,'string'))); %GABRIELE 2018/01/27; frequency is calculated as trialPoints/trialTime
    % since the SP has been likely changed, the frequencies of the
    % bandwidht filters should be taken care of
    handles = zebraguidata(hObject);
    setGain (hObject, handles, 1000);
    set (handles.microVck,'Value',1);
    handles = zebraguidata(hObject);
    %setCh (hObject, handles, 3);
    numCh=str2double(get(handles.chEdit,'String')); %GABRIELE 2018/01/27; allows to open file whith any number of channels
    setCh (hObject, handles, numCh);
    handles = zebraguidata(hObject);
    handles.dtaLen = 2048;
    set(handles.SPGwindow,'String','32');
    set(handles.SPGoverlap,'String','2');
    c = textscan(fid,'%s');
    handles.expID = [c{1}{1} ' ' c{1}{2}];
    % create an array with the indexes of all cells that contain the string
    % '2048'. The data will begin immediately after the last cell pointed by the 
    % index array. The length of the array is equal to the number of trials*3. 
    index = find(strcmp(c{:}, '2048'));
    l = size(index);
    handles.nTrials = l(1)/numCh;
    startC = index(handles.nTrials*numCh)+1;
    [endC,g] = size(c{:});

    % loop on all points of each trial
    for (j=1:handles.dtaLen)    %GABRIELE 2018/01/27; allows to open file whith any number of channels
        % loop on all trials
        for i=1:handles.nTrials
            % run on the entire cell array
            runningIndex = startC + (j-1)*handles.nTrials*numCh + (i-1)*numCh;
            % first, replace all the ',' with the decimal point
            for k=1:numCh
                handles.LFP (k,i,j) = str2num(strrep(c{1}{runningIndex+k-1},',','.'));
            end
        end    
    end
    % Close the data file
    fclose(fid);
end
%%

%% Load somatosensory potentials
if (flagSSP)
    % Better be safe...
    setSP (hObject, handles, 32768);
    handles = zebraguidata(hObject);
    setGain (hObject, handles, 1000);
    handles = zebraguidata(hObject);
    set (handles.microVck,'Value',1);
    c = textscan(fid,'%s');
    handles.expID = [c{1}{1} ' ' c{1}{2}];
    % create an array with the indexes of all cells that contain the string
    % '1639'. The data will begin immediately after the last cell pointed by the 
    % index array. The length of the array is equal to the number of trials*3. 
    index = find(strcmp(c{:}, '1639'));
    l = size(index);
    handles.nCh = l(1);
    handles.nTrials = 1;
    handles.dtaLen = 1639;

    % startC point to the beginning of the data
    startC = index(handles.nCh)+1;
    [endC,g] = size(c{:});

    % loop on all points of each trial
    for (j=1:handles.dtaLen)
        % loop on all channels
        for (i=1:handles.nCh)
            % run on the entire cell array
            runningIndex = startC + (j-1)*handles.nCh + (i-1);
            handles.LFP (i,1,j) = str2double(strrep(c{1}{runningIndex},',','.'));
        end    
    end
    % Now remove the artefact that sometimes pops out at the very beginning of the record
    clsW = int32(handles.dtaLen/10);
    fixIt(1:handles.nCh) = handles.LFP(1:handles.nCh,1,clsW);
    for (k=1:handles.nCh)
        handles.LFP(k,1,1:clsW) = fixIt(k);
    end    
    handles.workLFP = handles.LFP;    % Copied here since there is nothing to filter in these data    
    % Close the data file
    fclose(fid);
end
%%

%% Open IMM mea data
if (flagMEA)
    disp([num2str(toc) ' Begin loading of IMM MEA data.'])
    setGain (hObject, handles, 1000);
    handles = zebraguidata(hObject);
    % the data are contained in a structure called channelPlotData
    load(handles.dtaComplete)
    handles.nCh = 32;
    set (handles.chEdit,'String','32');
    
    name = fieldnames(channelPlotData);
    % extract the sampling period from channel0 data structure
    sp = channelPlotData.channel0.time(1,2);
    % NOTE: this SP is only approximate. THIS SHOULD BE FIXED!
    sp = 0.0001003394;
    zebraguidata(hObject,handles);
    setSP (hObject, handles, 1/sp);
    handles = zebraguidata(hObject);
    % preallocate the input array
    fullDataSet = zeros(handles.nCh,1,length(channelPlotData.(name{1}).value));
    % load all the channels
    for i=1:handles.nCh
        % run on the entire cell array
        tmp = channelPlotData.(name{i}).value;
        fullDataSet (i,1,:) = tmp(:);
    end
    
    if flagExSynch
        disp([num2str(toc) ' Extract single trials from a multi-trials data set.'])
        % extract the single trials
        tmp = [];
        stimOrigin = channelPlotData.Stimulation.time(1,6);
        traceOrigin = channelPlotData.channel0.time(1,1);
        stimTimes = channelPlotData.Stimulation.time(:,6) - traceOrigin;
        handles.nTrials = length(stimTimes);
        intervals = stimTimes(2:handles.nTrials)-stimTimes(1:handles.nTrials-1);
        handles.dtaLen = floor(min(intervals)/handles.sp + 1);
        pnt = handles.dtaLen;
        % prealocate LFP
        handles.LFP = zeros(handles.nCh,handles.nTrials,pnt);
        % now extract the data segments
        for j=1:handles.nTrials
           beginIndex = floor(stimTimes(j)/handles.sp + 1); 
           for i=1:handles.nCh
                handles.LFP(i,j,1:pnt) = fullDataSet (i,1,beginIndex:beginIndex+pnt-1);
                % filter the initial and final transient (trigger artifact)
                handles.LFP(i,j,1:10) = handles.LFP(i,j,20:-1:11);
                handles.LFP(i,j,pnt:-1:pnt-30) = handles.LFP(i,j,pnt-30:-1:pnt-60);                
           end
        end
    else
        handles.LFP = fullDataSet;
        handles.nTrials = 1;
        handles.dtaLen = length(handles.LFP);
    end
    handles.expID = '';
end

%% gab&Enrico trial concatenation
if handles.concat_ck.Value
    if isempty(handles.storedLFP) %ricordiamoci di inizializzarlo nella OpeningFunction
        handles.storedLFP = handles.LFP(:,flagSkipFirstTrial+1:end,:);
    else
        storeDim = size(handles.storedLFP);
        if storeDim(1) ~= handles.nCh %quality control
            error('Mismatch in number of channels')
        end
        storeLen = min([storeDim(3), handles.dtaLen]);
        handles.storedLFP = cat(2, handles.storedLFP(:,:,1:storeLen), handles.LFP(:,flagSkipFirstTrial+1:end,1:storeLen));
    end
else
    handles.storedLFP = [];
end

if handles.analyzeWhileConc_ck.Value
    concatenateTrials(hObject,handles);
    handles = zebraguidata(hObject);
elseif ~handles.concat_ck.Value
    openFile_partII(hObject,handles);
    handles = zebraguidata(hObject);
end
%% end of specific loading routines
% -----------------------------------------------------------------------------------------
zebraguidata(hObject, handles)
disp([num2str(toc) ' End of data loading.'])



function openFile_partII(hObject,handles) %Gab&Enrico 2018/11/02

%%--flag copied from OpenFile function-----------------------------------
exCh = str2double(get(handles.excludeCh,'String'));
sCh = str2double(get(handles.synchCh,'String'));

% Now the processing pipelines differentiate depending on the input mode.
flagN = get(handles.normalDataRB,'Value');
flagL = get(handles.lauraEEG,'Value');
flagMV = get(handles.meyerVepRB,'Value');
flagSSP = get(handles.meyerSspRB,'Value');
flagAxP = get(handles.axopatchCk,'Value');
flagTRC = get(handles.TRCimportRB,'Value');
flagMEA = get(handles.meaRB,'Value');
flagEx = get(handles.excludeChFlag,'Value');
flagExSynch =  get(handles.internalSynchCk,'Value');
%%-----------------------------------------------------------------------

% now set the values for the ch and trial sliders based on the values
% defined in the preceding code. If ch and trials are = 1 the control is
% not updated and is made not visible.
set (handles.trialsEdit,'String',handles.nTrials);
handles.currentTrial = 1;
%set (handles.chEdit,'String',handles.nCh);
handles.currentCh = 1;
if (handles.nTrials>1)
    set (handles.trialSlider, 'Max', handles.nTrials);
    set (handles.trialSlider, 'Value', 1);
    set (handles.trialSlider, 'SliderStep', [1/(handles.nTrials-1) 1/(handles.nTrials-1)]); %GAB 06/04/2018
    set (handles.trialSlider, 'Visible', 'on');
    set (handles.trialText, 'Visible', 'on');    
    set(handles.trialText, 'String', ['Trial ' num2str(handles.currentTrial)]);
else
    set (handles.trialSlider,'Visible', 'off');
    set (handles.trialText,'Visible', 'off');    
end

if (handles.nCh>1)
    set (handles.chSlider,'Max',handles.nCh);
    set (handles.chSlider, 'Value', 1);
    set (handles.chSlider, 'SliderStep', [1/(handles.nCh-1) 1/(handles.nCh-1)]);
    set (handles.chSlider,'Visible','on');
    set (handles.ChText,'Visible','on');    
    set (handles.ChText,'String',['Ch ' num2str(handles.currentCh)]);
    % set the slider control in the PWpower distribution figure
    if get(handles.analyzePWdistributionCk,'Value')
        set (handles.analyzePWdata.selectCh,'Max',handles.nCh);
        set (handles.analyzePWdata.selectCh, 'Value', 1);
        set (handles.analyzePWdata.selectCh, 'SliderStep', [1/(handles.nCh-1) 1/(handles.nCh-1)]);
        set (handles.analyzePWdata.selectCh,'Visible','on');
        set (handles.analyzePWdata.selectedChTxt,'Visible','on');    
        set (handles.analyzePWdata.selectedChTxt,'String',['Ch ' num2str(handles.currentCh)]);
    end
else
    set (handles.chSlider,'Visible','off');
    set (handles.ChText,'Visible','off');    
    if get(handles.analyzePWdistributionCk,'Value')
        set (handles.analyzePWdata.selectedChTxt,'Visible','off');
        set (handles.analyzePWdata.selectCh,'Visible','off');    
    end
end

zebraguidata(hObject, handles);

handles.fileInNum = handles.fileInNum + 1;
zebraguidata(hObject, handles);
drawnow    % refresh the GUI
cfactor = 1000/handles.gainAmp;      % conversion factor to express LFP in mV

% clear the inset windows
axes(handles.inset_WB);
cla
axes(handles.inset_HP);
cla
axes(handles.spectraInset);
cla

%dtaLen = size(handles.LFP);
%handles.dtaLen = dtaLen(3);


handles.x = 0:handles.sp:double(handles.sp*(handles.dtaLen-1));
handles.LFP = handles.LFP * cfactor;            % scale the LFP for the relative gain
% set t axes
handles.tmin = 0;
handles.tmax = handles.dtaLen*handles.sp;
set(handles.tminTxt,'String',handles.tmin);
set(handles.tmaxTxt,'String',handles.tmax);
zebraguidata(hObject, handles);

% break point here if you want to access the loaded data before computing the
% band-passed data set.
if (~flagSSP)
    disp([num2str(toc) ' Begin computation of band passed data.'])
    computeBPandNotch (hObject, handles);   % elements of handles will be modified at this call
    handles = zebraguidata(hObject);
    computePowerInTime (hObject, handles);  % temporal profile of power
end

% here I have learned an important lesson: handles are passed by values:
% any changes on the structure operated in called function IS NOT returned
% back to the calling function. If the changes are saved to the zebraguidata
% object, the new content of the object has to be loaded in order to update
% the handles.

handles = zebraguidata(hObject);    % this to refresh the local values of handles

disp([num2str(toc) ' Compute the data spectograms.'])
computeSpectrogram (hObject, handles);
handles = zebraguidata(hObject);    % this to refresh the local values of handles

% ---------------------------------------------------------------------------------------------------------------------

disp([num2str(toc) ' Passing data to the active external processor.'])
% now verify whether the external processors are active and, if they are, call
% them up passing them the proper data.

if get(handles.powerSpectraCk, 'Value')
    % call the ZebraSpectre code for spectra computation
    handles.spectreData = zebraguidata(handles.spectre);
    % zebraguidata(hObject, handles);
    % we are here after opening a new data. Set the proper default values
    % in the ZebraSpectre GUI
    set(handles.spectreData.fromTtxt,'string','0');      % beginning of the time window for spectre computation
    set(handles.spectreData.toTtxt,'string',handles.tmax);      
    % end of the time window for spectre computation
    
    % ZebraSpectre('computeSpectra_Callback',handles.spectreData.computeSpectra,[],handles.spectreData)
    %computePowerSpectra (hObject, handles);
    % handles = zebraguidata(hObject);    % this to refresh the local values of handles 
end

if get(handles.analyzePWdistributionCk, 'Value')
    analyzePowerDistribution (hObject, handles); 
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
end

if (handles.nTrials > 1)
    disp([num2str(toc) ' Compute averages on trials.'])
    % compute al relevant averages. First, mean LFP
    % compute mean band passed data
    % compute mean spectrogram and 'gamma' band power 
    computeMeanOnTrials (hObject,handles)    
end    
handles = zebraguidata(hObject);        % this to refresh the local values of handles

% Insert a breakpoint here if you want to access directly to the averaged
% trials. handles.meanLFP (nCh, 1:handles.dtaLen)


disp([num2str(toc) ' Plot data.'])
plotWideWindow (hObject, handles);
plotMultiBand(hObject, handles);
plotSpectrogram (hObject, handles);

handles = zebraguidata(hObject);    % this to refresh the local values of handles

disp([num2str(toc) ' Create summary table.'])
% SUMMARY TABLE
% compute the total RMS power in the five bands
if (~flagSSP)
    handles.lpRMS = 0;
    handles.bp1RMS = 0;
    handles.bp2RMS = 0;
    handles.bp3RMS = 0;
    handles.hpRMS = 0;
    for i = 1:handles.nCh
        if ~(flagExSynch && i==sCh) % if ext synch the channel is excluded from the output
            if ~(flagExSynch && flagEx && i==exCh) % exclude channel
                for j = 1:handles.nTrials
                    handles.cnt = handles.cnt + 1;
                    if handles.BPcomputed(1)
                        tmp(1:handles.dtaLen) = handles.bandPassed_LFP(1,i,j,1:handles.dtaLen); 
                        handles.lpRMS = rms(tmp);
                    end
                    if handles.BPcomputed(2)
                        tmp(1:handles.dtaLen) = handles.bandPassed_LFP(2,i,j,1:handles.dtaLen);
                        handles.bp1RMS = rms(tmp);
                    end                    
                    if handles.BPcomputed(3)
                        tmp(1:handles.dtaLen) = handles.bandPassed_LFP(3,i,j,1:handles.dtaLen);
                        handles.bp2RMS = rms(tmp);
                    end
                    if handles.BPcomputed(4)
                        tmp(1:handles.dtaLen) = handles.bandPassed_LFP(4,i,j,1:handles.dtaLen);
                        handles.bp3RMS = rms(tmp);
                    end
                    if handles.BPcomputed(5)
                        tmp(1:handles.dtaLen) = handles.bandPassed_LFP(5,i,j,1:handles.dtaLen);
                        handles.hpRMS = rms(tmp);
                    end
                    % now we place this data in a summary table that will be loaded with a
                    % number of fun fact about the data...
                    % data must be converted to the cell format in the following statements

                    handles.dtaSummary(handles.cnt,1) = cellstr(handles.dir_in);
                    handles.dtaSummary(handles.cnt,2) = cellstr(handles.file_in);
                    if handles.autoMode
                        handles.dtaSummary(handles.cnt,3) = cellstr('0');
                    else
                        handles.dtaSummary(handles.cnt,3) = cellstr(handles.fileTime);
                    end
                    handles.dtaSummary(handles.cnt,4) = num2cell(i);
                    handles.dtaSummary(handles.cnt,5) = num2cell(j);
                    handles.dtaSummary(handles.cnt,6) = num2cell(handles.lpRMS);
                    handles.dtaSummary(handles.cnt,7) = num2cell(handles.bp1RMS);
                    handles.dtaSummary(handles.cnt,8) = num2cell(handles.bp2RMS);
                    handles.dtaSummary(handles.cnt,9) = num2cell(handles.bp3RMS);
                    handles.dtaSummary(handles.cnt,10) = num2cell(handles.hpRMS);
                    handles.dtaSummary(handles.cnt,11) = cellstr('Comments here');
                end
            end
        end
    end 
    set(handles.zebraSaysData.outable,'Data',handles.dtaSummary);
else
    % Put here something if you want to compute 'some' metric of the SSP
    % power...
end
zebraguidata(hObject, handles);
%%

function plotWideWindow (hObject, handles)
%% The plot function requires a column vector. If the data is a matrix, the X axis moves
% along the row. If more than n>1 columns are present, the function will
% plot the various data in different colours.

% Find out the index of the left and right limits of the plot window
tMin = str2double(get(handles.tminTxt,'String'));
tMax = str2double(get(handles.tmaxTxt,'String'));
% convert the limits to array indexes
tMinI = floor(tMin/handles.sp+1);
tMaxI = floor(tMax/handles.sp+1);
% quality control
if tMinI <= 1, tMinI = 1;
end

if tMaxI > handles.dtaLen, tMaxI = handles.dtaLen;
end

handles.minLFP = min(handles.LFP(handles.currentCh, handles.currentTrial,tMinI:tMaxI));
handles.maxLFP = max(handles.LFP(handles.currentCh, handles.currentTrial,tMinI:tMaxI));
traceOut = [];

axes(handles.wide_plot_up);
cla
gca.ColorOrderIndex = 0;
hold on;

flagMd = get(handles.plotMeanFlag,'Value');     % plot mean trace
flagMP = get(handles.plotMultiCh,'Value');      % plot all channels
flagTr = get(handles.plotTrCk,'Value');         % plot single trials
flagAutoOff = get(handles.autoOffCk,'Value');
flagTempl = get(handles.plotTemplateCk,'Value');

if (flagAutoOff==0), range = 0;  % take this out if you want to auto-offset the averages
else
    range = (handles.maxLFP-handles.minLFP)/3;
end

if flagMP       % plot all channels
    % if multi channel no template!
    % if flagTr    Not sure what this flag stand for...
    traceOut(1:handles.nCh,1:handles.dtaLen) = handles.LFP(1:handles.nCh,handles.currentTrial,1:handles.dtaLen);
    for (i=1:handles.nCh)       % add a vertical offset
        traceOut(i,1:handles.dtaLen) = (i-1)*range + traceOut(i,1:handles.dtaLen);
    end
    % next: plot multiplot with the same colors
    gca.ColorOrderIndex = 1;
    gca.nextPlot = 'add';
    if flagTr   % plot the individual trials
        % the additional properties are necessary to stop the plot to
        % intercepts the mouse clicks
        plot(handles.x(tMinI:tMaxI),traceOut(:, tMinI:tMaxI), 'HitTest', 'off');
    end
    if (flagMd && handles.nTrials > 1)  % plot all means with offset
        traceOut(1:handles.nCh,1:handles.dtaLen) = handles.meanLFP(1:handles.nCh,1:handles.dtaLen);
        for (i=1:handles.nCh)
            traceOut(i,1:handles.dtaLen) = (i-1)*range + traceOut(i,1:handles.dtaLen);
        end
        gca.ColorOrderIndex = 1;
        plot(handles.x(tMinI:tMaxI),traceOut(:, tMinI:tMaxI),'LineWidth',2, 'HitTest', 'off');
    end
    
    % next: plot multiplot with the same colors
    %axx = gca;
else    
    if flagTr   % plot the individual trials
        traceOut(1:handles.dtaLen) = handles.LFP(handles.currentCh,handles.currentTrial,1:handles.dtaLen);
        plot(handles.x(tMinI:tMaxI),traceOut(tMinI:tMaxI),'Color','k', 'HitTest', 'off');
    end    
    if (flagMd && handles.nTrials > 1)
        plot(handles.x(tMinI:tMaxI),handles.meanLFP(handles.currentCh,(tMinI:tMaxI)),'Color','k',...
            'LineWidth', 2, 'HitTest', 'off');
    end
    if (flagTempl && handles.nTrials>1 && flagTr)
        % plot the template for each response with the correct scaling
        % factor and offset
        l = length(handles.templ);
        t0 = str2double(get(handles.templateFrom,'String'));
        for ii= 1:handles.ntemp
            a = handles.EPamplitude(handles.currentCh,handles.currentTrial,ii); 
            templ = a * handles.templ + handles.EPoffset(handles.currentCh,handles.currentTrial,ii);
            x0 = t0 + handles.dtaLen*handles.sp/2 * (ii-1);
            xend = x0 + l*handles.sp;
            xtempl = x0:handles.sp:xend;
            plot (xtempl(1:l), templ(1:l),'LineWidth',2, 'HitTest', 'off');
        end    
    end    
end
flagUnits = get (handles.microVck,'Value');
if flagUnits, ylabel('LFP (microV)');
else
    ylabel('LFP (mV)');
end    
axis ([tMin tMax -inf inf]);
%hold all;
%%hold off
zebraguidata(hObject, handles);

function wide_plot_up_ButtonDownFcn(hObject, eventdata, handles)
%% 2017, 12 27. 
% Executes on mouse press over the graph of the raw data. The code has two
% modes of operation.
% default: the mouse click select the center of the data segment to be maginified in the insets.
% reas screen:

persistent marker
persistent oldX oldY

if (ishandle(marker)), delete(marker)   % delete any previously drawn marker
end

flag = 1;
if flag
    % perform the callback in the screen data mode. Use the CurretPoint
    % property of the current axes to obtain the x, y coordinates.
    cp = get (gca, 'CurrentPoint'); %GAB: this should be avoided: if the current figure has no axes,
                                    % gca creates a new axes in that figure!
    x = cp(1,1);
    y = cp(1,2);
    mouseBtn = get(gcf,'selectiontype');    % read the button event. The string can assume the following values:
    % 'normal': left button, single click
    % 'extend': Shift-click left button
    % 'alt': Cntrl-click left or click right button
    % 'open' Double click any button
    switch mouseBtn
    case 'normal'
        % the axes coordinates are saved in the persistent variables
        % oldX, oldY. This call also clear the clipboard
        oldX = x;
        oldY = y;
        oldBtn = mouseBtn;
        clipboard('copy',' ');   % clear the clipboard content
    case 'extend'
        % if the previous event mouse button press was 'normal' deltaY
        % is computed as the difference between the new y and the old
        % y. This is used to compute manually the response amplitude
        % when th baseline has been identified with the most recent left mouse 
        % click. The resulting x and deltaY are saved to the clipboard.
        deltaY = abs(y - oldY);
        % read the clipboard content as a string
        dtaOut = clipboard('paste');
        % now append the latest data as a tab separated list. This
        % feature allows to measure multiple peaks in the waveform.
        if strcmp(dtaOut, ' ')
            dtaOut = sprintf('%f \t %f',x,deltaY);   % start a new clipboard string
        else    
            dtaOut = [dtaOut sprintf('\t %f \t %f',x,deltaY)];    % append to the clipboard
        end        
        clipboard('copy',dtaOut);   % copy the coordinates to the clipboard
    end        
else
    [x, y] = ginput(1);     % ginput returns the coordinates of the cliked point
    % now check if x is within the data limit and if the clicked field is
    % within the raw data window. hObject is the handle to the wide plot window (since 
    % it generated the ButtonDownFcn call) and gca is the handle of the current graphic object. 
    % For the program to function correctly the two handles must be identical.
    if (x >= 0) && (x <= handles.dtaLen*handles.sp) && (hObject == gca)
        handles.xnow = x;
        handles.ynow = y;
        x_from = handles.xnow - handles.insetW/2;
        if (x_from <= 0), x_from = 0;
            handles.xnow = handles.insetW/2;
        end    
        x_to = handles.xnow + handles.insetW/2;       % this is in time units, convert now in indexes of the arrays
        if (x_to > handles.dtaLen*handles.sp), x_to = handles.dtaLen*handles.sp;
            handles.xnow = x_to - handles.insetW/2;
            x_from = handles.xnow - handles.insetW/2;
        end
        % draw a rectangle marker around the magnified area
        hmarker = handles.maxLFP-handles.minLFP;
        if (hmarker <= 0), hmarker = 1; 
        end
        marker = rectangle('Position', [x_from .minLFP handles.insetW hmarker],'EdgeColor','red'); 
        %set(marker,'Color','red');
        %set(marker,'FaceAlpha',0.5);

        zebraguidata(hObject, handles);
        doInsets (hObject, handles);
    end
end
%%

function data_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function insetWnd_Callback(hObject, eventdata, handles)
handles.insetW = str2double(get(hObject,'String')); 
zebraguidata(hObject, handles);
doInsets (hObject, handles);

function insetWnd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spEdit_Callback(hObject, eventdata, handles)
% hObject    handle to spEdit (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

handles.sp = str2double(get(hObject,'String')); 
handles.acqF = 1/handles.sp;
handles.frequencyBand(5,2) = 0.5*handles.acqF;  %GAB 22/05/2018
zebraguidata(hObject, handles);

function spEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function chEdit_Callback(hObject, eventdata, handles)
handles.nCh = str2double(get(hObject,'String'));

if handles.nCh>1
    set (handles.chSlider, 'SliderStep', [1/(handles.nCh-1) 1/(handles.nCh-1)]);
    set (handles.ChText,'String',['Ch ' num2str(handles.currentCh)]);
    set (handles.chSlider, 'Max', handles.nCh);
    set (handles.chSlider, 'Value', 1);
else
    set (handles.chSlider,'Visible','off');
    set (handles.ChText,'Visible','off');
end

handles.currentCh = 1;
zebraguidata(hObject, handles);

function chEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gain_Callback(hObject, eventdata, handles)
handles.gainAmp = str2double(get(hObject,'String'));
zebraguidata(hObject, handles)

function gain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lpcheck_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    plotSpectrogram (hObject, handles)
    doInsets (hObject, handles);
end

function bpCheck1_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    plotSpectrogram (hObject, handles)
    doInsets (hObject, handles);
end

function bpCheck2_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    plotSpectrogram (hObject, handles)
    doInsets (hObject, handles);
end

function bpCheck3_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    plotSpectrogram (hObject, handles)
    doInsets (hObject, handles);
end

function hpCheck_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    plotSpectrogram (hObject, handles)
    doInsets (hObject, handles);
end

function alphaFrom_Callback(hObject, eventdata, handles)
handles.frequencyBand(2,1) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function alphaFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alphaTo_Callback(hObject, eventdata, handles)
handles.frequencyBand(2,2) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function alphaTo_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gammaFrom_Callback(hObject, eventdata, handles)
handles.frequencyBand(4,1) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function gammaFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gammaTo_Callback(hObject, eventdata, handles)
handles.frequencyBand(4,2) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function gammaTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function betaFrom_Callback(hObject, eventdata, handles)
handles.frequencyBand(3,1) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function betaFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lpTo_Callback(hObject, eventdata, handles)
handles.frequencyBand(1,2) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function lpTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hpFrom_Callback(hObject, eventdata, handles)
handles.frequencyBand(5,1) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function hpFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function betaTo_Callback(hObject, eventdata, handles)
handles.frequencyBand(3,2) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function betaTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotMultiBand(hObject, handles)
% Here we plot the multi band window according to the selected bands
flagSSP = get(handles.meyerSspRB,'Value');

if (~flagSSP)
    tmp = [];
    axes(handles.wide_plot_HP);
    cla
    hold all;
    flagTr = get(handles.plotTrCk,'Value');
    flagMd = get(handles.plotMeanFlag,'Value');
    flag = get(handles.lpcheck,'Value');
    if get(handles.lpcheck,'Value') && handles.BPcomputed(1)
        if flagTr
            tmp(1:handles.dtaLen) = handles.bandPassed_LFP(1,handles.currentCh, handles.currentTrial, 1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','k','LineWidth',0.5);
        end
        if flagMd
            tmp(1:handles.dtaLen) = handles.meanBP(1, handles.currentCh,1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','k','LineWidth',2);
        end      
    end
    
    if get(handles.bpCheck1,'Value') && handles.BPcomputed(2)
        if flagTr
            tmp(1:handles.dtaLen) = handles.bandPassed_LFP(2,handles.currentCh, handles.currentTrial, 1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','c','LineWidth',0.5);
        end
        if flagMd
            tmp(1:handles.dtaLen) = handles.meanBP(2, handles.currentCh,1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','c','LineWidth',2);
        end      
    end
    
    if get(handles.bpCheck2,'Value') && handles.BPcomputed(3)
        if flagTr  
            tmp(1:handles.dtaLen) = handles.bandPassed_LFP(3,handles.currentCh, handles.currentTrial, 1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','r','LineWidth',0.5);
        end
        if flagMd
            tmp(1:handles.dtaLen) = handles.meanBP(3, handles.currentCh,1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','r','LineWidth',2);
        end      
    end
    
    if get(handles.bpCheck3,'Value') && handles.BPcomputed(4)
        if flagTr
            tmp(1:handles.dtaLen) = handles.bandPassed_LFP(4,handles.currentCh, handles.currentTrial, 1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','g','LineWidth',0.5);
        end
        if flagMd
            tmp(1:handles.dtaLen) = handles.meanBP(4, handles.currentCh,1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','g','LineWidth',2);
        end      
    end
    
    if get(handles.hpCheck,'Value') && handles.BPcomputed(5)
        if flagTr    
            tmp(1:handles.dtaLen) = handles.bandPassed_LFP(5,handles.currentCh, handles.currentTrial, 1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','blue');
        end
        if flagMd
            tmp(1:handles.dtaLen) = handles.meanBP(5, handles.currentCh,1:handles.dtaLen);
            plot(handles.x(1:handles.dtaLen),tmp(1:handles.dtaLen),'Color','b','LineWidth',2);
        end      
    end
    axis ([handles.tmin handles.tmax -inf inf])
    set(gca, 'XTick', [])
    flagUnits = get (handles.microVck,'Value');
    if flagUnits, ylabel('LFP (microV)');
    else
        ylabel('LFP (mV)');
    end    
end

function spectraInset_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to spectraInset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x y] = ginput(1)
xl =[x x]
yl = [y-2 y+2]
h = line(xl,yl)

function inset_HP_ButtonDownFcn(hObject, eventdata, handles)

persistent csrFlag cursor1 cursor2 time1 time2
wsOut = [0 0 0 0; 1 1 1 1; 2 2 2 2];

[x y] = ginput(1);
% First, check if we are within the limits of the axes window
edges = xlim;
if (x >= edges(1)) && (x <= edges(2)) && (gca == hObject)
    if csrFlag >= 1
    else
       csrFlag = 0; % initialize csrFlag value
    end
    ylim('manual');

    if csrFlag == 3
        % this is not the first time that we measure an interval
        csrFlag = 1;
        % delete lines
        %if isvalid(cursor1,'handle')     %check if the cursor handle exists and is valid
        if ishandle(cursor1)
            delete (cursor1);            % if valid, there are cursors already present that 
            delete (cursor2);            % have to be removed.
        end
    end    
    ymin = min(handles.bandPassed_LFP(2,:));
    ymax = max(handles.bandPassed_LFP(2,:));

    if (x>0)        % make sure that you are not out of the window
        if (csrFlag <= 1)
            csrFlag = 2;
            cursor1 = 0;
            cursor2 = 0;
            cursor1 = line([x x],[ymin ymax]);
            time1 = x;
        else
            handles.eventNum = handles.eventNum + 1;
            csrFlag = 3;
            cursor2 = line([x x],[ymin ymax]);
            time2 = x;
            % compute the power in the three bands and in the baseline
            % first: load the sample in a subvector that contains only the required
            % segment. Second: compute the rms power. Third: compile the output table.
            timeCenter = (time1+time2)/2 ;
            deltaTime = abs(time2-time1);
            i1 = round(time1/handles.sp)+1;
            i2 = round(time2/handles.sp)+1;
            betaSegment = handles.bandPassed_LFP(2,handles.currentCh, handles.currentTrial,i1:i2);
            gammaSegment = handles.bandPassed_LFP(3,handles.currentCh, handles.currentTrial,i1:i2);
            hpSegment = handles.bandPassed_LFP(4,handles.currentCh, handles.currentTrial,i1:i2);
            baseT = str2double(get(handles.baselinePower,'string'));
            ib1=round((time1-baseT)/handles.sp)+1;
            ib2=round(time1/handles.sp);
            betaBase = handles.bandPassed_LFP(2,handles.currentCh, handles.currentTrial,ib1:ib2);
            gammaBase = handles.bandPassed_LFP(3,handles.currentCh, handles.currentTrial,ib1:ib2);
            hpBase = handles.bandPassed_LFP(4,handles.currentCh, handles.currentTrial,ib1:ib2);
            wsOut(:,4) = deltaTime;
            wsOut(1,2) = rms(betaSegment);
            wsOut(2,2) = rms(gammaSegment);
            wsOut(3,2) = rms(hpSegment);
            wsOut(1,1) = rms(betaBase);
            wsOut(2,1) = rms(gammaBase);
            wsOut(3,1) = rms(hpBase);
            wsOut(1,3) = wsOut(1,2)/wsOut(1,1);
            wsOut(2,3) = wsOut(2,2)/wsOut(2,1);
            wsOut(3,3) = wsOut(3,2)/wsOut(3,1);
            Eratio = (wsOut(1,2)-wsOut(1,1))/(wsOut(2,2)-wsOut(2,1)); 
            set(handles.outTable,'Data',wsOut);
            % now fill the summary table in ZebraSays.fig
            handles.eventSummary(handles.eventNum,1) = cellstr(handles.dir_in);
            handles.eventSummary(handles.eventNum,2) = cellstr(handles.file_in);
            handles.eventSummary(handles.eventNum,3) = num2cell(handles.eventNum);
            handles.eventSummary(handles.eventNum,4) = num2cell(deltaTime);
            handles.eventSummary(handles.eventNum,5) = num2cell(timeCenter);
            handles.eventSummary(handles.eventNum,6) = num2cell(wsOut(1,1));
            handles.eventSummary(handles.eventNum,7) = num2cell(wsOut(1,2));
            handles.eventSummary(handles.eventNum,8) = num2cell(wsOut(1,3));
            handles.eventSummary(handles.eventNum,9) = num2cell(wsOut(2,1));
            handles.eventSummary(handles.eventNum,10) = num2cell(wsOut(2,2));
            handles.eventSummary(handles.eventNum,11) = num2cell(wsOut(2,3));
            handles.eventSummary(handles.eventNum,13) = num2cell(wsOut(3,1));
            handles.eventSummary(handles.eventNum,14) = num2cell(wsOut(3,2));
            handles.eventSummary(handles.eventNum,15) = num2cell(wsOut(3,3));
            handles.eventSummary(handles.eventNum,12) = num2cell(Eratio);
            handles.eventSummary(handles.eventNum,16) = cellstr('Notes');    
            % push the data in the table
            set(handles.zebraSaysData.outEventsSummary,'Data',handles.eventSummary);
        end
        ylim('auto');
    end
end
zebraguidata(hObject, handles);

function baselinePower_Callback(hObject, eventdata, handles)

function baselinePower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function viewOutData_Callback(hObject, eventdata, handles)
v = get(handles.OutTable,'visible');
if strcmp(v,'on'), set(handles.OutTable,'visible','off')
else
    set(handles.OutTable,'visible','on')
end

function notchFilter_Callback(hObject, eventdata, handles)
computeBP_Callback(hObject, eventdata, handles);
handles = zebraguidata(hObject);
computePowerSpectra (hObject, handles);


function notchFrom_Callback(hObject, eventdata, handles)
handles.notch(1) = str2double(get(hObject,'String')); 
zebraguidata(hObject, handles);

function notchFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function notchTo_Callback(hObject, eventdata, handles)
handles.notch(2) = str2double(get(hObject,'String')); 
zebraguidata(hObject, handles);


function notchTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function selectWnd_CreateFcn(hObject, eventdata, handles)

function selectWnd_ButtonDownFcn(hObject, eventdata, handles)

function selectWnd_SelectionChangeFcn(hObject, eventdata, handles)

%p=get(handles.dataProcessing,'Parent')
%t=get(p,'Tag')

y = get(hObject,'Tag');     % obtain the tag to the selected button
clearMainWindow (hObject, handles, y);

function clearMainWindow (hObject, handles,y)

if strcmp(y, 'dataProcessingCk')
    scButtonOff(hObject, handles)
    set(get(handles.wide_plot_up,'children'),'visible','on') %hide the axes contents 
    set(handles.wide_plot_up,'visible','on') %hide the axes contents 
    set(handles.autoProcessing,'visible','off');
    set(handles.wide_plot_HP,'visible','off');
    set(get(handles.wide_plot_HP,'children'),'visible','off') %hide the axes contents 
    set(handles.spectraWide,'visible','off');
    set(get(handles.spectraWide,'children'),'visible','off') %hide the current axes contents 
    set(handles.inset_HP,'visible','off');
    set(get(handles.inset_HP,'children'),'visible','off') %hide the axes contents 
    set(handles.spectraInset,'visible','off');
    set(get(handles.spectraInset,'children'),'visible','off') %hide the axes contents 
    
    set(handles.dataProcessing,'visible','on');
end

if strcmp(y, 'autoprocessingCk')
    scButtonOn(hObject, handles)
    set(handles.dataProcessing,'visible','off')
    set(get(handles.wide_plot_up,'children'),'visible','off') %hide the axes contents 
    set(handles.wide_plot_up,'visible','off') %hide the axes contents 
    set(handles.wide_plot_HP,'visible','on');
    set(handles.spectraWide,'visible','on');
    set(get(handles.spectraWide,'children'),'visible','on') %hide the current axes contents 
    set(handles.inset_HP,'visible','on');
    set(get(handles.inset_HP,'children'),'visible','on') %hide the axes contents 
    set(handles.spectraInset,'visible','on');
    set(get(handles.spectraInset,'children'),'visible','on') %hide the axes contents 
    set(handles.autoProcessing,'visible','on')
end

if strcmp(y, 'mainVisualizer')
    scButtonOn(hObject, handles)
    set(get(handles.wide_plot_up,'children'),'visible','on') %hide the axes contents 
    set(handles.wide_plot_up,'visible','on') %hide the axes contents 
    set(handles.dataProcessing,'visible','off')
    set(handles.autoProcessing,'visible','off')
    set(handles.wide_plot_HP,'visible','on');
    set(get(handles.wide_plot_HP,'children'),'visible','on');
    set(handles.spectraWide,'visible','on');
    set(handles.spectraWide,'visible','on');
    set(get(handles.spectraWide,'children'),'visible','on') 
    set(handles.inset_HP,'visible','on');
    set(get(handles.inset_HP,'children'),'visible','on')  
    set(handles.spectraInset,'visible','on');
    set(get(handles.spectraInset,'children'),'visible','on') 
end    

function predefinedBW_SelectionChangeFcn(hObject, eventdata, handles)
% select default values for bandpass filters
y = get(hObject,'Tag');     % obtain the tag to the selected button
leak = str2double(get(handles.leakageLP,'String'));

if strcmp(y, 'zebraCk')
    handles.frequencyBand = [0 8; 0.4 4; 4 8; 8 12; 100 0.5/handles.sp; 0 leak];
    % now update the edit fields
end

if strcmp(y, 'EEGck')
    handles.frequencyBand = [0 4; 0.5 4; 6 12; 12 30; 100 0.5/handles.sp; 0 leak];
end

if strcmp(y, 'mouseCk')
    handles.frequencyBand = [0 8; 0.5 4;6 12; 30 80; 100 0.5/handles.sp; 0 leak];
end

if strcmp(y, 'NewcastleBand')
    handles.frequencyBand = [0 8; 0.5 4; 20 60; 70 300; 300 0.5/handles.sp; 0 leak]; 
end

% update the edit fields
set (handles.lpTo,'String',handles.frequencyBand (1,2));
set (handles.alphaFrom,'String',handles.frequencyBand (2,1));
set (handles.alphaTo,'String',handles.frequencyBand (2,2));
set (handles.betaFrom,'String',handles.frequencyBand (3,1));
set (handles.betaTo,'String',handles.frequencyBand (3,2));
set (handles.gammaFrom,'String',handles.frequencyBand (4,1));
set (handles.gammaTo,'String',handles.frequencyBand (4,2));
set (handles.hpFrom,'String',handles.frequencyBand (5,1));

zebraguidata(hObject, handles);

function computeBP_Callback(hObject, eventdata, handles)
%% This function computes and plots the filtered data. It must be evoked after
% changing any filter settings to update the filtered data.
% May 7, 2017. 
% I have added the temporal profile of the power in BP1, BP2, BP3 and HP

if (handles.fileInNum) % of course the following code is executed only if a data has been loaded
    computeBPandNotch (hObject, handles);   % handles will be modified at this call
    handles = zebraguidata(hObject);             % this to refresh the local values of handles
    computePowerInTime (hObject, handles);  % temporal profile of power
    handles = zebraguidata(hObject);             % this to refresh the local values of handles
    plotMultiBand(hObject, handles);
    computeSpectrogram (hObject, handles)
    handles = zebraguidata(hObject);             % this to refresh the local values of handles
    plotSpectrogram (hObject, handles)
    doInsets (hObject, handles);
end
%%


function computeBPandNotch(hObject, handles)
%%
% March 12, 2018
% code has been modified to able/enable the BP computations

% first: check if the notch filter is selected.
if (get(handles.notchFilter,'Value'))
    [b,a] = butter(3,2*handles.notch*handles.sp,'stop');    %handles.notch is an array with the two param
    for i = 1:handles.nCh
        for j= 1:handles.nTrials
            tmp (1:handles.dtaLen) = handles.LFP (i,j,1:handles.dtaLen);    %.'
            tmp2 = filtfilt(b,a,tmp.');
            %handles.workLFP (i,j,1:handles.dtaLen) = filtfilt(b,a,tmp.');       
            %handles.workLFP (i,j,1:end) = filtfilt(b,a,handles.LFP.');
            handles.workLFP(i,j,1:handles.dtaLen) = tmp2(1:handles.dtaLen);
        end
    end    
else
    handles.workLFP = handles.LFP;    
end
% computations are performed starting from the workLFP file

handles.filterOrder = str2double(get(handles.filtOrder,'String'));
% compute band passed data
if get(handles.LPenableCK,'Value')
    [b,a] = butter(handles.filterOrder,2*handles.frequencyBand(1,2)*handles.sp,'low');  %Low pass filter
    computeFilteredData (hObject, handles, b, a, 1);
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
    handles.BPcomputed(1) = 1;  %GAB 06/04/2018 this line has to stay after the zebraguidata call
end

if get(handles.BP1enableCK,'Value')
    [b,a] = butter(handles.filterOrder,2*handles.frequencyBand(2,:)*handles.sp);        % band pass 1
    computeFilteredData (hObject, handles, b, a, 2)
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
    handles.BPcomputed(2) = 1;
end

if get(handles.BP2enableCK,'Value')
    [b,a] = butter(handles.filterOrder,2*handles.frequencyBand(3,:)*handles.sp);        % band pass 2
    computeFilteredData (hObject, handles, b, a, 3)
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
    handles.BPcomputed(3) = 1;
end

if get(handles.BP3enableCK,'Value')
    [b,a] = butter(handles.filterOrder,2*handles.frequencyBand(4,:)*handles.sp);        % band pass 3
    computeFilteredData (hObject, handles, b, a, 4)
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
    handles.BPcomputed(4) = 1;
end

if get(handles.HPenableCK,'Value')
    [b,a] = butter(handles.filterOrder,2*handles.frequencyBand(5,1)*handles.sp, 'high');% high pass and single units
    computeFilteredData (hObject, handles, b, a, 5)
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
    handles.BPcomputed(5) = 1;
end
% compute the LP data for spectral leakage correction
if get(handles.leakageCk,'Value')
    [b,a] = butter(handles.filterOrder,2*handles.frequencyBand(6,2)*handles.sp,'low');  %Low pass filter
    computeFilteredData (hObject, handles, b, a, 6)
    handles = zebraguidata(hObject);    % this to refresh the local values of handles
    handles.BPcomputed(6) = 1;
end

% envelope computation in beta+gamma band for the automatic ID of events
%tmp = abs(handles.bandPassed_LFP(2,:)+handles.bandPassed_LFP(3,:));
%handles.envelope = abs(hilbert(tmp));
% computer rolling average 
%handles.rollingAV = sgolayfilt(tmp,3,1001);

%zebraguidata(gcbo, handles);
zebraguidata(hObject,handles);
%%

function computeFilteredData (hObject, handles, b, a, k)
%% b, a are the descriptor of the filter. k is the index of filtered data
% 1 LP; 2-4 three BP; 5 HP; 6 low pass for spectral leakage
for i = 1:handles.nCh
    tmp1 = [];
    tmp2 = [];
    for j= 1:handles.nTrials
        tmp1(1:handles.dtaLen) = handles.workLFP(i,j,1:handles.dtaLen);    %.'
        tmp2 = filtfilt(b,a,tmp1.');
        handles.bandPassed_LFP(k,i,j,:)=tmp2(:);
    end
end    
%zebraguidata(gcbo, handles);
zebraguidata(hObject,handles);
%%
function plotEnvelopeCk_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    doInsets (hObject, handles);
end

function plotBPpowerCk_Callback(hObject, eventdata, handles)
if (handles.fileInNum), plotMultiBand(hObject, handles)
    doInsets (hObject, handles);
    plotSpectrogram (hObject, handles)
end

function computePowerInTime (hObject, handles)
%% May 8, 2017
% Here we compute the short time RMS power for the four band passed data
% The computation is performed on a window which amplitude is given by the 
% larger of 1) the time window of the SPG computation (ChronuxRW) or
% 2) the lower period of the bandwidth. 
% The overlap is one fifth of this.
% Therefore each STRMS plot has its own sampling, so each array contains
% both the x and the power

chronuxWin = str2double(get (handles.ChronuxRW,'string'));

% First, preallocation of the array handles.BPpower
% This array has different lengths for the different bands. The max
% length is given by the window size determined by the value of
% ChronuxWin. 
pntWin = floor((chronuxWin/handles.sp)/2);        % Cunning strategy! See later!  
pntOvl = floor(pntWin/5);
pntWin = 2* pntOvl * 5;                           % In this way pntWin is a multiple of 10
centerBegin = pntWin/2+1;
BPindexes = centerBegin:pntOvl:(handles.dtaLen-pntWin/2);  % this array is defined to preallocate it and to compute 
maxDim = length(BPindexes);                               % the size of the longest array

% Each BP power data has its own length that is stored in the array
% handles.BPpowerL
handles.BPpowerL = zeros(5);
% indexes of BPpower: 1) 1-> x, 2-> power; 2) band pass; 3) channel; 4)
% trial; 6) time.
handles.BPpower = zeros(2, 5, handles.nCh, handles.nTrials, maxDim);

% Lets cycle on the filtered data
for iCh = 1:handles.nCh
    for jTr= 1:handles.nTrials
        for BPsel=2:5
            if handles.BPcomputed(BPsel)   % perform the computation only if that BP has been computed
                % first compute the x
                BPindexes = [1];    % initialize the index array
                BPindLeft = [1];
                BPindRight = [1];
                % compute the width of the window for power computation
                lp = 1/handles.frequencyBand(BPsel,1);      % lower limit of the bandpass
                winT = max([lp chronuxWin]);                % for higher frequency I use the window used for SPG computation
                pntWin = floor((winT/handles.sp)/2);        % Cunning strategy! See later!  
                pntOvl = floor(pntWin/5);
                pntWin = 2* pntOvl * 5;                     % In this way pntWin is a multiple of 10
                pntWin2 = pntWin/2;
                centerBegin = pntWin/2+1;
                %pntOvl = pntWin;  %GAB!! PROVA MOMENTANEA!!!!
                BPindexes = centerBegin:pntOvl:(handles.dtaLen-pntWin2);  % this array contains the indexes of the 
                                                                          % BP data around which are computed the RMS powers
                xBPpowerEnd = length(BPindexes);
                handles.BPpowerL(BPsel) = xBPpowerEnd;
                if xBPpowerEnd<1, xBPpowerEnd=1;
                    % this condition can only occur when the data is too short
                    % for the specified bandwidth. This is done to prevent an
                    % ensuing error condition, but, of course, no power can be
                    % computed in this bandwidth.
                    BPindexes(1) = handles.dtaLen/2;
                    % the power plot contains a single point at the middle of
                    % the trial.
                end
                BPindLeft(1:xBPpowerEnd) = BPindexes(1:xBPpowerEnd)-pntWin2;    % left and right limits (in sampling points)
                BPindRight(1:xBPpowerEnd) = BPindexes(1:xBPpowerEnd)+pntWin2;   % of the integration window
                if (BPindLeft(1)<1), BPindLeft(1) = 1;
                end
                if (BPindRight(xBPpowerEnd)>handles.dtaLen), BPindRight(xBPpowerEnd)=handles.dtaLen;
                end
                % compute X axis
                handles.BPpower (1,BPsel,iCh,jTr,1:xBPpowerEnd) = handles.sp * (BPindexes(1:xBPpowerEnd)-1);  % time coordinate
                if (xBPpowerEnd>1)         % again, more controls to protect bad bandwidths
                    handles.BPpowerSP(BPsel) = handles.sp * (BPindexes(2)-BPindexes(1));
                else    
                    handles.BPpowerSP(BPsel) = 0;
                end    
                for timI = 1:xBPpowerEnd
                    % this loop computes the power
    %                 handles.BPpower (2,BPsel,iCh,jTr,1:xBPpowerEnd) = rms...
    %                     (handles.bandPassed_LFP(BPsel,iCh,jTr,BPindLeft(1:xBPpowerEnd):BPindRight(1:xBPpowerEnd)),4);
                    handles.BPpower (2,BPsel,iCh,jTr,timI) = rms...
                        (handles.bandPassed_LFP(BPsel,iCh,jTr,BPindLeft(timI):BPindRight(timI)),4);
                end
            end
        end
    end
end
zebraguidata(hObject, handles);


function spectraWide_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to spectraWide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x, y] = ginput(1)     % ginput returns the coordinates of the clicked point

% now check if x is within the data limit and if the clicked field is
% within the raw data window. Boy this is great!
if (x >= 0) && (x <= handles.dtaLen*handles.sp) && (hObject == gca)
    
end    


function wide_plot_HP_ButtonDownFcn(hObject, eventdata, handles)
% --- Executes on mouse press over axes background.
% hObject    handle to wide_plot_HP (see GCBO)

[x y] = ginput(1);     % ginput returns the coordinates of the cliked point

% now check if x is within the data limit and if the clicked field is
% within the raw data window. Boy this is great!
if (x >= 0) && (x <= handles.dtaLen*handles.sp) && (hObject == gca)
    x_from = x - 0.25;
    if (x_from <= 0), x_from = 0;
        x = 0.25;
    end    
    x_to = x + 0.25;       % this is in time units (sec)
    if (x_to > handles.dtaLen*handles.sp), x_to = handles.dtaLen*handles.sp;
        x = x_to - .25;
        x_from = x - .25;
    end
    % ok! The interval is well defined. Extract the data segment and
    % Fourierize it.
    i_from = round(x_from/handles.sp) + 1;
    i_to = round(x_to/handles.sp) + 1;
    segment = handles.workLFP(i_from:1:i_to);
    [pxx,f] = pwelch(segment,[],[],[],handles.acqF);
    axes(handles.spectreData.powerSpectra);
    hold all
    plot(f,10*log10(pxx));
    axis ([1 100 -70 -30]);
end   

function trialsEdit_Callback(hObject, eventdata, handles)

function trialsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function chSlider_Callback(hObject, eventdata, handles)
% --- Executes on slider movement.
% hObject    handle to chSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.currentCh = get(hObject,'Value');
set(handles.ChText,'String',['Ch ' num2str(handles.currentCh)]);

zebraguidata(hObject, handles);
plotWideWindow (hObject, handles);
plotMultiBand(hObject, handles);
plotSpectrogram (hObject, handles);
if get(handles.powerSpectraCk,'Value')
    startSpectraPlot (handles);
end
doInsets (hObject, handles);

function chSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function trialSlider_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if (handles.fileInNum >= 1 && handles.nTrials > 1)
    handles.currentTrial = int16(get(hObject,'Value'));
    set(handles.trialText,'String',['Trial ' num2str(handles.currentTrial)]);
    zebraguidata(hObject, handles);
    plotWideWindow (hObject, handles);
    plotMultiBand(hObject, handles);
    plotSpectrogram (hObject, handles);
    if get(handles.powerSpectraCk, 'Value')
        startSpectraPlot (handles);
    end
    doInsets (hObject, handles);
end

function trialSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function decimateFlag_Callback(hObject, eventdata, handles)

function plotMeanFlag_Callback(hObject, eventdata, handles)
if (handles.fileInNum>0), plotWideWindow (hObject, handles)
    plotMultiBand(hObject, handles);
    plotSpectrogram (hObject, handles)
    if get(handles.powerSpectraCk, 'Value')
        startSpectraPlot (handles);
    end    
end

function envSPGfrom_Callback(hObject, eventdata, handles)

function envSPGfrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function envSPGto_Callback(hObject, eventdata, handles)

function envSPGto_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SPGfs_Callback(hObject, eventdata, handles)

function SPGfs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filterBPpower_Callback(hObject, eventdata, handles)

function lowerDB_Callback(hObject, eventdata, handles)
handles.SPGlowDB = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);
if (handles.fileInNum>0), plotSpectrogram (hObject, handles)
end

function lowerDB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function upperDB_Callback(hObject, eventdata, handles)
handles.SPGhighDB = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);
if (handles.fileInNum>0), plotSpectrogram (hObject, handles)
end

function upperDB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tmaxTxt_Callback(hObject, eventdata, handles)
handles.tmax = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);
if (handles.tmin < handles.tmax), plotWideWindow (hObject, handles);
    plotMultiBand(hObject, handles);
    plotSpectrogram (hObject, handles);
end 

function tmaxTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tminTxt_Callback(hObject, eventdata, handles)
handles.tmin = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);
if (handles.tmin < handles.tmax), plotWideWindow (hObject, handles);
    plotMultiBand(hObject, handles);
    plotSpectrogram (hObject, handles);
end

function tminTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function baselineFrom_Callback(hObject, eventdata, handles)

function baselineFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function baselineTo_Callback(hObject, eventdata, handles)

function baselineTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tResetBtn_Callback(hObject, eventdata, handles)
handles.tmin = 0;
set(handles.tminTxt,'String',handles.tmin);
handles.tmax = handles.dtaLen*handles.sp;
set(handles.tmaxTxt,'String',handles.tmax);
zebraguidata(hObject, handles);
plotWideWindow (hObject, handles);
plotMultiBand(hObject, handles);
plotSpectrogram (hObject, handles);

function excludeChFlag_Callback(hObject, eventdata, handles)

function excludeCh_Callback(hObject, eventdata, handles)
 
function excludeCh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotMultiCh_Callback(hObject, eventdata, handles)
if (handles.fileInNum>0), plotWideWindow (hObject, handles)
    plotSpectrogram (hObject, handles)
    if get(handles.powerSpectraCk, 'Value')
        startSpectraPlot (handles);
    end    
end

function respFrom_Callback(hObject, eventdata, handles)

function respFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function respTo_Callback(hObject, eventdata, handles)

function respTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTrCk_Callback(hObject, eventdata, handles)
if (handles.fileInNum>0), plotWideWindow (hObject, handles)
    plotMultiBand(hObject, handles);
    plotSpectrogram (hObject, handles)
    if get(handles.powerSpectraCk, 'Value')
        startSpectraPlot (handles);
    end    
end

function autoOffCk_Callback(hObject, eventdata, handles)
if (handles.fileInNum>0), plotWideWindow (hObject, handles)
end

function dataFormatPanel_SelectionChangeFcn(hObject, eventdata, handles)
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object

leak = str2double(get(handles.leakageLP,'String'));
flag = get (hObject,'Tag');
if strcmpi(flag,'lauraEEG')
    setSP (hObject, handles, 512);
    setCh (hObject, handles, 2);
    setGain (hObject, handles, 5000);
    set (handles.EEGck,'Value',1);    
    handles.SPGlowDB = -80;
    set(handles.lowerDB,'String',handles.SPGlowDB);
    handles.SPGhighDB = -30;
    set(handles.upperDB,'String',handles.SPGhighDB);
    
    % set the correct spectral bands
    handles.frequencyBand = [0 4; 0.5 4; 12 16; 30 80; 100 0.5/handles.sp; 0 leak];
end

if strcmpi(flag,'MeyerVepRB')
    % set all required parameters
    %setSP (hObject, handles, 4096);
    setSP (hObject, handles, 2048/str2double(get(handles.SweepTimeEdit,'string'))); %gab
    handles = zebraguidata(hObject);
    setCh (hObject, handles, 3);
    handles = zebraguidata(hObject);
    setGain (hObject, handles, 1000);
    handles = zebraguidata(hObject);
    set (handles.microVck,'Value',1);
    set (handles.EEGck,'Value',1);    
    handles.SPGlowDB = -30; %gab
    set(handles.lowerDB,'String',handles.SPGlowDB);
    handles.SPGhighDB = -10; %gab
    set(handles.upperDB,'String',handles.SPGhighDB);
    
    % set the correct spectral bands
    handles.frequencyBand = [0 4; 0.5 4; 12 16; 30 80; 100 0.5/handles.sp; 0 leak];
    
      %GABRIELE 2018/01/27
    set(handles.spEdit,'Enable','off');
    set(handles.SweepTimeEdit,'visible','on');
    set(handles.SweepTimeTxt,'visible','on');
else
    if strcmpi(eventdata.OldValue.Tag,'meyerVepRB') %meyerVEP has been deselected
        set(handles.spEdit,'Enable','on');
        set(handles.SweepTimeEdit,'visible','off');
        set(handles.SweepTimeTxt,'visible','off');
    end
end

if strcmpi(flag,'MeyerSspRB')
    % set all required parameters
    setSP (hObject, handles, 32768);
    handles = zebraguidata(hObject);
    setCh (hObject, handles, 2);
    handles = zebraguidata(hObject);
    set (handles.microVck,'Value',1);
    setGain (hObject, handles, 1000);
    handles = zebraguidata(hObject);
    set (handles.EEGck,'Value',1);    
    handles.SPGlowDB = -80;
    set(handles.lowerDB,'String',handles.SPGlowDB);
    handles.SPGhighDB = -30;
    set(handles.upperDB,'String',handles.SPGhighDB);

    % this mode does not perform spectral analysis
end
% now update the edit fields
set (handles.lpTo,'String',handles.frequencyBand (1,2));
set (handles.alphaFrom,'String',handles.frequencyBand (2,1));
set (handles.alphaTo,'String',handles.frequencyBand (2,2));
set (handles.betaFrom,'String',handles.frequencyBand (3,1));
set (handles.betaTo,'String',handles.frequencyBand (3,2));
set (handles.gammaFrom,'String',handles.frequencyBand (4,1));
set (handles.gammaTo,'String',handles.frequencyBand (4,2));
set (handles.hpFrom,'String',handles.frequencyBand (5,1));

zebraguidata(hObject, handles);

function internalSynchCk_Callback(hObject, eventdata, handles)

function synchCh_Callback(hObject, eventdata, handles)

function synchCh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ymaxTxt_Callback(hObject, eventdata, handles)
%        str2double(get(hObject,'String')) returns contents of ymaxTxt as a double
handles.ymax = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);
if (handles.tmin < handles.tmax), plotWideWindow (hObject, handles);
    plotMultiBand(hObject, handles);
    plotSpectrogram (hObject, handles);
end

function ymaxTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yminTxt_Callback(hObject, eventdata, handles)

function yminTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function autoScaleCk_Callback(hObject, eventdata, handles)

function saveCurrentCh_Callback(hObject, eventdata, handles)
% Save as an ascii file the current view of the wide window
% If multi plot is on, all channels are saved as separate columns,
% otherwise it is saved only the current channel.
% This can be used to save only a small portion of the record depending on
% the X limits of the visualisation.

fileOut = [handles.dir_in 'currentTrace.dat'];
% check what data to save (single trial or mean over trials, single Ch or multi Ch)
flagMd = get(handles.plotMeanFlag,'Value');
flagMulti = get(handles.plotMultiCh,'Value');
tmpOut = [];
% determine the initial and final index of the visualized data
lim = axis(handles.wide_plot_up);
ifrom = lim(1)/handles.sp + 1;
ito = lim(2)/handles.sp + 1;
if (ito>handles.dtaLen), ito = handles.dtaLen;
end
if flagMd
    % save all channels of the mean trace. Reserve one column for the
    % ordinates.
    tmpOut (2:handles.nCh+1,1:(ito-ifrom+1)) = handles.meanLFP (1:handles.nCh,ifrom:ito);        
else
    % save the spectra of all channels of the current trial
%    tmpOut (2:handles.nCh+1,1:(ito-ifrom+1)) = handles.workLFP (1:handles.nCh,handles.currentTrial,ifrom:ito);        
    tmpOut (2:handles.nCh+1,1:(ito-ifrom+1)) = handles.workLFP (1:handles.nCh,handles.currentTrial,ifrom:ito);        
end
tmpOut(1,:) = handles.x (1:(ito-ifrom+1));
% create the format descriptor
fmtSt='';
for ii=1:handles.nCh      % add 1 for the x axis.
   fmtSt=[fmtSt '%8.6f '] ;
end
fmtSt=[fmtSt '%8.6f\n'];       % start new line

fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut);
fclose(fid);

function doubleStimCk_Callback(hObject, eventdata, handles)

function setSP (hObject, handles, freq)
% This function set the sampling period, frequency and all the relative numers in the GUI
handles.sp = 1/freq;        % sampling period in s
handles.acqF = freq;        % sampling frequency in Hz
set (handles.spEdit,'String',handles.sp);
handles.frequencyBand(5,2) = 0.5*freq;  %GAB 22/05/2018
zebraguidata(hObject, handles);
    
function setGain (hObject, handles, gain)
handles.gainAmp = gain;
set (handles.gain,'String',handles.gainAmp);
zebraguidata(hObject, handles);

function setCh (hObject, handles,nCh)
handles.nCh = nCh;
set (handles.chEdit,'String',handles.nCh);
handles.currentCh = 1;
zebraguidata(hObject, handles);

function SPGwindow_Callback(hObject, eventdata, handles)

function SPGwindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SPGoverlap_Callback(hObject, eventdata, handles)

function SPGoverlap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ChronuxCk_Callback(hObject, eventdata, handles)

function SPGlogCk_Callback(hObject, eventdata, handles)
doInsets (hObject, handles);
plotSpectrogram (hObject, handles)

function SPGmode_SelectionChangeFcn(hObject, ~, handles)
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
flag = get (hObject,'Tag');
doInsets (hObject, handles);
plotSpectrogram (hObject, handles)

function ChronuxAutoConfigCk_Callback(hObject, eventdata, handles)

function ChronuxRW_Callback(hObject, eventdata, handles)
 
function ChronuxRW_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ChronuxOvl_Callback(hObject, eventdata, handles)

function ChronuxOvl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function recomputeSPG_Callback(hObject, eventdata, handles)
computeSpectrogram (hObject, handles);
handles = zebraguidata(hObject);
if (handles.nTrials>1), computeMeanOnTrials (hObject, handles);
end
handles = zebraguidata(hObject);
% warning: if trial >1 the means have not been recomputed!
plotSpectrogram (hObject, handles);
doInsets (hObject, handles);

function SPGfromFreq_Callback(hObject, eventdata, handles)

function SPGfromFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SPGtoFreq_Callback(hObject, eventdata, handles)

function SPGtoFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SPGresolution_Callback(hObject, eventdata, handles)

function SPGresolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function templateCk_Callback(hObject, eventdata, handles)

function templateFrom_Callback(hObject, eventdata, handles)

function templateFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function templateTo_Callback(hObject, eventdata, handles)

function templateTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function autoscaleYck_Callback(hObject, eventdata, handles)

function plotTemplateCk_Callback(hObject, eventdata, handles)

function saveSPG_Callback(hObject, eventdata, handles)
% Save as an ascii file the current view of the SPG, envelope, RMS power and E/I index
% If plotMean is on it saves the computed means otherwise it saves the current trace.
% This can be used to save only a small portion of the record depending on
% the X limits of the visualisation.

% check what data to save (single trial or mean over trials, single Ch or multi Ch)
flagMd = get(handles.plotMeanFlag,'Value');
%flagMulti = get(handles.plotMultiCh,'Value'); yet unusued 

% determine the initial and final index of the visualized data
lim = axis(handles.wide_plot_up);
timeFrom = lim(1);
timeTo = lim(2);

% process the RMS powers of the BP1-3, HP.
% since each file may have different time axis they are saved as separate
% files

outBPpowers (handles, 2, timeFrom, timeTo, 'BP1pw.dat');
outBPpowers (handles, 3, timeFrom, timeTo, 'BP2pw.dat');
outBPpowers (handles, 4, timeFrom, timeTo, 'BP3pw.dat');
outBPpowers (handles, 5, timeFrom, timeTo, 'HPpw.dat');

% process the E/I ratio first, the envelope second
SPGmin = handles.spgt(1);
SPGmax = handles.spgt(handles.spgl);
resolutionSPG = (SPGmax-SPGmin)/(handles.spgl-1);
SPGfrom = int32(timeFrom/resolutionSPG)+1;
SPGto =  int32(timeTo/resolutionSPG)+1;

if (SPGto>handles.spgl), SPGto = handles.spgl;
end
tmpOut = [];

if flagMd
    % save all channels of the mean E/I ratio. Reserve one column for the
    % time coordinate.
    tmpOut (2:handles.nCh+1,1:(SPGto-SPGfrom+1)) = handles.meanEIrat (1:handles.nCh,SPGfrom:SPGto);        
else
    % save the E/I ratio of all channels of the current trial
    tmpOut (2:handles.nCh+1,1:(SPGto-SPGfrom+1)) = handles.EIrat (1:handles.nCh,handles.currentTrial,SPGfrom:SPGto);        
end
c(1,:) = handles.spgt (:);
% create the format descriptor
fmtSt='';
for ii=1:handles.nCh      % add 1 for the x axis.
   fmtSt=[fmtSt '%10.8f '] ;
end
fmtSt=[fmtSt '%10.8f\n'];       % start new line
fileOut = [handles.dir_in 'EIratio.dat'];
fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut);
fclose(fid);

% process the envelope
if flagMd
    % save all channels of the mean E/I ratio. Reserve one column for the
    % time coordinate.
    tmpOut (2:handles.nCh+1,1:(SPGto-SPGfrom+1)) = handles.meanSpgPlot (1:handles.nCh,SPGfrom:SPGto);        
else
    % save the E/I ratio of all channels of the current trial
    tmpOut (2:handles.nCh+1,1:(SPGto-SPGfrom+1)) = handles.spgPlot (1:handles.nCh,handles.currentTrial,SPGfrom:SPGto);        
end

fileOut = [handles.dir_in 'SPGplot.dat'];
fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut);
fclose(fid);

function outBPpowers (handles, bpSel, timeFrom, timeTo, fileName);

% process BP1
i_from = floor((timeFrom-5*handles.BPpowerSP(bpSel)) / handles.BPpowerSP(bpSel));
if (i_from <=1), i_from = 1;
end
i_to = ceil((timeTo-5*handles.BPpowerSP(bpSel)) / handles.BPpowerSP(bpSel));
if (i_to > length(handles.BPpower)), i_to = length(handles.BPpower);
end

tmpOut = [];
fmtSt = [];
% the time axis is identical for all channels
tmpOut(1,1:i_to-i_from+1) = ...
    handles.BPpower (1,bpSel,1,handles.currentTrial,i_from:i_to);  
% iterate on all channels
for ii=1:handles.nCh      % add 1 for the x axis.
    tmpOut(ii+1,1:i_to-i_from+1) = ...
        handles.BPpower (2,bpSel,ii,handles.currentTrial,i_from:i_to);
end
% prepare for output
for ii=1:handles.nCh      % add 1 for the x axis.
   fmtSt=[fmtSt '%10.8f '] ;
end
fmtSt=[fmtSt '%10.8f\n'];       % start new line
fileOut = [handles.dir_in fileName];
fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut);
fclose(fid);

function saveBP_Callback(hObject, eventdata, handles)
% Save as an ascii file the current view of the wide window
% If multi plot is on, all channels are saved as separate columns,
% otherwise it is saved only the current channel.
% This can be used to save only a small portion of the record depending on
% the X limits of the visualisation.

% Data are saved with the following order (GAB):
% 1) Single trial, multi channels
% col 1. time   ( NOT  Full data (with notch filter, if applied) )
% col 2-5. the four band width
% repeat for all channels if the multi channel flag is selected

% 2) Multi trial, multi channels (GAB, not checked)
% December 8, 2017
% The button saves only the bandpassed data that are visualized and that refer 
% to the current channel. I have implemented two different functioning modes:
% 2.1) if Plot mean checkbox is selected, the program saves only the mean 
% 2.2) if Plot mean is checked the program saves all trials as separate
% columns. If more than one bandwidth is specifiedthey are alla grouped
% together.

fileOut = [handles.dir_in 'currentTraceBP.dat'];
% check what data to save (single trial or mean over trials, single Ch or multi Ch)
flagMd = get(handles.plotMeanFlag,'Value');
flagMulti = get(handles.plotMultiCh,'Value');
tmpOut = [];
% determine the initial and final index of the visualized data
lim = axis(handles.wide_plot_up);
ifrom = lim(1)/handles.sp + 1;
ito = lim(2)/handles.sp + 1;
if (ito>handles.dtaLen), ito = handles.dtaLen;
end
tmpOut(1,:) = handles.x (1:(ito-ifrom+1));
cntOut = 1;
ct = handles.currentTrial;
for i=1:handles.nCh
    % plot on all available channels.
    if (handles.currentCh == i || flagMulti)
        % this channel has to be saved out        
        if (handles.nTrials>1 && flagMd)
            % save the mean data and BP data. If flagMd is false save only
            % the current trace.
            cntOut=cntOut+1;
            dimMeanBP=size(handles.meanBP);
            %tmpOut (cntOut,1:(ifrom-ito+1)) = handles.meanPWS (1:handles.nCh,1:(ifrom-ito+1));
            tmpOut (cntOut:cntOut+dimMeanBP(1)-1,1:(ito-ifrom+1)) = handles.meanBP(:,handles.currentCh,ifrom:ito); %GAB
        else
            % save data from currentTrial
            cntOut=cntOut+1;
            tmpOut (cntOut,1:(ito-ifrom+1)) = handles.workLFP(i,ct,ifrom:ito);
            % save the band passed data
            tmpOut (cntOut+1:cntOut+4,1:(ito-ifrom+1)) = handles.bandPassed_LFP(1:4,i, ct, ifrom:ito);
            cntOut=cntOut+4;
        end
    end
end    
% create the format descriptor
fmtSt='';
for ii=1:cntOut-1
   fmtSt=[fmtSt '%10.8f '] ;
end
fmtSt=[fmtSt '%10.8f\n'];       % start new line

fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut);
fclose(fid);

function leakageCk_Callback(hObject, eventdata, handles)

function leakageLP_Callback(hObject, eventdata, handles)
handles.frequencyBand(6,2) = str2double(get(hObject,'String'));
zebraguidata(hObject, handles);

function leakageLP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function displayLeakCk_Callback(hObject, eventdata, handles)
doInsets (hObject, handles)
plotSpectrogram (hObject, handles)

function winPWleft_Callback(hObject, eventdata, handles)

function winPWleft_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function winPWright_Callback(hObject, eventdata, handles)

function winPWright_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filtOrder_Callback(hObject, eventdata, handles)

function filtOrder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scrollLeft_Callback(hObject, eventdata, handles)

delta = handles.tmax - handles.tmin; 
handles.tmin = handles.tmin - delta;
if handles.tmin<0, handles.tmin=0;
end
handles.tmax = handles.tmin + delta;
set(handles.tminTxt,'String',handles.tmin);
set(handles.tmaxTxt,'String',handles.tmax);
zebraguidata(hObject, handles);
plotWideWindow (hObject, handles);
plotMultiBand(hObject, handles);
plotSpectrogram (hObject, handles);

function scrollRight_Callback(hObject, eventdata, handles)

delta = handles.tmax - handles.tmin; 
handles.tmin = handles.tmax; 
handles.tmax = handles.tmax + delta;
set(handles.tminTxt,'String',handles.tmin);
set(handles.tmaxTxt,'String',handles.tmax);
zebraguidata(hObject, handles);
plotWideWindow (hObject, handles);
plotMultiBand(hObject, handles);
plotSpectrogram (hObject, handles);

function autoclear_Callback(hObject, eventdata, handles)

function singleUnitsThr_Callback(hObject, eventdata, handles)
% hObject    handle to singleUnitsThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleUnitsThr as text
%        str2double(get(hObject,'String')) returns contents of singleUnitsThr as a double

function singleUnitsThr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderBPpower_Callback(hObject, eventdata, handles)
handles.BP4power = get(hObject,'Value');
set(handles.BPtext,'String',['Select BP channel ' num2str(handles.BP4power)]);
zebraguidata(hObject, handles);

function sliderBPpower_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.BP4power = get(hObject,'Value');
zebraguidata(hObject, handles);

function skip1trialCk_Callback(hObject, eventdata, handles)

% --- Executes on button press in sendToClipboard.
function sendToClipboard_Callback(hObject, eventdata, handles)
h = ZebraExplore;
hgexport(h,'-clipboard');

function sc4_Callback(hObject, eventdata, handles)
% send to the clipboard 
scButtonOff(hObject, handles)

gcfPosition = get(gcf,'position');      % get the position of the ZE figure expressed in pixels
sc4Position = get(handles.inset_WB,'Position'); % get the position of 'inset_WB' expressed in normalized units
        % 'OuterPosition' includes the axes margin, labels et al. 
sc4margin = get(handles.inset_WB,'TightInset'); % get the position of 'inset_WB' expressed in normalized units

sc4PositionPixels (1) = (sc4Position(1)+0.5*sc4margin(1))*(gcfPosition(3));
sc4PositionPixels (2) = (sc4Position(2)-1.5*sc4margin(2))*(gcfPosition(4));
sc4PositionPixels (3) = (sc4Position(3)-0.5*sc4margin(1))*(gcfPosition(3));
sc4PositionPixels (4) = (sc4Position(4)+sc4margin(2))*(gcfPosition(4));
%imageData = screencapture(handles.inset_WB, [], 'clipboard');
imageData = screencapture(gcf, sc4PositionPixels, 'clipboard');
scButtonOn(hObject, handles)


% --- Executes on button press in sc5.
function sc5_Callback(hObject, eventdata, handles)
% send to the clipboard 
scButtonOff(hObject, handles)
gcfPosition = get(gcf,'position');      % get the position of the ZE figure expressed in pixels
sc5Position = get(handles.inset_HP,'Position'); % get the position of 'inset_WB' expressed in normalized units
        % 'OuterPosition' includes the axes margin, labels et al. 
sc5margin = get(handles.inset_HP,'TightInset'); % get the position of 'inset_WB' expressed in normalized units
sc4Position = get(handles.inset_WB,'Position'); % get the position of 'inset_WB' expressed in normalized units
sc4margin = get(handles.inset_WB,'TightInset'); % get the position of 'inset_WB' expressed in normalized units

sc5PositionPixels (1) = (sc4Position(1)+0.5*sc4margin(1))*(gcfPosition(3));
sc5PositionPixels (2) = (sc5Position(2)-1.5*sc5margin(2))*(gcfPosition(4));
sc5PositionPixels (3) = (sc4Position(3)-0.5*sc4margin(1))*(gcfPosition(3));
sc5PositionPixels (4) = (sc5Position(4)+sc5margin(2))*(gcfPosition(4));
%imageData = screencapture(handles.inset_WB, [], 'clipboard');
imageData = screencapture(gcf, sc5PositionPixels, 'clipboard');
scButtonOn(hObject, handles)

function sc45_Callback(hObject, eventdata, handles)
% send to the clipboard 
scButtonOff(hObject, handles)

gcfPosition = get(gcf,'position');      % get the position of the ZE figure expressed in pixels
sc5Position = get(handles.inset_HP,'Position'); % get the position of 'inset_WB' expressed in normalized units
        % 'OuterPosition' includes the axes margin, labels et al. 
sc5margin = get(handles.inset_HP,'TightInset'); % get the position of 'inset_WB' expressed in normalized units
sc4Position = get(handles.inset_WB,'Position'); % get the position of 'inset_WB' expressed in normalized units
sc4margin = get(handles.inset_WB,'TightInset'); % get the position of 'inset_WB' expressed in normalized units

sc5PositionPixels (1) = (sc4Position(1)+0.5*sc4margin(1))*(gcfPosition(3));
sc5PositionPixels (2) = (sc5Position(2)-1.5*sc5margin(2))*(gcfPosition(4));
sc5PositionPixels (3) = (sc4Position(3)-0.5*sc4margin(1))*(gcfPosition(3));
sc5PositionPixels (4) = (sc5Position(4)+sc5margin(2))*(gcfPosition(4))...
    +(sc4Position(4)+1.5*sc4margin(2))*(gcfPosition(4));
%imageData = screencapture(handles.inset_WB, [], 'clipboard');
imageData = screencapture(gcf, sc5PositionPixels, 'clipboard');
scButtonOn(hObject, handles)

function sc56_Callback(hObject, eventdata, handles)

scButtonOff(hObject, handles)

gcfPosition = get(gcf,'position');      % get the position of the ZE figure expressed in pixels
sc5Position = get(handles.inset_HP,'Position'); % get the position of 'inset_WB' expressed in normalized units
        % 'OuterPosition' includes the axes margin, labels et al. 
sc5margin = get(handles.inset_HP,'TightInset'); % get the position of 'inset_WB' expressed in normalized units
sc6Position = get(handles.spectraInset,'Position'); % get the position of 'inset_WB' expressed in normalized units
sc6margin = get(handles.spectraInset,'TightInset'); % get the position of 'inset_WB' expressed in normalized units
sc4Position = get(handles.inset_WB,'Position'); % get the position of 'inset_WB' expressed in normalized units
sc4margin = get(handles.inset_WB,'TightInset'); % get the position of 'inset_WB' expressed in normalized units

% remember that all screen captures are aligned to sc4
sc56PositionPixels (1) = (sc4Position(1)+0.5*sc4margin(1))*(gcfPosition(3));
sc56PositionPixels (2) = (sc6Position(2)-sc6margin(2))*(gcfPosition(4));
sc56PositionPixels (3) = (sc4Position(3)-0.5*sc4margin(1))*(gcfPosition(3));
sc56PositionPixels (4) = (sc5Position(4)+sc6margin(2))*(gcfPosition(4))...
    +(sc6Position(4)+sc6margin(2))*(gcfPosition(4));

imageData = screencapture(gcf, sc56PositionPixels, 'clipboard');
scButtonOn(hObject, handles)

function sc6_Callback(hObject, eventdata, handles)
scButtonOff(hObject, handles)
gcfPosition = get(gcf,'position');      % get the position of the ZE figure expressed in pixels
sc6Position = get(handles.spectraInset,'Position'); % get the position of 'inset_WB' expressed in normalized units
        % 'OuterPosition' includes the axes margin, labels et al. 
sc6margin = get(handles.spectraInset,'TightInset'); % get the position of 'inset_WB' expressed in normalized units
sc4Position = get(handles.inset_WB,'Position'); % get the position of 'inset_WB' expressed in normalized units
sc4margin = get(handles.inset_WB,'TightInset'); % get the position of 'inset_WB' expressed in normalized units

sc6PositionPixels (1) = (sc4Position(1)+0.5*sc4margin(1))*(gcfPosition(3));
sc6PositionPixels (2) = (sc6Position(2)-sc6margin(2))*(gcfPosition(4));
sc6PositionPixels (3) = (sc4Position(3)-0.5*sc4margin(1))*(gcfPosition(3));
sc6PositionPixels (4) = (sc6Position(4)+1.5*sc6margin(2))*(gcfPosition(4));
imageData = screencapture(gcf, sc6PositionPixels, 'clipboard');
scButtonOn(hObject, handles)

function sc456_Callback(hObject, eventdata, handles)
scButtonOff(hObject, handles)
gcfPosition = get(gcf,'position');      % get the position of the ZE figure expressed in pixels
sc6Position = get(handles.spectraInset,'Position'); % get the position of 'inset_WB' expressed in normalized units
        % 'OuterPosition' includes the axes margin, labels et al. 
sc6margin = get(handles.spectraInset,'TightInset'); % get the position of 'inset_WB' expressed in normalized units
sc4Position = get(handles.inset_WB,'Position'); % get the position of 'inset_WB' expressed in normalized units
sc4margin = get(handles.inset_WB,'TightInset'); % get the position of 'inset_WB' expressed in normalized units

sc456PositionPixels (1) = (sc4Position(1)+0.5*sc4margin(1))*(gcfPosition(3));
sc456PositionPixels (2) = (sc6Position(2)-sc6margin(2))*(gcfPosition(4));
sc456PositionPixels (3) = (sc4Position(3)-0.5*sc4margin(1))*(gcfPosition(3));
%sc456PositionPixels (4) = (sc4Position(4)+sc4margin(2)+(sc4Position(2)-sc6Position(2)+sc6margin(2)))...
sc456PositionPixels (4) = (sc4Position(4)+(sc4Position(2)-sc6Position(2)+0.5*sc6margin(2)))...    
    *(gcfPosition(4));
imageData = screencapture(gcf, sc456PositionPixels, 'clipboard');
scButtonOn(hObject, handles)

function scButtonOn(hObject, handles)
set (handles.sc4,'visible','on');
set (handles.sc5,'visible','on');
set (handles.sc6,'visible','on');
set (handles.sc45,'visible','on');
set (handles.sc56,'visible','on');
set (handles.sc456,'visible','on');

function scButtonOff(hObject, handles)
set (handles.sc4,'visible','off');
set (handles.sc5,'visible','off');
set (handles.sc6,'visible','off');
set (handles.sc45,'visible','off');
set (handles.sc56,'visible','off');
set (handles.sc456,'visible','off');

function screenCaptureBtn_Callback(hObject, eventdata, handles)
flag = get(hObject,'Value');
if flag, scButtonOn(hObject, handles);
else
    scButtonOff(hObject, handles)
end

function powerSpectraCk_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    % open the child process
    handles.spectre = ZebraSpectre(hObject,handles);  % handle to the power spectra window
    handles.spectreData = zebraguidata(handles.spectre);
else
    close(handles.spectre);    
end
zebraguidata(hObject, handles);

function autoSavePSck_Callback(hObject, eventdata, handles)

function analyzePWdistributionCk_Callback(hObject, eventdata, handles)

handles = zebraguidata(hObject);  % refresh handles. It should not be necessary, but better be safe than sorry.
if get(hObject,'Value')
    % open the child process
    handles.analyzePW = analyzePWdistribution(gcf,handles);      % handle to the child figure
                                            % arguments: parent object and
                                            % parent handles
    handles.analyzePWdata = zebraguidata(handles.analyzePW);    % and its controls
else
    % close the child process
    close(handles.analyzePW);
    handles.analyzePW = [];
end
zebraguidata(hObject, handles);

function autoSavePWck_Callback(hObject, eventdata, handles)

function startSpectraPlot (handles)
% this function causes the refresh of the spectra plots by activating a
% callback in the ZebraSpectre figure

%GAB 11/03/2018
%   if the condition is not true, it means that Spectre window is open but
%   spectra have not been calculated (and so plotted) already.
if ~isempty(handles.spectreData.powerSpectra.Children)   
    handles.spectreData = zebraguidata(handles.spectre);
    ZebraSpectre('maxFreq_Callback',handles.spectreData.maxFreq,[],handles.spectreData)
end

function triggerMaxLength_Callback(hObject, eventdata, handles)

function triggerMaxLength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LPenableCK_Callback(hObject, eventdata, handles)

function HPenableCK_Callback(hObject, eventdata, handles)

function BP1enableCK_Callback(hObject, eventdata, handles)

function BP3enableCK_Callback(hObject, eventdata, handles)

function BP2enableCK_Callback(hObject, eventdata, handles)

function AutoOpening_Callback(hObject, eventdata, handles)

handles = zebraguidata(hObject);  % refresh handles. It should not be necessary, but better be safe than sorry.
if get(hObject,'Value')
    % open the child process
    handles.auto = ZebraAuto(gcf,handles);      % handle to the child figure
                                            % arguments: parent object and
                                            % parent handles
    handles.autoData = zebraguidata(handles.auto);    % and its controls
else
    % close the child process
    close(handles.auto);
    handles.autoData = [];
end
zebraguidata(hObject, handles);

function SweepTimeEdit_Callback(hObject, eventdata, handles)

%GABRIELE 2018/01/27: useful for setting trial time in MeyerVEP mode

handles.sp = str2double(get(hObject,'String'))/2048;
set(handles.spEdit,'string',num2str(handles.sp));
handles.acqF = 1/handles.sp;
zebraguidata(hObject, handles);


function SweepTimeEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TrRej_ck_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    handles.trialRej=trialRejection;
    handles.trialRejData=zebraguidata(handles.trialRej);
else
    close(handles.trialRej)
    handles.trialRejData=[];
end
zebraguidata(hObject,handles);

function multiChAnalyzer_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    % open the child process
    handles.xCorr = multiChannelAnalyzer(hObject,handles);  % handle to the multi channel analyzer
    handles.xCorrData = zebraguidata(handles.xCorr);
else
    close(handles.xCorr);    
end
zebraguidata(hObject, handles);

function singleUnitExtraction_ck_Callback(hObject, eventdata, handles)
% hObject    handle to singleUnitExtraction_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Didi October 2018: open single units program to extract single units from
% LFP
 if get(hObject, 'Value')
     % open the single unit program
     handles.SU = singleunits2(hObject, handles);
     handles.SUData = zebraguidata(handles.SU);
     
 else
     close(handles.SU);
 end
 zebraguidata(hObject, handles);



% --- Executes on button press in templExp_bt.
function templExp_bt_Callback(hObject, eventdata, handles)
% hObject    handle to templExp_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

templ = handles.templ;
uisave('templ','template.mat')


% --- Executes on button press in templImp_bt.
function templImp_bt_Callback(hObject, eventdata, handles)
% hObject    handle to templImp_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
templTemp = uigetfile;
handles.templ = templTemp;
zebraguidata(hObject,handles)

% --- Executes on button press in templRecalc_bt.
function templRecalc_bt_Callback(hObject, eventdata, handles)
% hObject    handle to templRecalc_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in smooth_ck.
function smooth_ck_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_ck


% --- Executes on button press in ROI4template_bt.
function ROI4template_bt_Callback(hObject, eventdata, handles)
% hObject    handle to ROI4template_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'hRoi')
    delete(handles.hRoi);
end
handles.hRoi = drawrectangle(handles.wide_plot_up);
templUpdateFunction(handles);
zebraguidata(hObject,handles);

function templUpdateFunction(handles)
inf_i = round(handles.hRoi.Vertices(1,1)/handles.sp)+1; %index of the roi interval left edge
sup_i = round(handles.hRoi.Vertices(3,1)/handles.sp)+1; %index of the roi interval right edge
[~,p100_i] = min(handles.meanLFP(handles.currentCh,inf_i:sup_i));
p100_i = p100_i +inf_i-1;
p100_t = (p100_i-1)*handles.sp;
bef = str2double(handles.p100bef.String);
aft = str2double(handles.p100after.String);
handles.templateFrom.String = num2str(p100_t-bef);
handles.templateTo.String = num2str(p100_t+aft);



function p100after_Callback(hObject, eventdata, handles)
% hObject    handle to p100after (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p100after as text
%        str2double(get(hObject,'String')) returns contents of p100after as a double


% --- Executes during object creation, after setting all properties.
function p100after_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p100after (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p100bef_Callback(hObject, eventdata, handles)
% hObject    handle to p100bef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p100bef as text
%        str2double(get(hObject,'String')) returns contents of p100bef as a double


% --- Executes during object creation, after setting all properties.
function p100bef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p100bef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hit_ck.
function hit_ck_Callback(hObject, eventdata, handles)
% hObject    handle to hit_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hit_ck
if get(hObject,'Value')
    handles.responsometer = Responsometer(gcf,handles);
    %handles.responsometerData = zebraguidata(handles.responsometer);
else
    % close the child process
    close(handles.responsometer);
    %handles.responsometerData = [];
    handles.responsometer = [];
end
zebraguidata(hObject, handles);



function trTimeOffset_txt_Callback(hObject, eventdata, handles)
% hObject    handle to trTimeOffset_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trTimeOffset_txt as text
%        str2double(get(hObject,'String')) returns contents of trTimeOffset_txt as a double


% --- Executes during object creation, after setting all properties.
function trTimeOffset_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trTimeOffset_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveParams_btn.
function saveParams_btn_Callback(hObject, eventdata, handles)
% hObject    handle to saveParams_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj(ancestor(hObject,'figure','toplevel'),'Style','Edit','-or','Style','Checkbox');
tableMain = TabFromObj(h);
settingTab.main = tableMain;

spec = findobj('Tag','ZebraSpectre');
if max(size(spec))
    tableSpectre = TabFromObj(findobj(spec,'Style','Edit','-or','Style','Checkbox'));
    settingTab.spectre = tableSpectre;
end

auto = findobj('Tag','ZebraAuto');
if max(size(auto))
    tableAuto = TabFromObj(findobj(auto,'Style','Edit','-or','Style','Checkbox'));
    settingTab.auto = tableAuto;
end


uisave('settingTab','parameters.mat')


% --- Executes on button press in loadParams_btn.
function loadParams_btn_Callback(hObject, eventdata, handles)
% hObject    handle to loadParams_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiopen('*.mat');
for i = 1 : length(settingTab.main)
    i
    if i==24
        disp('eccoci')
    end
    eval(['set(handles.' settingTab.main{i,1} ',"String","' settingTab.main{i,2} '")'])
    eval(['set(handles.' settingTab.main{i,1} ',"Value",' num2str(settingTab.main{i,3}) ')'])
    try
        eval([settingTab.main{i,1} '_Callback(handles.' settingTab.main{i,1} ',[],handles)'])
        handles = zebraguidata(hObject);
    catch ME
        warning('off','backtrace')
        warning(['Error while executing ' settingTab.main{i,1} '_Callback'])
        warning(ME.message)
        disp(newline) 
    end
end

if isfield(settingTab, 'spectre')
    disp(['spectre debugging' newline])
    for i = 1 : length(settingTab.spectre)
        disp(i)
        eval(['set(handles.spectreData.' settingTab.spectre{i,1} ',''String'',''' settingTab.spectre{i,2} ''')'])
        eval(['set(handles.spectreData.' settingTab.spectre{i,1} ',''Value'',' num2str(settingTab.spectre{i,3}) ')'])
        try
            eval(['ZebraSpectre(''' settingTab.spectre{i,1} '_Callback'',handles.spectreData.' settingTab.spectre{i,1} ',[],handles.spectreData)'])
%             handles = zebraguidata(hObject); %this does not update local variables!
            handles.spectreData = zebraguidata(handles.spectre); %this updates only local variables;
%                                                             if something is changed in ZebraMain handles,
%                                                             it will be lost.
                                                           
        catch ME
            warning('off','backtrace')
            warning(['Error while executing ' settingTab.spectre{i,1} '_Callback'])
            warning(ME.message)
            disp(newline)           
        end
%         %%% debugging        
%         xxx = findobj('Tag','ZebraSpectre');
%         xxxHand = zebraguidata(xxx);
%         disp(xxxHand.parentObject)
%         %%%
    end
    disp([newline newline newline])

end

if isfield(settingTab, 'auto')
    disp(['auto debugging' newline])
    for i = 1 : length(settingTab.auto)
        eval(['set(handles.autoData.' settingTab.auto{i,1} ',''String'',''' settingTab.auto{i,2} ''')'])
        eval(['set(handles.autoData.' settingTab.auto{i,1} ',''Value'',' num2str(settingTab.auto{i,3}) ')'])
        eval(['ZebraAuto(''' settingTab.auto{i,1} '_Callback'',handles.autoData.' settingTab.auto{i,1} ',[],handles.autoData)'])
%         handles = zebraguidata(hObject); %this does not update local variables!
        handles.autoData = zebraguidata(handles.auto); %this updates only local variables;
%                                                     if something is changed in ZebraMain handles,
%                                                     it will be lost.
%         %%% debugging
%         xxx = findobj('Tag','ZebraAuto');
%         xxxHand = zebraguidata(xxx);
%         disp(xxxHand.parentObject)
%         %%%
    end
    ZebraAuto('Synch_Callback',handles.autoData.Synch,[],handles.autoData);
end

warning('on','backtrace')
% load handel
% sound(y,Fs)
disp('params loaded!')

function tab = TabFromObj(h)

tab = cell(length(h),3);
for i=1:length(h)
    tab{i,1} = h(i).Tag;
    tab{i,2} = h(i).String;
    tab{i,3} = h(i).Value;
end


% --- Executes on button press in concat_ck.
function concat_ck_Callback(hObject, eventdata, handles)
% hObject    handle to concat_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    handles.concat_btn.Enable = 'on';
else
    handles.concat_btn.Enable = 'off';
end

% Hint: get(hObject,'Value') returns toggle state of concat_ck

% --- Executes on button press in analyzeWhileConc_ck.
function analyzeWhileConc_ck_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeWhileConc_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analyzeWhileConc_ck
if hObject.Value
    handles.concat_ck.Value = 1;
end

% --- Executes on button press in concat_btn.
function concat_btn_Callback(hObject, eventdata, handles)
% hObject    handle to concat_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
concatenateTrials(hObject,handles)


function concatenateTrials(hObject,handles)
handles.LFP = handles.storedLFP;
handles.dtaLen = size(handles.LFP,3);
handles.nTrials = size(handles.LFP,2);
openFile_partII(hObject,handles);


% --- Executes on button press in clearConc_btn.
function clearConc_btn_Callback(hObject, eventdata, handles)
% hObject    handle to clearConc_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.storedLFP = [];
zebraguidata(hObject,handles)


% --- Executes on button press in pushbutton_evDetUpdate.
function pushbutton_evDetUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_evDetUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.evDetData = zebraguidata (handles.evDetection);


% --- Executes on button press in checkbox_sleepscoring.
function checkbox_sleepscoring_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sleepscoring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Didi October 2018: open the sleep scoring program to score NREM sleep
% episodes in the australian data, based on delta/theta ratio, and
% determine power in delta band in these episodes

if get(hObject, 'Value')
     % open the program
     handles.S_SC = sleepscoring(hObject, handles);
     handles.S_SCData = zebraguidata(handles.S_SC);
     
 else
     close(handles.S_SC);
 end
 zebraguidata(hObject, handles);
