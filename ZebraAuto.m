function varargout = ZebraAuto(varargin)
% ZEBRAAUTO MATLAB code for ZebraAuto.fig
%      ZEBRAAUTO, by itself, creates a new ZEBRAAUTO or raises the existing
%      singleton*.
%
%      H = ZEBRAAUTO returns the handle to a new ZEBRAAUTO or the handle to
%      the existing singleton*.
%
%      ZEBRAAUTO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZEBRAAUTO.M with the given input arguments.
%
%      ZEBRAAUTO('Property','Value',...) creates a new ZEBRAAUTO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ZebraAuto_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ZebraAuto_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ZebraAuto

% Last Modified by GUIDE v2.5 04-Nov-2018 16:02:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ZebraAuto_OpeningFcn, ...
                   'gui_OutputFcn',  @ZebraAuto_OutputFcn, ...
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


% --- Executes just before ZebraAuto is made visible.
function ZebraAuto_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ZebraAuto (see VARARGIN)

% Choose default command line output for ZebraAuto
handles.output = hObject;

handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

handles.autoMode=0;

failThreshold_edit_Callback(handles.failThreshold_edit, [], handles)
handles = zebraguidata(hObject);

spectraFrom_edit_Callback(handles.spectraFrom_edit,[],handles)
handles = zebraguidata(hObject);

spectraTo_edit_Callback(handles.spectraTo_edit,[],handles)
handles = zebraguidata(hObject);

handles.fileNameIn=handles.fileNameIn_edit.String;

handles.selCh = str2double(get(handles.channel_txt,'String'));

% Update handles structure
zebraguidata(hObject, handles);

% UIWAIT makes ZebraAuto wait for user response (see UIRESUME)
% uiwait(handles.ZebraAuto);


% --- Outputs from this function are returned to the command line.
function varargout = ZebraAuto_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in DialBox_bt.
function DialBox_bt_Callback(hObject, eventdata, handles)
% hObject    handle to DialBox_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.dataSet={'0.105', '0.123', '0.141', '0.173', '0.300', '0.555'}; %contrast values
% handles.nCont=length(handles.dataSet);
handles.nCont = str2double(handles.nConds_txt.String);

handles.dataSet={};
conds = 1:handles.nCont;
for i = 1:length(conds)
handles.dataSet{i}=num2str(conds(i));
end

prompt=cell(handles.nCont,1);
preset=cell(handles.nCont,1);
for i=1:handles.nCont
    prompt{i}=['specify folders for contrast level ' handles.dataSet{i}];
    if isfield(handles,'dlgIn')
        if ~isempty(handles.dlgIn{i})
            preset{i}=handles.dlgIn{i};
        else
            preset{i}='';
        end
    else
        preset{i}='';
    end
end

title='specify folders for analysis (space separated)';
nlines=1;
handles.dlgIn=inputdlg(prompt, title, nlines, preset, 'on');
if isempty(handles.dlgIn)
    return
end

for i=1:length(handles.dlgIn)
    handles.folders{i}=str2num(handles.dlgIn{i});
end

zebraguidata(hObject,handles);


function mainFolder_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mainFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folderIn = hObject.String;
zebraguidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function mainFolder_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fileNameIn_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fileNameIn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fileNameIn = hObject.String;
zebraguidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function fileNameIn_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileNameIn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outFileName_Callback(hObject, eventdata, handles)
% hObject    handle to outFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fileNameOut=hObject.String;
zebraguidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function outFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Synch.
function Synch_Callback(hObject, eventdata, handles)
% hObject    handle to Synch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.parentHandles=zebraguidata(handles.parentObject);
handles.parentHandles.autoData=[];
zebraguidata(hObject,handles);


% --- Executes on button press in table_ck.
function table_ck_Callback(hObject, eventdata, handles)
% hObject    handle to table_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in spectra_ck.
function spectra_ck_Callback(hObject, eventdata, handles)
% hObject    handle to spectra_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectra_ck


% --- Executes on button press in spg_ck.
function spg_ck_Callback(hObject, eventdata, handles)
% hObject    handle to spg_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spg_ck


% --- Executes on button press in go_bt.
function go_bt_Callback(hObject, eventdata, handles)
% hObject    handle to go_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



%for checking where eventual errors occur, i want to obtain back here the h
%and i of the automatization2 loop. so i define hObj and iObj and pass them
%by reference: they're updated in real time!
% hObj=HandleObjectGab('hObj', 0);
% iObj=HandleObjectGab('iObj', 0);


handles.parentHandles.path_in=handles.folderIn;
handles.parentHandles.selection=handles.fileNameIn;

nCond=length(handles.folders);

if handles.tableAvg_ck.Value
    outTab=struct('AmpAvg',[],'AmpSEM',[],'TpeakAvg',[],'TpeakSEM',[],'BasePWAvg',[],'BasePWSEM',[],...
        'RespPWAvg',[],'RespPWSEM',[],'RatioPWAvg',[],'RatioPWSEM',[],'TemplPWAvg',[],'TemplPwSEM',[]);
else
    outTab=struct('Table',[]);
end
    
s=[];

if handles.spg_ck.Value
    spg=struct('spgMat',[],'spgT',[],'spgW',[]);
end

for i_Cond=1:nCond
    currFolds=handles.folders{i_Cond}; %currFolds=folders of the current condition
    nTraces=length(currFolds);
    skip = 0; %useful for "concatenate" condition. If skip = 1, skips the storage of data form a trace
    if nTraces>0
        for i_Traces=1:nTraces
                        
            iCond = i_Cond
            iTraces = i_Traces
            
            if currFolds(iTraces)>9 %we need this because folder names as Dir_03 cannot be returned, i.g.: num2string(03)=3 
                handles.parentHandles.path_in=[handles.folderIn '\Dir_' num2str(currFolds(iTraces))];
            else
                handles.parentHandles.path_in=[handles.folderIn '\Dir_0' num2str(currFolds(iTraces))];
            end
            
            zebraguidata(handles.parentHandles.ZebraMainFig,handles.parentHandles);
            ZebraExplore('select_data_Callback',handles.parentHandles.select_data,[],handles.parentHandles);
                        
            handles.parentHandles=zebraguidata(handles.parentHandles.ZebraMainFig);
            if isfield(handles.parentHandles,'autoData')
                handles.parentHandles.autoData=[];
            end %this was to avoid exponential propagation of handles copies from parents to children and vice versa!!
            
            %for the moment always select ch1
            %handles.selCh=1;
            if handles.selCh ~= 1
                handles.parentHandles.chSlider.Value = handles.selCh;
                zebraguidata(handles.parentObject,handles.parentHandles);
                ZebraExplore('chSlider_Callback',handles.parentHandles.chSlider,[],handles.parentHandles);
                handles.parentHandles = zebraguidata(handles.parentObject);
            end
            %isn't there an easier way?
            
            %2018/11/04: if "concatenation is selected", just open all the
            %"nTraces" files and extract data only at the end of a "condition".
            if handles.parentHandles.concat_ck.Value
                if iTraces < nTraces
                    skip = 1;
                else
                    skip = 0;
                    iTraces = 1; %so that the result is stored as it was the "fist trace"
                    %cooncatenate trials and compute everything in ZebraMain
                    ZebraExplore('concat_btn_Callback',handles.parentHandles.concat_btn,[],handles.parentHandles); 
                    handles.parentHandles = zebraguidata(handles.parentObject);
                end
            end
            
            if ~skip
                %% tables
                if handles.table_ck.Value
                    CheckCh=cell2mat(handles.parentHandles.zebraSaysData.outEvochedResponses.Data(:,2))==handles.selCh;
                    tmpTab=handles.parentHandles.zebraSaysData.outEvochedResponses.Data(CheckCh,:);
                    l=size(tmpTab);
                    if iTraces==1 %if we're in the 1st folder of a new dataset, clear all
                        if ~handles.parentHandles.concat_ck.Value
                            currCondTab=cell(l(1), l(2), nTraces); %let's prepare a table to store a table from each of the nTraces files.
                        else
                            currCondTab=cell(l(1), l(2), 1); %we only have 1 table summarizing all the concatenated trials
                        end
                    end

                    %it may happen that different traces have different number of epochs, so:
                    l_cct=size(currCondTab);
                    if l(1)==l_cct(1)
                        currCondTab(:,:,iTraces)=tmpTab;
                    else
                        if l(1)>l_cct(1) 
                            %fill previous columns with some NaN s
                            temp=currCondTab(end-1:end,:,:);
                            filling=NaN(l(1)-l_cct(1),l_cct(2), nTraces);
                            currCondTab=[currCondTab(1:end-2,:,:); num2cell(filling); temp];
                            currCondTab(:,:,iTraces)=tmpTab;
                        else
                            %copy the table but the last 2 rows, fill the gap, copy
                            %the last 2 rows.
                            currCondTab(1:l(1)-2,:,iTraces)=tmpTab(1:end-2,:);
                            currCondTab(l(1)-1:end-2,:,iTraces)=num2cell(NaN(l_cct(1)-l(1),l_cct(2),1));
                            currCondTab(end-1:end,:,iTraces)=tmpTab(end-1:end,:);
                        end
                    end
                    outTab(iCond).templPeakSgolay(iTraces) = handles.parentHandles.templatePeakSigned(handles.selCh);
                end

                %% spectra
                if handles.spectra_ck.Value
                    %--set parameters on ZebraSpectre----------------------
                    if handles.winSpectra_ck.Value
                        handles.parentHandles = zebraguidata(handles.parentObject); %is this really necessary?
                        handles.parentHandles.spectreData.fromTtxt.String = num2str(handles.spectraFrom);
                        handles.parentHandles.spectreData.toTtxt.String = num2str(handles.spectraTo);
                        handles.parentHandles.spectreData.windowCk.Value = 1;
                        handles.parentHandles.spectreData.subMean_ck.Value = 1;
                        if handles.winSpectra2_ck.Value
                            handles.parentHandles.spectreData.fromT2txt.String = num2str(handles.spectraFrom2);
                            handles.parentHandles.spectreData.toT2txt.String = num2str(handles.spectraTo2);
                            handles.parentHandles.spectreData.window2Ck.Value = 1;
                        end
                        %zebraguidata(handles.parentHandles.spectre,handles.parentHandles.spectreData); %not necessary
                    end
                    ZebraSpectre('computeSpectra_Callback',handles.parentHandles.spectre,[],handles.parentHandles.spectreData)
                    handles.parentHandles = zebraguidata(handles.parentObject);
                    %zebraguidata(handles.parentHandles.spectre,handles.parentHandles.spectreData);

                    %--Initialize variables at the beginning---------------------
                    if (iTraces==1 && iCond==1) %ho messo iTraces al posto di nTraces, che non aveva senso (2018/11/04)
                        if handles.spectraTr_ck.Value %save spectra of single traces
                            specTemp=cell(nCond,nTraces,2); %It's a cell array because different traces can have different number of trials.
                            if handles.winSpectra2_ck.Value %second window
                                specTemp2=cell(nCond,nTraces,2);
                            end
                        end
                        if handles.spectraAvg_ck.Value %save mean power spectra
                            specTempMean=zeros(nCond,nTraces,2,length(handles.parentHandles.pwsf));
                            if handles.winSpectra2_ck.Value %second window
                                specTempMean2=zeros(nCond,nTraces,2,length(handles.parentHandles.pwsf2));
                            end
                        end
                    end

                    %--Store spectra-------------------------------------------
                    if handles.spectraAvg_ck.Value %save mean power spectra
                        specTempMean(iCond,iTraces,1,:)=handles.parentHandles.pwsf; %frequency vector
                        specTempMean(iCond,iTraces,2,:)=handles.parentHandles.meanPWS(handles.selCh,:);
                        if handles.winSpectra2_ck.Value %second window
                            specTempMean2(iCond,iTraces,1,:)=handles.parentHandles.pwsf2; %frequency vector
                            specTempMean2(iCond,iTraces,2,:)=handles.parentHandles.meanPWS2(handles.selCh,:);
                        end
                    end
                    if handles.spectraTr_ck.Value %save spectra of single trials
                        specTemp{iCond,iTraces,1}=handles.parentHandles.pwsf;
                        specTemp{iCond,iTraces,2}=handles.parentHandles.powerS(handles.selCh,:,:);
                        if handles.winSpectra2_ck.Value %second window
                            specTemp2{iCond,iTraces,1}=handles.parentHandles.pwsf2;
                            specTemp2{iCond,iTraces,2}=handles.parentHandles.powerS2(handles.selCh,:,:);
                        end
                    end
                end
                %% spg
                if handles.spg_ck.Value
                    flagLeak = get (handles.parentHandles.displayLeakCk,'Value');
                    if iTraces == 1
                        if (handles.parentHandles.deleakedFlag && flagLeak)
                            dim = size(handles.parentHandles.meanSpgDeleaked(handles.selCh,:,:));
                        else
                            dim = size(handles.parentHandles.meanSpg(handles.selCh,:,:));
                        end
                        dim2 = size(handles.parentHandles.spgt);
                        dimMin = min(dim(3),dim2(2));
                        dim(3) = dimMin;

                    else
                        dim = size(handles.spg(iCond).spgMat(iTraces-1,:,:));
                        if (handles.parentHandles.deleakedFlag && flagLeak)
                            dim1 = size(handles.parentHandles.meanSpgDeleaked(handles.selCh,:,:));
                        else
                            dim1 = size(handles.parentHandles.meanSpg(handles.selCh,:,:));
                        end

                        %---
                        dim2 = size(handles.parentHandles.spgt);
                        dimMin = min(dim1(3),dim2(2));
                        %---

                        if dimMin<dim(3)
                            handles.spg(iCond).spgMat = handles.spg(iCond).spgMat(:,1:dim1(2),1:dimMin);
                            handles.spg(iCond).spgT = handles.spg(iCond).spgT(:,1:dimMin);
                            dim = size(handles.spg(iCond).spgMat(iTraces-1,:,:));
                        end
                    end
                    handles.spg(iCond).spgT(iTraces,:) = handles.parentHandles.spgt(1:dim(3));
                    handles.spg(iCond).spgW(iTraces,:) = handles.parentHandles.spgw(1:dim(2));
                    if (handles.parentHandles.deleakedFlag && flagLeak)
                        handles.spg(iCond).spgMat(iTraces,:,:) = handles.parentHandles.meanSpgDeleaked(handles.selCh,1:dim(2),1:dim(3));
                    else
                        handles.spg(iCond).spgMat(iTraces,:,:) = handles.parentHandles.meanSpg(handles.selCh,1:dim(2),1:dim(3));
                    end

                    %check that time and freq are the same among trials
                    if iTraces > 1
                        checkT = isequal(handles.spg(iCond).spgT(iTraces,:),handles.spg(iCond).spgT(iTraces-1,:));
                        checkW = isequal(handles.spg(iCond).spgW(iTraces,:),handles.spg(iCond).spgW(iTraces-1,:));
                        if checkT * checkW == 0
                            error('SPG: time or frequency mismatch')
                        end
                    end
                end

                %% traces
                if handles.trace_ck.Value
                    %Traces can have different length, but they will be saved
                    %in the same matrix "trace" (trace*time). so:
                    %   1- if iTraces==1 we're storing the first trace of this
                    %       condition, so it's ok
                    %   2- if the current trace is longer than the ones already stored,
                    %       save only part of it.
                    %   3- if the current trace is shorted than the ones already
                    %       stored, shorten the matrix "trace" and then store
                    %       the trace

                    if iTraces == 1     
                        handles.traces(iCond).trace = [];
                        handles.traces(iCond).trace = handles.parentHandles.meanLFP(handles.selCh,:);
    %                     dimTrace = size(handles.parentHandles.meanLFP(handles.selCh,:));
                    else
                        dimTraces = size(handles.traces(iCond).trace);
                        dimThis = size(handles.parentHandles.meanLFP(handles.selCh,:));

                        if dimThis(2) >= dimTraces(2)
                            handles.traces(iCond).trace(iTraces,1:dimTraces(2)) = handles.parentHandles.meanLFP(handles.selCh,1:dimTraces(2));
                        else
                            handles.traces(iCond).trace = handles.traces(iCond).trace(:,1:dimThis(2));
                            handles.traces(iCond).trace(iTraces,1:dimThis(2)) = handles.parentHandles.meanLFP(handles.selCh,1:dimThis(2));
                        end
                    end

    %                 handles.traces(iCond).trace(iTraces,1:dimTraces(2)) = handles.parentHandles.meanLFP(handles.selCh,1:dimTrace(2));
                end    
            end
        end
        
        %% this is done at the end of each "condition"
        if handles.table_ck.Value
            tabSubset=currCondTab(:,5:10,:); %collect only numeric data
            if handles.parentHandles.doubleStimCk.Value
                tabSubset(end-1:end,6,:)=num2cell(0); %this is needed because the last 2 cells of the 6th colums are empty and cannot be corverted into numbers later
            else
                tabSubset(end,6,:)=num2cell(0);
            end

            %convert the cell array into a matrix
            finalDim_a=size(currCondTab); %final dimensions of matrix a
            dim_subset=size(tabSubset);
%             tabMat=zeros(dim_subset);
%             for i=1:nTraces %runs over "analysed folders" CI DEVO METTERE i NON iTraces!!!!!
%                 for j=1:6 %runs over currCondTab colums, that are supposed to be 6(ampl,time,4*power) in the current version of Zebra
%                     tabMat(:,j,i)=cell2mat(tabSubset(:,j,i));
%                     if j==1
%                         tabMat(:,j,i)=tabMat(:,j,i)*(-1);
% %                         if handles.parentHandles.doubleStimCk.Value
% %                             if c(end-1,j,i)<0          
% %                                 c(end-1,j,i)=0; 
% %                             end
% %                         end
% %                         if c(end,j,i)<0 
% %                             c(end,j,i)=0; 
% %                         end
%                     end
% 
%                 end
%             end
            tabMat = cell2mat(tabSubset);
            tabMat(:,1,:) = -1*tabMat(:,1,:); %this command works even if ntraces = 1, and so tabMat is a 2D-matrix and not a 3D-matrix
            
            
            %let's do some statistics if AverageConditions is checked:
            if handles.tableAvg_ck.Value
                if handles.parentHandles.doubleStimCk.Value
                    fromRow=finalDim_a(1)-1;
                else
                    fromRow=finalDim_a(1);
                end
                [outTab(iCond).AmpAvg,     outTab(iCond).AmpSEM    ] = MeanAndSEM(tabMat,1,fromRow,finalDim_a(1));
                [outTab(iCond).TpeakAvg,   outTab(iCond).TpeakSEM  ] = MeanAndSEM(tabMat,2,fromRow,finalDim_a(1));
                [outTab(iCond).BasePWAvg,  outTab(iCond).BasePWSEM ] = MeanAndSEM(tabMat,3,1,finalDim_a(1)-2);
                [outTab(iCond).RespPWAvg,  outTab(iCond).RespPWSEM ] = MeanAndSEM(tabMat,4,1,finalDim_a(1)-2);
                [outTab(iCond).RatioPWAvg, outTab(iCond).RatioPWSEM] = MeanAndSEM(tabMat,5,1,finalDim_a(1)-2);
                [outTab(iCond).TemplPWAvg, outTab(iCond).TemplPwSEM] = MeanAndSEM(tabMat,6,1,finalDim_a(1)-2);
            %else, let's save all the "tables" .
            else
                for i=1:nTraces
                    outTab(iCond).Table=tabMat;
                end
            end
        end
        
        if handles.failures_ck.Value
            dimC = size(tabMat); %c is the "table" matrix relative to this condition.
            if handles.parentHandles.doubleStimCk.Value
                l = dimC(1)-2;
            else
                l = dimC(1)-1;
            end
            
            fails = sum(tabMat(:,1,:) < -1*handles.failThreshold,1); %number of fails for each trace; previously all values haves been multiplied by -1!!
            nans = sum(isnan(tabMat(:,1,:)),1); %number of NaNs eventually introduced before
            actTrials = l - nans; %actual number of trials
            
            %overall stats for this condition
            totFails = sum(fails);
            totTrials = sum(actTrials);
            
            %finalFails is a matrix conditions * twice(nTraces+1)
            if iCond == 1
                finalFails = zeros(nCond,2*(nTraces+1));
                for i=1:dimC(3)
                    finalFails(iCond,2*i-1) = fails(i);
                    finalFails(iCond,2*i) = actTrials(i);
                end
            else  %adjust if different number of traces
                dim_prec = size(finalFails);
                
                if dim_prec(2) < 2*(dimC(3)+1) %if this condition has more traces, fill the previous matrix with zeroes
                    diff = 2*(dimC(3)+1) - dim_prec(2);
                    lastCol = finalFails(1:iCond-1,end-1:end);
                    finalFails(:,end-1:end) = [];
                    finalFails(:,end+1:end+2) = 0;
                    finalFails(1:iCond-1,end+1:end+2) = lastCol;
                end
                
                for i=1:dimC(3) %copy the results for the traces of this condition
                finalFails(iCond,2*i-1) = fails(i);
                finalFails(iCond,2*i) = actTrials(i);
                end
                    
                if dim_prec(2) > 2*(dimC(3)+1)  %if this conditions has few traces, fill it with zeroes  
                    for i=dimC(3)+1:dim_prec(2)
                    finalFails(iCond,2*i-1) = 0;
                    finalFails(iCond,2*i) = 0;
                    end
                end
            end
            finalFails(iCond,end-1) = totFails;
            finalFails(iCond,end) = totTrials;
        end
        
%         if handles.spectra_ck.Value
%             %for the moment just collect spectra
%         end
        
        if handles.spg_ck.Value
            %fare spg stessa lunghezza
%             if nTraces>1
%                 checkT=1;
%                 for i=2:nTraces
%                     checkT=checkT*isequal(handles.spg(1).spgT,handles.spg(i).spgT);
%                 end
%                 if ~checkT
%                     error(['sgpT are different: length(spgT(1))=' num2str(length(spg(1).spgT)) '; length(spgT(' num2str(i) '))=' num2str(length(spg(i).spgT))])
%                 end
%             end
            handles.spgOut(iCond).spgT = handles.spg(iCond).spgT(1,:);
            handles.spgOut(iCond).spgW = handles.spg(iCond).spgW(1,:);
            handles.spgOut(iCond).spgMat = squeeze(mean(handles.spg(iCond).spgMat,1));
        end
        
        if handles.trace_ck.Value
            %fare tracce stessa lunghezza
            if iCond == 1
                l = length(handles.traces(iCond).trace(iTraces,:));
            else
                %if current condition has longer trace, then it is limited by
                %stopping at l;
                %if current condition has shorter trace, then shorten the
                %meanTraces and update l
                l = length(handles.meanTraces(iCond-1,:));
                l1 = length(handles.traces(iCond).trace(iTraces,:));
                if l1 < l
                    handles.meanTraces = handles.meanTraces(:,1:l1);
                    l = length(handles.meanTraces(iCond-1,:));
                end
            end
            
            tmp = mean(handles.traces(iCond).trace(iTraces,:),1);
            handles.meanTraces(iCond,:) = tmp(1:l);
        end   
        if handles.parentHandles.concat_ck.Value
            ZebraExplore('clearConc_btn_Callback',handles.parentHandles.clearConc_btn,[],handles.parentHandles)
            handles.parentHandles = zebraguidata(handles.parentObject);
        end
    else
        reminder = ['no selected folders for dataset ' handles.dataSet{iCond}]
    end
end

%% Save

if handles.table_ck.Value
    if handles.tableAvg_ck.Value
        z=struct2cell(outTab);
        dim=size(z);
        for i=1:dim(1)
            for j=1:dim(2)
                for k=1:dim(3)
                    if isempty(z{i,j})
                        z{i,j,k}=0;
                    end
                end
            end
        end
        z=squeeze(z)';
        z=cell2mat(z);
        save(handles.fileNameOut, 'z');
    else
        save(handles.fileNameOut, 'outTab');
    end
end

if handles.failures_ck.Value
    save([handles.fileNameOut(1:end-4) '_failures.mat'], 'finalFails');
end

if handles.spectra_ck.Value
    
    
    if handles.spectraAvg_ck.Value %save mean power spectra
        MeanSpectra.spectra1 = specTempMean;
        if handles.winSpectra2_ck.Value %second window
            MeanSpectra.spectra2 = specTempMean2;
        end
        save([handles.fileNameOut(1:end-4) '_SpectraMean.mat'], 'MeanSpectra');
    end
    if handles.spectraTr_ck.Value %save spectra of single trials
        TrialsSpectra.spectra1 = specTemp;
        if handles.winSpectra2_ck.Value %second window
            TrialsSpectra.spectra2 = specTemp2;
        end
        save([handles.fileNameOut(1:end-4) '_SpectraTrials.mat'], 'TrialsSpectra');
    end
end

if handles.spg_ck.Value
%     for i=1:length(handles.dataSet)
%         save([handles.fileNameOut(1:end-4) '_spgT_' handles.dataSet{i} '.mat'  ], 'handles.spgOut(i).spgT');
%         save([handles.fileNameOut(1:end-4) '_spgW_' handles.dataSet{i} '.mat'  ], 'handles.spgOut(i).spgW');
%         save([handles.fileNameOut(1:end-4) '_spgMat_' handles.dataSet{i} '.mat'], 'handles.spgOut(i).spgMat');
%     end
    tmp2=handles.spgOut;
    save([handles.fileNameOut(1:end-4) '_spg.mat'], 'tmp2');
end

if handles.trace_ck.Value
    tmp1 = [];
    tmp1(:,2:nCond+1) = handles.meanTraces';
    tmp1(:,1) =(0:length(tmp1)-1)*handles.parentHandles.sp;
    save([handles.fileNameOut(1:end-4) '_meanTrace.mat'],'tmp1')

    for i=1:length(handles.traces)
        handles.traces(i).time = (0:length(handles.traces(i).trace)-1)*handles.parentHandles.sp;
    end
    tmp2 = handles.traces;
    save([handles.fileNameOut(1:end-4) '_allTraces.mat'],'tmp2')
end
load gong
sound(y,Fs)
disp('ho finito!')



% --- Executes on button press in autoMode_bt.
function autoMode_bt_Callback(hObject, eventdata, handles)
% hObject    handle to autoMode_bt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.parentHandles=zebraguidata(handles.parentObject);

if ~handles.autoMode
    handles.autoMode=1;
    handles.go_bt.Enable='on';
    
    handles.parentHandles.autoMode=1;
    handles.parentHandles.select_data.Enable='off';
    
    hObject.BackgroundColor='g';
    hObject.String='Auto Mode ON';
else
    handles.autoMode=0;
    handles.go_bt.Enable='off';
    
    handles.parentHandles.autoMode=0;
    handles.parentHandles.select_data.Enable='on';
    
    hObject.BackgroundColor='r';
    hObject.String='Auto Mode OFF';
end
zebraguidata(hObject,handles);
zebraguidata(handles.parentObject,handles.parentHandles);
    
    

% --- Executes on button press in tableAvg_ck.
function tableAvg_ck_Callback(hObject, eventdata, handles)
% hObject    handle to tableAvg_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tableAvg_ck


% --- Executes on button press in trace_ck.
function trace_ck_Callback(hObject, eventdata, handles)
% hObject    handle to trace_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trace_ck


% --- Executes on button press in failures_ck.
function failures_ck_Callback(hObject, eventdata, handles)
% hObject    handle to failures_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of failures_ck



function failThreshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to failThreshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of failThreshold_edit as text
%        str2double(get(hObject,'String')) returns contents of failThreshold_edit as a double
handles.failThreshold = str2double(handles.failThreshold_edit.String);
zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function failThreshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to failThreshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in winSpectra_ck.
function winSpectra_ck_Callback(hObject, eventdata, handles)
% hObject    handle to winSpectra_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of winSpectra_ck



function spectraFrom_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spectraFrom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectraFrom_edit as text
%        str2double(get(hObject,'String')) returns contents of spectraFrom_edit as a double
handles.spectraFrom = str2double(hObject.String);
zebraguidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function spectraFrom_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectraFrom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spectraTo_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spectraTo_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectraTo_edit as text
%        str2double(get(hObject,'String')) returns contents of spectraTo_edit as a double
handles.spectraTo = str2double(hObject.String);
zebraguidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function spectraTo_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectraTo_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function channel_txt_Callback(hObject, eventdata, handles)
% hObject    handle to channel_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_txt as text
%        str2double(get(hObject,'String')) returns contents of channel_txt as a double
handles.selCh = str2double(get(hObject,'String'));
zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function channel_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nConds_txt_Callback(hObject, eventdata, handles)
% hObject    handle to nConds_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nConds_txt as text
%        str2double(get(hObject,'String')) returns contents of nConds_txt as a double


% --- Executes during object creation, after setting all properties.
function nConds_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nConds_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in winSpectra2_ck.
function winSpectra2_ck_Callback(hObject, eventdata, handles)
% hObject    handle to winSpectra2_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of winSpectra2_ck



function spectraFrom2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spectraFrom2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectraFrom2_edit as text
%        str2double(get(hObject,'String')) returns contents of spectraFrom2_edit as a double
handles.spectraFrom2 = str2double(hObject.String);
zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function spectraFrom2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectraFrom2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spectraTo2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spectraTo2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spectraTo2_edit as text
%        str2double(get(hObject,'String')) returns contents of spectraTo2_edit as a double
handles.spectraTo2 = str2double(hObject.String);
zebraguidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function spectraTo2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectraTo2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in spectraTr_ck.
function spectraTr_ck_Callback(hObject, eventdata, handles)
% hObject    handle to spectraTr_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectraTr_ck


% --- Executes on button press in spectraAvg_ck.
function spectraAvg_ck_Callback(hObject, eventdata, handles)
% hObject    handle to spectraAvg_ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectraAvg_ck
