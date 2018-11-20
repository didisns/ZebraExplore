function varargout = ZebraSpectre(varargin)

%GAB, 16/05/2018: subtract the mean before computing the spectrum.

% ZEBRASPECTRE MATLAB code for ZebraSpectre.fig
%      ZEBRASPECTRE, by itself, creates a new ZEBRASPECTRE or raises the existing
%      singleton*.
%
%      H = ZEBRASPECTRE returns the handle to a new ZEBRASPECTRE or the handle to
%      the existing singleton*.
%
%      ZEBRASPECTRE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZEBRASPECTRE.M with the given input arguments.
%
%      ZEBRASPECTRE('Property','Value',...) creates a new ZEBRASPECTRE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ZebraSpectre_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ZebraSpectre_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ZebraSpectre

% Last Modified by GUIDE v2.5 04-Nov-2018 12:46:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ZebraSpectre_OpeningFcn, ...
                   'gui_OutputFcn',  @ZebraSpectre_OutputFcn, ...
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


% --- Executes just before ZebraSpectre is made visible.
function ZebraSpectre_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for ZebraSpectre
handles.output = hObject;
% the next two lines were payed by spilling (my) blood!  
handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

handles.whWinLeft = 0;      % to be used for spectral peak computation
handles.whWinRight = 0;
handles.whiteFlag = 0;

axes(handles.powerSpectra);
xlabel('');
ylabel('Power (Db)');

% Update handles structure
zebraguidata(hObject, handles);

% UIWAIT makes ZebraSpectre wait for user response (see UIRESUME)
% uiwait(handles.ZebraSpectre);


% --- Outputs from this function are returned to the command line.
function varargout = ZebraSpectre_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function minPower_Callback(hObject, eventdata, handles)
plotPowerSpectra(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minPower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxPower_Callback(hObject, eventdata, handles)
plotPowerSpectra(hObject, handles);

function maxPower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minFreq_Callback(hObject, eventdata, handles)
plotPowerSpectra(hObject, handles);

function minFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxFreq_Callback(hObject, eventdata, handles)
plotPowerSpectra(hObject, handles);

function maxFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function clearSpectre_Callback(hObject, eventdata, handles)
axes (handles.powerSpectra);
cla;

function normLowFreq_Callback(hObject, eventdata, handles)

function normLowFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function normHighFreq_Callback(hObject, eventdata, handles)

function normHighFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function saveSpectra_Callback(hObject, eventdata, handles)

% Save the displayed and whitened spectra in the current directory

handles.parentHandles = zebraguidata(handles.parentObject);    % this to refresh the local values of parentHandles

fileOut = [handles.parentHandles.dir_in 'PWS.dat'];
% check what data to save (single trial or mean over trials)
flagMd = get(handles.parentHandles.plotMeanFlag,'Value');
tmpOut = [];
if flagMd
    % save all channels of the mean spectra. Reserve one column for the
    % ordinates.
    tmpOut (2:handles.parentHandles.nCh+1,1:handles.parentHandles.pwl) = handles.parentHandles.meanPWS (1:handles.parentHandles.nCh,1:handles.parentHandles.pwl);        
else
    % save the spectra of all channels of the current trial
    tmpOut (2:handles.parentHandles.nCh+1,1:handles.parentHandles.pwl) = handles.parentHandles.powerS (1:handles.parentHandles.nCh,handles.parentHandles.currentTrial,1:handles.parentHandles.pwl);        
end
tmpOut(1,:) = handles.parentHandles.pwsf (:);
tmpOut(handles.parentHandles.nCh+2,:) = handles.whitenedPS (:);

% creat the format descriptor
fmtSt='';
for ii=1:handles.parentHandles.nCh+1      % add 2 for the x axis and the whitened spectra.
   fmtSt=[fmtSt '%12.10f '] ;
end
fmtSt=[fmtSt '%12.10f\n'];       % start new line

fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut);
fclose(fid);

%ENRICO 29/10/2018 ADDED THE MATLAB VARIABLE
pws.meanPWS = squeeze(handles.parentHandles.meanPWS(handles.parentHandles.currentCh,1:handles.parentHandles.pwl))';        
pws.trialsPWS = squeeze(handles.parentHandles.powerS(handles.parentHandles.currentCh,1:handles.parentHandles.currentTrial,1:handles.parentHandles.pwl));
pws.pwsf = handles.parentHandles.pwsf(:);
save([handles.parentHandles.dir_in 'PWS.mat'],'pws')

function leftLimitWh_Callback(hObject, eventdata, handles)

function leftLimitWh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rightLimitWh_Callback(hObject, eventdata, handles)

function rightLimitWh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rightLimitWh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function whiteOut_Callback(hObject, eventdata, handles)

% Compute the linear fit to the log log spectra and divide the spectra by
% this function. This process allows a better visualization of local peaks
% in the power spectra

handles.parentHandles = zebraguidata(handles.parentObject);    % this to refresh the local values of parentHandles

fileOut = [handles.parentHandles.dir_in 'PWS.dat'];
% check what data is displayied (single trial or mean over trials)
flagMd = get(handles.parentHandles.plotMeanFlag,'Value');
tmpOut = [];

% if flagMd
%     % Load the power of each channels in a temporary array.
%     spectra (2:handles.parentHandles.nCh+1,1:handles.parentHandles.pwl) = handles.parentHandles.meanPWS (1:handles.parentHandles.nCh,1:handles.parentHandles.pwl);        
% else
%     % save the spectra of all channels of the current trial
% %    tmpOut (2:handles.parentHandles.nCh+1,1:handles.parentHandles.pwl) = handles.parentHandles.powerS (1:handles.parentHandles.nCh,handles.parentHandles.currentTrial,1:handles.parentHandles.pwl);        
% end

ch = handles.parentHandles.currentCh;
tr = handles.parentHandles.currentTrial;

frequency(1,:) = handles.parentHandles.pwsf (:);
if get(handles.whitenByFreq,'Value'),
    % whitening by multiplication of the power by the frequency
    if flagMd,
        % whiteout the mean spectra
        spectra(1:handles.parentHandles.pwl) = handles.parentHandles.meanPWS(ch,1:handles.parentHandles.pwl);
    else
        % white out the current trial
        %handles.parentHandles.powerS (1:handles.parentHandles.nCh,handles.parentHandles.currentTrial,1:handles.parentHandles.pwl);
        spectra(1:handles.parentHandles.pwl) = handles.parentHandles.powerS (ch, tr, 1:handles.parentHandles.pwl);
    end
    % now divide element by element the power by the frequency
    handles.whitenedPS = 10*log10(spectra .* frequency);
else
    % whitening by linear fit of the log log spectra 
    freqMin = str2double(get (handles.leftLimitWh, 'String'));
    freqMax = str2double(get (handles.rightLimitWh, 'String'));
    freqMinI = 0;
    freqMaxI = 0;
    for i=1:handles.parentHandles.pwl
        % I know, this is lousy code...
        f = handles.parentHandles.pwsf(i);
        if f >= freqMin,
            if freqMinI == 0,
                freqMinI = i;
            end
        end 
        if f <= freqMax,
            freqMaxI = i;
        end
    end
    % now the freqMinI and freqMaxI contain the limits for the slope fitting
    % the computation is performed only for the current channel and trial

    fitF = [];
    fitF(1:freqMaxI-freqMinI+1,1:2) = 1;
    fitF(1:freqMaxI-freqMinI+1,2) = log10(handles.parentHandles.pwsf (freqMinI:freqMaxI));
    if flagMd,
        % whiteout the mean spectra
        fitS = 10*log10(handles.parentHandles.meanPWS (ch, freqMinI:freqMaxI));
        spectra = 10*log10(handles.parentHandles.meanPWS(ch,1:handles.parentHandles.pwl));
    else
        % white out the current trial
        %handles.parentHandles.powerS (1:handles.parentHandles.nCh,handles.parentHandles.currentTrial,1:handles.parentHandles.pwl);
        fitS(1:freqMaxI-freqMinI+1) = 10*log10(handles.parentHandles.powerS (ch, tr, freqMinI:freqMaxI));
        spectra(1:handles.parentHandles.pwl) = 10*log10(handles.parentHandles.powerS (ch, tr, 1:handles.parentHandles.pwl));
    end    
    % fit of the power expressed in Db
    fitTerms = fitF\fitS'; 
    % compute the linear regression
    PWSlinearFit = fitTerms(2)*log10(handles.parentHandles.pwsf) + fitTerms(1);
    handles.whitenedPS = spectra - PWSlinearFit';
end

handles.whiteFlag = 1;
zebraguidata(hObject, handles);

function whitenedSpectra_ButtonDownFcn(hObject, eventdata, handles)
% persistent leftF rightF

% Computation of the peak frequency in a given interval
% [x,y,button] = ginput
% C = get (gca, 'CurrentPoint');
% txOut = sprintf('Frequency (Hz): %0.1f  Power (Db): %0.2f', C(1,1),C(1,2));
% set (handles.infoString,'String',txOut);
%if (handles.whWinLeft = 0),
    % first call. Define first edge of the interval
    
% % --- Executes on mouse motion over figure - except title and menu.
% function ZebraSpectre_WindowButtonMotionFcn(hObject, eventdata, handles)
% % % hObject    handle to ZebraSpectre (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    structure with handles and user data (see GUIDATA)
% C = get (gca, 'CurrentPoint');
% txOut = sprintf('Frequency (Hz): %0.1f  Power (Db): %0.2f', C(1,1),C(1,2));
% set (handles.infoString,'String',txOut);
% % sel = get(gcf,'SelectionType');
% % switch sel
% %     case 'normal'
% %     case 'open'
% %         if (handles.whWinLeft == 0),
% %             % first time we enter here
% %             sel
%             clock
% end

function freqLog_Callback(hObject, eventdata, handles)

function sendToClipboard_Callback(hObject, eventdata, handles)

% send to the clipboard 

% set(gcf,'units','pixels');
% gcfPosition = get(gcf,'position');      % get the position of the Zebra Spectre figure expressed in pixels
% %imageData = screencapture(handles.inset_WB, [], 'clipboard');
% imageData = screencapture(gcf, gcfPosition, 'clipboard');
hgexport(gcf,'-clipboard');

% --- Executes on mouse motion over figure - except title and menu.
function ZebraSpectre_WindowButtonMotionFcn(hObject, eventdata, handles)
C = get (gca, 'CurrentPoint');
txOut = sprintf('Frequency (Hz): %0.1f  Power (Db): %0.2f', C(1,1),C(1,2));
set (handles.infoString,'String',txOut);

function  plotPowerSpectra (hObject, handles)

% This function has been moved to ZebraSpectre on December 2017
% hObject and handles refer to ZebraSpectre.

% First thing first: define the page title
handles.parentHandles = zebraguidata(handles.parentObject);
dirs = regexp(handles.parentHandles.dir_in,'\','split');
l = size(dirs);
set(handles.spectreTitle,'String',dirs(l(1,2)-1));

tmp = [];
axes(handles.powerSpectra);
hold on;
lim = axis;
cla

% for the time being, I use the flags handled by the parent process
% (ZebraExplore)

flagMean = get(handles.parentHandles.plotMeanFlag,'Value');
flagMP = get(handles.parentHandles.plotMultiCh,'Value');
flagTr =  get(handles.parentHandles.plotTrCk,'Value');
flagN = get(handles.normSpectraCk,'Value');
nch = handles.parentHandles.nCh;
pwl = handles.parentHandles.pwl;
ctr = handles.parentHandles.currentTrial;
cch = handles.parentHandles.currentCh;
if flagMP
    if flagMean
        tmp(1:nch,1:pwl) = handles.parentHandles.meanPWS (1:nch,1:pwl);
        % plot the mean traces
        gca.ColorOrderIndex = 1;
        if flagN
            % compute etc must be moved here from ZebraExplore
            norm = computeSpectraNormalisationFactors (handles,tmp);
        else
            norm(1:nch) = 1;
        end
        for i=1:nch
            tmp(i,:) = norm(i) * tmp(i,:);
        end    
        plot(handles.parentHandles.pwsf,10*log10(tmp.'),'LineWidth',2);            % spectra in decibels
    end
    if flagTr
        % plot the current trace
        tmp(1:nch,1:pwl) = handles.parentHandles.powerS (1:nch,ctr,1:pwl);
        % plot the mean traces
        gca.ColorOrderIndex = 1;
        if flagN
            norm = computeSpectraNormalisationFactors (handles,tmp);
        else
            norm(1:nch) = 1;
        end
        for i=1:nch
            tmp(i,:) = norm(i) * tmp(i,:);
        end    
        plot(handles.parentHandles.pwsf,10*log10(tmp.'));            % spectra in decibels
    end
else
    if flagMean
        tmp(1:pwl) = handles.parentHandles.meanPWS (cch,1:pwl);        
        plot(handles.parentHandles.pwsf,10*log10(tmp.'),'LineWidth',2);            % spectra in decibels
    end
    if flagTr
        tmp(1:pwl) = handles.parentHandles.powerS (cch,ctr,1:pwl);
        plot(handles.parentHandles.pwsf,10*log10(tmp.'));            % spectra in decibels
    end
end
lim(1) = str2double(get(handles.minFreq,'String'));
lim(2) = str2double(get(handles.maxFreq,'String'));
lim(3) = str2double(get(handles.minPower,'String'));
lim(4) = str2double(get(handles.maxPower,'String'));
axis (lim);

if get(handles.freqLog,'Value')
    set(handles.powerSpectra,'XScale','log');   
else
    set(handles.powerSpectra,'XScale','lin');       
end
ylabel('Power (dB)');

% Now plot the whitened spectra (if computed...)

if handles.whiteFlag,
    axes(handles.whitenedSpectra);
    plot(handles.parentHandles.pwsf,handles.whitenedPS,'LineWidth',2);
    minF = str2double(get (handles.minFreq, 'String'));
    maxF = str2double(get (handles.maxFreq, 'String'));
    maxP = str2double(get (handles.maxPower, 'String'));
    minP = str2double(get (handles.minPower, 'String'));
    if get(handles.whitenByFreq,'Value'),
        axis([minF maxF lim(3) lim(4)]);
    else
        axis([minF maxF -5 5]);
    end
    if get(handles.freqLog,'Value')
        set(handles.whitenedSpectra,'XScale','log');   
    else
        set(handles.whitenedSpectra,'XScale','lin');   
    end
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Delta power (Db)');
end

function whitenByFreq_Callback(hObject, eventdata, handles)

function fromTtxt_Callback(hObject, eventdata, handles)

function fromTtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function toTtxt_Callback(hObject, eventdata, handles)

function toTtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function computeSpectra_Callback(hObject, eventdata, handles)
% Here we compute the power spectra after reading all the options specfied
% in the GUI. This function is called by ZebraExplore upon loading of a new
% data and it is evoked locally within the ZebraSpectre GUI everytime that
% we change the computation paramenters

% One choice should be made: where are we going to store the spectra? In
% the parentHandles structure or in the local handles? Or both?

% remember: parentHandles point to ZebraExplore

% The first thing to do is to align the parentHandles structure with its current content
handles.parentHandles = zebraguidata(handles.parentObject);

handles.parentHandles.powerS = [];
handles.parentHandles.powerS2 = [];
if (handles.parentHandles.fileInNum)
    tmp = [];
    tmp2 = [];
    if get(handles.windowCk,'Value')
        % read the temporal limits of spectra computation converting from sec to points
        tin = 1+handles.parentHandles.acqF * str2double(get(handles.fromTtxt,'string'));     
        tfin = handles.parentHandles.acqF * str2double(get(handles.toTtxt,'string'));
        if handles.window2Ck.Value %gab 2018/11/04: second window
            tin2 = 1+handles.parentHandles.acqF * str2double(get(handles.fromT2txt,'string'));     
            tfin2 = handles.parentHandles.acqF * str2double(get(handles.toT2txt,'string'));
        end
    else
        tin = 1;
        tfin = handles.parentHandles.dtaLen;
    end
    for i = 1:handles.parentHandles.nCh
        for j= 1:handles.parentHandles.nTrials
            tmp(1:(tfin-tin)+1) = handles.parentHandles.workLFP(i,j,tin:tfin);
            if handles.window2Ck.Value %gab 2018/11/04: second window
                tmp2(1:(tfin-tin)+1) = handles.parentHandles.workLFP(i,j,tin2:tfin2);
            end
            if handles.subMean_ck.Value
                tmp = tmp - mean(tmp);      %GAB: subtract the mean
                if handles.window2Ck.Value %gab 2018/11/04: second window
                    tmp2 = tmp2 - mean(tmp2);
                end 
            end
            flagChronux = get(handles.parentHandles.ChronuxCk,'Value');
            if flagChronux
                params.tapers = [5 9] ;%[5 9];
                params.Fs = handles.parentHandles.acqF;
                params.fpass = [0 1000];
%                zebraguidata(hObject, handles); % why here? what is good for?
                [pxx,f] = mtspectrumc(tmp, params);
                f = f';
                if handles.window2Ck.Value %gab 2018/11/04: second window
                    [pxx2,f2] = mtspectrumc(tmp2, params);
                    f2 = f2';
                end
            else
                [pxx,f] = pwelch(tmp,[],[],[],handles.parentHandles.acqF);
                if handles.window2Ck.Value %gab 2018/11/04: second window
                    [pxx2,f2] = pwelch(tmp2,[],[],[],handles.parentHandles.acqF);
                end
            end
            l = size(f);
            handles.parentHandles.powerS (i,j,1:l) = pxx(1:l) ;
            if handles.window2Ck.Value %gab 2018/11/04: second window
                l2 = size(f2);
                handles.parentHandles.powerS2 (i,j,1:l) = pxx2(1:l2) ;
            end
        end
    end
    tmp = [];
    handles.parentHandles.pwsf = f;
    handles.parentHandles.pwl = l;
    if handles.window2Ck.Value %gab 2018/11/04: second window
        handles.parentHandles.pwsf2 = f2;
        handles.parentHandles.pwl2 = l2;
    end
end

handles.parentHandles.meanPWS = [];
% If ntrial > 1 than compute the mean over trials.
if get(handles.parentHandles.skip1trialCk,'value'), trialBegins = 2;
else
    trialBegins = 1;
end
if ((handles.parentHandles.nTrials-trialBegins)>1)
    for (i=1:handles.parentHandles.nCh)
        % compute mean power spectra
        tmp = [];
        tmp (1:handles.parentHandles.nTrials-trialBegins+1,1:handles.parentHandles.pwl) =...
            handles.parentHandles.powerS(i,trialBegins:handles.parentHandles.nTrials,1:handles.parentHandles.pwl);
        handles.parentHandles.meanPWS (i,1:handles.parentHandles.pwl) = mean (tmp,1);
        if handles.window2Ck.Value %gab 2018/11/04: second window
            tmp2 = [];
            tmp2(1:handles.parentHandles.nTrials-trialBegins+1,1:handles.parentHandles.pwl2) =...
                handles.parentHandles.powerS2(i,trialBegins:handles.parentHandles.nTrials,1:handles.parentHandles.pwl2);
            handles.parentHandles.meanPWS2 (i,1:handles.parentHandles.pwl2) = mean (tmp2,1);
        end
    end
end

zebraguidata(handles.parentObject, handles.parentHandles);
% zebraguidata(gcbo, handles);
% compute whitened spectra, WTF
whiteOut_Callback(hObject, [], handles);
handles = zebraguidata(hObject);
plotPowerSpectra (hObject, handles);        % finally, plot the spectra
disp('Spectra computed!')

function normF = computeSpectraNormalisationFactors (handles,tmp)
% first: compute the indexes that define the normalisation interval
lowF = str2double(get(handles.normLowFreq,'String'));
highF = str2double(get(handles.normHighFreq,'String'));
i = 1;
f = handles.handlesParent.pwsf(1); 
while f < lowF
    i = i+1;
    f = handles.handlesParent.pwsf(i);         
end
iStart = i;
while f < highF
    i = i+1;
    f = handles.handles.Parent.pwsf(i);
end
iEnd = i-1;
% second: compute the mean value of spectral power for channel 1
%        tmp(1:handles.nCh,1:handles.pwl)
%norm = sum(handles.powerS(1,handles.currentTrial,iStart:iEnd));
normF(1) = sum(tmp(1,iStart:iEnd));
% compute the integral of the others channels and compute the
% array of the normalization factors.
for (j=2:handles.handlesParent.nCh)
    normF(j) = normF(1)/sum(tmp (j,iStart:iEnd));
end
normF(1) = 1;

function windowCk_Callback(hObject, eventdata, handles)
if hObject.Value
    handles.window2Ck.Enable = 'on';
else
    handles.window2Ck.Enable = 'off';
end

function subMean_ck_Callback(hObject, eventdata, handles)


function normSpectraCk_Callback(hObject, eventdata, handles)



function fromT2txt_Callback(hObject, eventdata, handles)
% hObject    handle to fromT2txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fromT2txt as text
%        str2double(get(hObject,'String')) returns contents of fromT2txt as a double


% --- Executes during object creation, after setting all properties.
function fromT2txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fromT2txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function toT2txt_Callback(hObject, eventdata, handles)
% hObject    handle to toT2txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of toT2txt as text
%        str2double(get(hObject,'String')) returns contents of toT2txt as a double


% --- Executes during object creation, after setting all properties.
function toT2txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toT2txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in window2Ck.
function window2Ck_Callback(hObject, eventdata, handles)
% hObject    handle to window2Ck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of window2Ck
