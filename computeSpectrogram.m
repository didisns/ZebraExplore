function computeSpectrogram (hObject, handles)

% Hello, my name is me and this is March 3, 2017
% This code computes the spectograms of all data.
% If so desired, the SPG is computed after removal of the (large) low
% frequency transients in order to attenuate the associated spectral
% leakage.

SPGfrom = str2double(get(handles.SPGfromFreq,'String'));
SPGto = str2double(get(handles.SPGtoFreq,'String'));
freqStep = str2double(get(handles.SPGresolution,'String'));
freq = SPGfrom:freqStep:SPGto;
handles.freqN = (SPGto-SPGfrom)/freqStep + 1;
window = str2double(get(handles.SPGwindow,'String'));
noverlap = str2double(get(handles.SPGoverlap,'String'));
leakageFlag = get (handles.leakageCk,'Value');
handles.spg=[];
for i = 1:handles.nCh
    for j= 1:handles.nTrials
        flag = get(handles.decimateFlag, 'Value');
        if flag
            % decimate the input data in case of really long data set or
            % when using a slow pc...
            LFPsmall = decimate(handles.workLFP(i,j,1:handles.dtaLen),10);
            LFPsmall(i,j,1:handles.dtaLen) = decimate(handles.workLFP(i,j,1:handles.dtaLen),10);
            [s, w, t, ps] = spectrogram(LFPsmall,hamming(window),noverlap, freq, handles.acqF/10, 'yaxis');        
        else
            tmp = handles.workLFP(i,j,1:handles.dtaLen);
            flagChronux = get(handles.ChronuxCk,'Value');
            if flagChronux
                params.tapers = [3 5] ;%[10 0.5 1];
                params.Fs = handles.acqF;
                params.fpass = [SPGfrom SPGto];
                [rollingWin, overlap] = setSPGwindow (handles, handles.dtaLen*handles.sp);
                [ps, t, w] = mtspecgramc(tmp,[rollingWin overlap],params);
                l = size (ps);
                handles.spgl = l(1);
                handles.freqN = l(2);
                ps = ps';
                handles.spg(i,j,1:handles.freqN,1:handles.spgl) = ps (1:handles.freqN,1:handles.spgl);
                if leakageFlag
                    % first: subtract the LP filtered data from the workLFP
                    % data. Second: compute its SPG
                    tmp = squeeze(handles.workLFP(i,j,:))-squeeze(handles.bandPassed_LFP(6,i,j,:));
                    [ps, t, w] = mtspecgramc(tmp,[rollingWin overlap],params);
                    ps = ps';
                    handles.spgDeleaked(i,j,1:handles.freqN,1:handles.spgl) = ps (1:handles.freqN,1:handles.spgl);
                    handles.deleakedFlag = 1;
                end
            else
                [s, w, t, ps] = spectrogram(tmp,hamming(window),noverlap, freq, handles.acqF, 'yaxis');
                l = size (ps);
                handles.spgl=l(2);
                handles.spg(i,j,1:handles.freqN,1:handles.spgl) = ps (1:handles.freqN,1:handles.spgl);
                if leakageFlag
                    tmp = squeeze(handles.workLFP(i,j,:))-squeeze(handles.bandPassed_LFP(6,i,j,:));
                    [s, w, t, ps] = spectrogram(tmp,hamming(window),noverlap, freq, handles.acqF, 'yaxis');
                    handles.spgDeleaked(i,j,1:handles.freqN,1:handles.spgl) = ps (1:handles.freqN,1:handles.spgl);
                    handles.deleakedFlag = 1;
                end
            end
        end
        % compute the short time E/I index
        % First: create the array containing the indexes of the spectral density in the interval 20-50 Hz
        [segmentI, segmentF, freqCnt] = findSPGelements (w, handles.freqN, 20, 50);
        % Second: loop on all time point of the spg to compute the linear
        % fit to the spectra segment
        for ti=1:handles.spgl
            spectraSegment = handles.spg(i,j,segmentI(1):segmentI(freqCnt),ti);
            % compute now the linear fit to the log log segment
            segmentF = squeeze(segmentF);   % remove the singleton dimensions
            spectraSegment = squeeze(spectraSegment);
            % compute slope and offset of the fit
            %X = [ones(length(segmentF),1) segmentF'];
            X = [ones(length(segmentF),1) log(segmentF)'];
            fit = X\log(spectraSegment);
            handles.EIrat (i,j,ti) = fit(2);
            % the linear regression operators '\' needs row vectors.
        end    
        % compute the integral of the SPG in 'gamma' band
        PWfrom = str2double(get(handles.envSPGfrom,'String'));
        PWto = str2double(get(handles.envSPGto,'String'));        
        [segmentI, segmentF, freqCnt] = findSPGelements (w, handles.freqN, PWfrom, PWto);
        handles.spgPlot(i,j,1:handles.spgl) = sum(handles.spg(i,j,segmentI(1):segmentI(freqCnt),1:handles.spgl));
        if leakageFlag
            handles.spgPlotDeleaked(i,j,1:handles.spgl) = sum(handles.spgDeleaked(i,j,segmentI(1):segmentI(freqCnt),1:handles.spgl));            
        end    
    end    
end
handles.spgt = t;
handles.spgw = w;
%zebraguidata(gcbo, handles);
zebraguidata(hObject, handles);

function [segmentI, segmentF, freqCnt] = findSPGelements (w, freqNum, lowFreq, highFreq)
% Given a frequency range (lowFreq and highFreq) identify the elements of
% the spectrogram that fall within this range. The function returns an array
% (freqIndex) containing the indexes of the elements, and an array
% containing the relative frequencies.

segmentI = [];
segmentF = [];

freqCnt = 0;
for fj=1:freqNum
    if w(fj)>highFreq      % interupt the iterations if we pass the frequency limit
        break
    end
    if w(fj)>=lowFreq
        freqCnt = freqCnt + 1;
        segmentI(freqCnt) = fj;
        segmentF(freqCnt) = w(fj);
    end
end   

function [wnd ovl] = setSPGwindow (handles, deltaT)
% This function reads (from GUI) or generates the values of window size and
% overlap to be used for the spectrogram computation.
% The automatic computation uses the width of the data to adapt the width
% of the rolling window.

flagAuto = get(handles.ChronuxAutoConfigCk,'Value');
if flagAuto
    if deltaT<5, wnd = 0.1;
    else
        if deltaT<10, wnd = 0.2;
        else
            if deltaT<30, wnd = 0.2;
            else
                if deltaT<100, wnd = 0.5;
                else
                    wnd = 0.5;
                end
            end
        end
    end
    ovl = wnd/5;
else
    wnd = str2double(get(handles.ChronuxRW,'String'));
    ovl = str2double(get(handles.ChronuxOvl,'String'));
end    
