
function doInsets (hObject, handles)
% Dec 26, 2016. Modified for the new data structure that include Ch and Trial numbers.
% March 06, 2017. Modified for the dual display of normal and
% leak-corrected SPG
% May 20, 2017. Modified to plot the temporal profile of band-limited
% spectral power.

% restore the visibility of all clipboard copy buttons
set (handles.sc4,'visible','on')

if (handles.xnow>0)
    % now compute the x limits of the inset plot
    % beginning and end of the magnified window
    x_from = handles.xnow - handles.insetW/2;
    x_to = handles.xnow + handles.insetW/2;       % this is in time units, convert now in indexes of the arrays
    i_from = int32(x_from / handles.sp +1);
    i_to = int32(x_to / handles.sp +1);           % indexes in sampling periods
    np = i_to - i_from+1;
    x = x_from:handles.sp:x_to;
    %dimx = size(x)

    axes(handles.inset_WB);
    hold off;
    i_to = i_from + length(x) - 1;
    tmp(i_from:1:i_to) = handles.LFP(handles.currentCh,handles.currentTrial,i_from:1:i_to);
    pl = plot(x,tmp(i_from:1:i_to));
    axis([x_from x_to -inf inf]);
    hold off;
    flagDisplayL = get (handles.displayLeakCk,'Value');
    
    plotMultiBandInset(handles)
    
    % Compute inset spectrogram only if width > 1 s
    if (handles.insetW>1)
        SPGfrom = str2double(get(handles.SPGfromFreq,'String'));
        SPGto = str2double(get(handles.SPGtoFreq,'String'));
        freqStep = str2double(get(handles.SPGresolution,'String'));
        if (flagDisplayL && handles.deleakedFlag)
            segment = squeeze(handles.workLFP(handles.currentCh,handles.currentTrial,i_from:1:i_to));
            segment = segment - squeeze(handles.bandPassed_LFP(6,handles.currentCh,handles.currentTrial,i_from:1:i_to));
        else    
            segment = handles.workLFP(handles.currentCh,handles.currentTrial,i_from:1:i_to);
        end

        % First: compute the SPGs of the data segments.
        flagChronux = get(handles.ChronuxCk,'Value');
        if flagChronux
            params.tapers = [3 5] ;%[10 0.5 1];
            params.Fs = handles.acqF;
            params.fpass = [SPGfrom SPGto];
            [rollingWin, overlap] = setSPGwindow (handles, handles.insetW);
            [ps, t, f] = mtspecgramc(segment,[rollingWin overlap],params);
            ps = ps';
        else
            freq = [SPGfrom:freqStep:SPGto];
            handles.freqN = (SPGto-SPGfrom)/freqStep+1;
            % window is expressed in points
            window = str2double(get(handles.SPGwindow,'String'));
            noverlap = str2double(get(handles.SPGoverlap,'String'));
            [s, f, t, ps]=spectrogram(segment,hamming(window),noverlap, freq, handles.acqF);        
        end
        sz = size(t);
        sz = sz (2);
        % Second, compute and plot the desired output
        flagBPpw = get(handles.gammaPower,'Value');
        flagEI = get(handles.EIbutton,'Value');
        flagSPG = get(handles.SPG,'Value');
        flagBP = get(handles.BP,'Value');
        flagLog = get(handles.SPGlogCk,'Value');
        flagMd = get(handles.plotMeanFlag,'Value');     % plot average trace
        flagTr = get(handles.plotTrCk,'Value');         % plot single trace

        axes(handles.spectraInset);

        if (flagBPpw)
            PWfrom = str2double(get(handles.envSPGfrom,'String'));
            PWto = str2double(get(handles.envSPGto,'String'));        
            [segmentI, segmentF, freqCnt] = findSPGelements (f, handles.freqN, PWfrom, PWto);
            y = sum(ps(segmentI(1):segmentI(freqCnt),:));
            plot(t,y);
            axis([t(1,1) t(1,sz) -inf inf]);
            if flagLog, set(gca,'yscale','log');
            end    
        end

        if flagSPG
            % the spectrogram is represented in decibels
            imagesc(t,-f,10*log10(ps))
            set(get(gca,'YLabel'),'String','')
            set(gca,'YScale','lin');
            axis ([t(1,1) t(1,sz) -SPGto -SPGfrom]);
            caxis ([handles.SPGlowDB handles.SPGhighDB]);
        end
    
        if flagBP
            axis xy;    % this function set the polarity of the plot. Increasing Y -> top
            hold all;   % we need this to plot multiple data
            % collect the BP data to plot
            flagBP1 = get(handles.bpCheck1,'Value');
            flagBP2 = get(handles.bpCheck2,'Value');
            flagBP3 = get(handles.bpCheck3,'Value');
            flagHP = get(handles.hpCheck,'Value');

            % BP power data are plotted one at the time 
            if flagBP1
                % find the x axis limits
                i_from = floor((x_from-5*handles.BPpowerSP(2)) / handles.BPpowerSP(2));
                i_to = ceil((x_to-5*handles.BPpowerSP(2)) / handles.BPpowerSP(2));
                tmpX = [];
                tmpY = [];
                if flagMd
                    % plot average of BP power
            %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
                end
                if flagTr
                    % plot single trial
                    tmpX(i_from:i_to) = ...
                        handles.BPpower (1,2,handles.currentCh,handles.currentTrial,i_from:i_to);  
                    tmpY(i_from:i_to) = ...
                        handles.BPpower (2,2,handles.currentCh,handles.currentTrial,i_from:i_to);
                    plot(tmpX,tmpY,'Color','c','LineWidth',1);
                end
            end
            if flagBP2
                i_from = floor((x_from-5*handles.BPpowerSP(3)) / handles.BPpowerSP(3));
                i_to = floor((x_to-5*handles.BPpowerSP(3)) / handles.BPpowerSP(3));
                tmpX = [];
                tmpY = [];
                if flagMd
%                    % plot average of BP power
            %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
                end
                if flagTr
                    % plot single trial
                    tmpX(i_from:i_to) = ...
                        handles.BPpower (1,3,handles.currentCh,handles.currentTrial,i_from:i_to);  
                    tmpY(i_from:i_to) = ...
                        handles.BPpower (2,3,handles.currentCh,handles.currentTrial,i_from:i_to);
                    plot(tmpX,tmpY,'Color','r','LineWidth',1);
                end
            end
            if flagBP3
                i_from = floor((x_from-5*handles.BPpowerSP(4)) / handles.BPpowerSP(4));
                i_to = ceil((x_to-5*handles.BPpowerSP(4)) / handles.BPpowerSP(4));
                tmpX = [];
                tmpY = [];
                if flagMd
 %                   % plot average of BP power
            %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
                end
                if flagTr
                    % plot single trial
                    tmpX(i_from:i_to) = ...
                        handles.BPpower (1,4,handles.currentCh,handles.currentTrial,i_from:i_to);  
                    tmpY(i_from:i_to) = ...
                        handles.BPpower (2,4,handles.currentCh,handles.currentTrial,i_from:i_to);
                    plot(tmpX,tmpY,'Color','g','LineWidth',1);
                end
            end
            if flagHP
                i_from = floor((x_from-5*handles.BPpowerSP(5)) / handles.BPpowerSP(5));
                i_to = ceil((x_to-5*handles.BPpowerSP(5)) / handles.BPpowerSP(5));
                tmpX = [];
                tmpY = [];
                if flagMd
                    % plot average of BP power
            %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
                end
                if flagTr
                    % plot single trial
                    tmpX(i_from:i_to) = ...
                        handles.BPpower (1,5,handles.currentCh,handles.currentTrial,i_from:i_to);  
                    tmpY(i_from:i_to) = ...
                        handles.BPpower (2,5,handles.currentCh,handles.currentTrial,i_from:i_to);
                    plot(tmpX,tmpY,'Color','blue','LineWidth',1);
                end
            end

            axis([x_from x_to -inf inf]);
            if flagLog, set(gca,'yscale','log');
            end
            hold off    
        end

        if flagEI
            [segmentI, segmentF, freqCnt] = findSPGelements (f, handles.freqN, 20, 50);
            segmentF = squeeze(segmentF);   % remove the singleton dimensions
            for ti=1:sz
                spectraSegment = ps(segmentI(1):segmentI(freqCnt),ti);
                % compute now the linear fit to the log log segment                
                spectraSegment = squeeze(spectraSegment);
                % the linear fit must allow for offset <> 0!!!
                y (ti) = log(segmentF)'\log(spectraSegment);
            end
            plot(t,y);
            axis([t(1,1) t(1,sz) -inf inf]);
            % the linear regression operators '\' needs row vectors.
            if flagLog, set(gca,'yscale','log');
            end    
        end
    end
end
zebraguidata(hObject, handles);

function plotMultiBandInset(handles)
% Here we plot the multi band inset window according to the selected bands
% Dec 26, 2016. Modified for the new data structure.

x_from = handles.xnow - handles.insetW/2;
x_to = handles.xnow + handles.insetW/2;       % this is in time units, convert now in indexes of the arrays
i_from = floor(x_from / handles.sp +1);
i_to = floor(x_to / handles.sp +1);
np = i_to - i_from+1;
x = x_from:handles.sp:x_to;

axes(handles.inset_HP);
axis([x_from x_to -inf inf]);
cla
set(gca,'XTickLabel',[])
hold all;
% make sure x and the tmp segment have the same lenght
i_to = i_from + length(x) - 1;
flag = get(handles.lpcheck,'Value');
if flag, tmp(i_from:1:i_to) = handles.bandPassed_LFP(1,handles.currentCh, handles.currentTrial,i_from:1:i_to);
    in1 = plot(x,tmp(i_from:1:i_to));
    set(in1,'Color','black');
end
flag = get(handles.bpCheck1,'Value');
if flag, tmp(i_from:1:i_to) = handles.bandPassed_LFP(2,handles.currentCh, handles.currentTrial,i_from:1:i_to);
    in2 = plot(x,tmp(i_from:1:i_to));
    set(in2,'Color','c');
end
flag = get(handles.bpCheck2,'Value');
if flag, tmp(i_from:1:i_to) = handles.bandPassed_LFP(3,handles.currentCh, handles.currentTrial,i_from:1:i_to);
    in3 = plot(x,tmp(i_from:1:i_to));
    set(in3,'Color','r');
end
flag = get(handles.bpCheck3,'Value');
if flag, tmp(i_from:1:i_to) = handles.bandPassed_LFP(4,handles.currentCh, handles.currentTrial,i_from:1:i_to);
    in3 = plot(x,tmp(i_from:1:i_to));
    set(in3,'Color','green');
end
flag = get(handles.hpCheck,'Value');
if flag, tmp(i_from:1:i_to) = handles.bandPassed_LFP(5,handles.currentCh, handles.currentTrial,i_from:1:i_to); 
    in4 = plot(x,tmp(i_from:1:i_to));
    set(in4,'Color','blue');
end
flag = get(handles.plotEnvelopeCk,'Value');
if flag, in5 = plot(x,handles.envelope(i_from:1:i_to));
    set(in5,'Color','black');
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
