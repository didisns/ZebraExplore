function plotSpectrogram (hObject, handles)
% This function performs the plot in the SPG window.
% Depending on the selected radio button different data are plotted down
% there.

tmp = [];
axes(handles.spectraWide);      % select the SPG axes
cla
% load in all necessary parameters from the handles structure
freqLow = str2double(get(handles.SPGfromFreq,'String'));
freqHigh = str2double(get(handles.SPGtoFreq,'String'));
handles.SPGlowDB = str2double(get(handles.lowerDB,'String'));
handles.SPGhighDB = str2double(get(handles.upperDB,'String'));
flagMd = get(handles.plotMeanFlag,'Value');     % plot average trace
flagTr = get(handles.plotTrCk,'Value');         % plot single trace
flagFlt = get(handles.filterBPpower,'Value');
flagBPpower = get(handles.gammaPower,'Value');
flagEI = get(handles.EIbutton,'Value');
flagSPG = get(handles.SPG,'Value');
flagBP = get(handles.BP,'Value');
flagLog = get(handles.SPGlogCk,'Value');
flagLeak = get (handles.displayLeakCk,'Value');

flagSkip1st = get(handles.skip1trialCk,'Value'); %gab&enrico 2018/11/02 provvisorio

if flagSPG
    if flagMd
        if (handles.deleakedFlag && flagLeak)
            temp(1:handles.freqN,1:handles.spgl) = handles.meanSpgDeleaked(handles.currentCh,1:handles.freqN,1:handles.spgl);
        else
            temp(1:handles.freqN,1:handles.spgl) = handles.meanSpg(handles.currentCh,1:handles.freqN,1:handles.spgl);    
        end
    else    
        if (handles.deleakedFlag && flagLeak)
            temp(1:handles.freqN,1:handles.spgl) = handles.spgDeleaked(handles.currentCh,handles.currentTrial,1:handles.freqN,1:handles.spgl);
        else
            temp(1:handles.freqN,1:handles.spgl) = handles.spg(handles.currentCh,handles.currentTrial,1:handles.freqN,1:handles.spgl);
        end
    end
    % the spectrogram is represented in decibels
    imagesc(handles.spgt,-handles.spgw,10*log10(temp))
    axis ([handles.tmin handles.tmax -freqHigh -freqLow]);
    caxis ([handles.SPGlowDB handles.SPGhighDB]);
end

if flagBP
    axis xy;
    hold all;
    % collect the BP data to plot
    flagBP1 = get(handles.bpCheck1,'Value');
    flagBP2 = get(handles.bpCheck2,'Value');
    flagBP3 = get(handles.bpCheck3,'Value');
    flagHP = get(handles.hpCheck,'Value');
    
    % BP power data are plotted one at the time 
    if flagBP1
        if flagMd
            % plot average of BP power
%             tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
              tmpX(1:handles.BPpowerL(2)) = handles.BPpower (1,2,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(2));
              tmpY(1:1:handles.BPpowerL(2)) = ...
            squeeze(mean(handles.BPpower (2,2,handles.currentCh,flagSkip1st+1:end,1:handles.BPpowerL(2)),4));
            plot(tmpX,tmpY,'Color','c','LineWidth',1);
        end
        if flagTr
            % plot single trial
            tmpX(1:handles.BPpowerL(2)) = ...
                handles.BPpower (1,2,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(2));  
            tmpY(1:1:handles.BPpowerL(2)) = ...
                handles.BPpower (2,2,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(2));
            plot(tmpX,tmpY,'Color','c','LineWidth',1);
        end
    end
    if flagBP2
        if flagMd
            % plot average of BP power
    %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
            tmpX(1:handles.BPpowerL(3)) = handles.BPpower (1,3,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(3));
            tmpY(1:1:handles.BPpowerL(3)) = ...
                squeeze(mean(handles.BPpower (2,3,handles.currentCh,flagSkip1st+1:end,1:handles.BPpowerL(3)),4));
            plot(tmpX,tmpY,'Color','r','LineWidth',1);
        end
        if flagTr
            % plot single trial
            tmpX(1:handles.BPpowerL(3)) = ...
                handles.BPpower (1,3,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(3));  
            tmpY(1:1:handles.BPpowerL(3)) = ...
                handles.BPpower (2,3,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(3));
            plot(tmpX,tmpY,'Color','r','LineWidth',1);
        end
    end
    if flagBP3
        if flagMd
            % plot average of BP power
    %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
            tmpX(1:handles.BPpowerL(4)) = handles.BPpower (1,4,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(4));
            tmpY(1:1:handles.BPpowerL(4)) = ...
                squeeze(mean(handles.BPpower (2,4,handles.currentCh,flagSkip1st+1:end,1:handles.BPpowerL(4)),4));
            plot(tmpX,tmpY,'Color','g','LineWidth',1);
        end
        if flagTr
            % plot single trial
            tmpX(1:handles.BPpowerL(4)) = ...
                handles.BPpower (1,4,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(4));  
            tmpY(1:1:handles.BPpowerL(4)) = ...
                handles.BPpower (2,4,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(4));
            plot(tmpX,tmpY,'Color','g','LineWidth',1);
        end
    end
    if flagHP
        if flagMd
            % plot average of BP power
    %        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
            tmpX(1:handles.BPpowerL(5)) = handles.BPpower (1,5,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(5));
              tmpY(1:1:handles.BPpowerL(5)) = ...
            squeeze(mean(handles.BPpower (2,5,handles.currentCh,flagSkip1st+1:end,1:handles.BPpowerL(5)),4));
            plot(tmpX,tmpY,'Color','blue','LineWidth',1);
        end
        if flagTr
            % plot single trial
            tmpX(1:handles.BPpowerL(5)) = ...
                handles.BPpower (1,5,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(5));  
            tmpY(1:1:handles.BPpowerL(5)) = ...
                handles.BPpower (2,5,handles.currentCh,handles.currentTrial,1:handles.BPpowerL(5));
            plot(tmpX,tmpY,'Color','blue','LineWidth',1);
        end
    end
%         if flagFlt
%             % adapt the filter frame and order to the length of the data
%             frame = int32(handles.spgl / 64) * 2 + 1;   % frame must be odd. Thats is why *2+1 !
%             if (frame > 3), tmp = sgolayfilt(tmp,3,frame);
%             else
%                 if (frame > 2), tmp = sgolayfilt(tmp,2,frame);
%                 end
%             end
%         end

    axis ([handles.tmin handles.tmax 0 inf]);
    if flagLog, set(gca,'yscale','log');
    end
    hold off    
end
    
if flagBPpower
    if flagMd
        tmp(1:handles.spgl) = handles.meanSpgPlot(handles.currentCh,1:handles.spgl);
    else
        if (handles.deleakedFlag && flagLeak)
            tmp(1:handles.spgl) = handles.spgPlotDeleaked(handles.currentCh,handles.currentTrial,1:handles.spgl);
        else
            tmp(1:handles.spgl) = handles.spgPlot(handles.currentCh,handles.currentTrial,1:handles.spgl);
        end
        if flagFlt
            % adapt the filter frame and order to the length of the data
            frame = int32(handles.spgl / 64) * 2 + 1;   % frame must be odd. Thats is why *2+1 !
            if (frame > 3), tmp = sgolayfilt(tmp,3,frame);
            else
                if (frame > 2), tmp = sgolayfilt(tmp,2,frame);
                end
            end
        end
    end    
    plot (handles.spgt,tmp);
    axis ([handles.tmin handles.tmax -inf inf]);
    if flagLog, set(gca,'yscale','log');
        end
end

if flagEI
    if flagMd
        tmp(1:handles.spgl) = handles.meanEIrat(handles.currentCh,1:handles.spgl);
    else
        tmp(1:handles.spgl) = handles.EIrat(handles.currentCh,handles.currentTrial,1:handles.spgl);
        if flagFlt
            % adapt the filter frame and order to the length of the data
            frame = int32(handles.spgl / 64) * 2 + 1;   % frame must be odd. Thats is why *2+1 !
            if (frame > 3), tmp = sgolayfilt(tmp,3,frame);
            else
                if (frame > 2), tmp = sgolayfilt(tmp,2,frame);
                end
            end
        end
    end    
    plot (handles.spgt,tmp);
    axis ([handles.tmin handles.tmax -inf inf]);
    if flagLog, set(gca,'yscale','log');
    end
    %set(gca,'yscale','log');
end