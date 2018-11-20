function computeMeanOnTrials (hObject, handles)
% This function is evoched in the presence of multiple trials. 
% Data are averaged on the trial number
% Finally, it computes metrics of the evoched response
% February 5, 2017
% Debug on January, 3 2018.
% Made some check to improve reliability for operation in the MeyerVEP mode
display ('Compute mean on trials');
computeMeanData (hObject, handles)
handles = zebraguidata(hObject);     % refresh handles

skipFlag = get (handles.skip1trialCk,'value');
if skipFlag, trialBegins = 2;
else
    trialBegins = 1;
end

% load from the GUI all the necessary parameters

syCh = str2double(get(handles.synchCh,'String'));
exCh = str2double(get(handles.excludeCh,'String'));
flagExCh = get(handles.excludeChFlag,'Value');
flagSynch = get(handles.internalSynchCk,'Value');
flagTemplate = get(handles.templateCk,'Value');
flagDouble = get(handles.doubleStimCk,'Value');

sweepDuration = handles.sp * (handles.dtaLen-1);    % duration of each trial

% fix the following for meyer veps
%             handles.ntemp = 1;
%             if flagDouble, handles.ntemp = 2;
%             end

halfSweep = int32 (handles.dtaLen/2);       % this can be easily generalize to multiple stim per sweep
halfSweepSPG = int32(handles.spgl/2);

% window for template computation and fit
templFromI = int32(str2double(get(handles.templateFrom,'String'))/handles.sp)+1;
templToI = int32(str2double(get(handles.templateTo,'String'))/handles.sp)+1;
PWfromI = int32(str2double(get(handles.winPWleft,'String'))/handles.sp)+1;  % window for power computation
PWtoI = int32(str2double(get(handles.winPWright,'String'))/handles.sp)+1;
% now check that all thee pointers have legal values, i.e. are compatible
% with the sweep length
if (templToI > handles.dtaLen), templToI = handles.dtaLen;
end
if (PWtoI > handles.dtaLen), PWtoI = handles.dtaLen;
end

templateLen = templToI - templFromI + 1;
tempTempl = zeros(1,templateLen);    
templateConvolutionLen = PWtoI - PWfromI +1;
tempTemplConv = zeros(1,templateConvolutionLen);
handles.templateConv = [];

% windows for baseline computation and response
blFrom = int32(str2double(get(handles.baselineFrom,'String'))/handles.sp)+1;
blTo = int32(str2double(get(handles.baselineTo,'String'))/handles.sp)+1;
respFrom = int32(str2double(get(handles.respFrom,'String'))/handles.sp)+1;
respTo = int32(str2double(get(handles.respTo,'String'))/handles.sp)+1;

% processing of spectrograms
% compute the deltaT of the spectrogram
dlt = (handles.spgt(end)-handles.spgt(1))/(size(handles.spgt,2));
blFromSPG = int32((str2double(get(handles.baselineFrom,'String'))-handles.spgt(1))/dlt) + 1;
if (blFromSPG<1), blFromSPG=1;
end
blToSPG = int32((str2double(get(handles.baselineTo,'String'))-handles.spgt(1))/dlt) + 1;
if (blToSPG<1), blToSPG=1;
end
rsFromSPG= int32((str2double(get(handles.respFrom,'String'))-handles.spgt(1))/dlt) + 1;
if (rsFromSPG<1), rsFromSPG=1;
end
rsToSPG = int32((str2double(get(handles.respTo,'String'))-handles.spgt(1))/dlt) + 1;
bltemp = [];
rstemp = [];

if get(handles.autoclear,'Value')
    handles.respSummary = {};
    handles.EPcnt = 0;
end

%----enrico+gab-2018/10/09
%set sgolay filter window
wSgolay = 2*round(0.5*(0.0101/handles.sp-1))+1;
%----
for (i=1:handles.nCh)
    % Computation of the metrics of the responses
    % Here we exclude the trigger channel and, if it is so marked, the
    % excluded channel.
    
    if (i ~= syCh || ~flagSynch) % if flagSynch is false always do the following
        if ~(flagExCh && i == exCh) 
            % compute evoked response from the mean trials and SPGs
            % first: computation in the time domain
            bltemp = [];
            rstemp = [];
            % compute the template all the way to the end of the
            % convolution window for the computation of the normalized
            % spectral power.
            handles.ntemp = 1;
            if flagDouble, handles.ntemp = 2;
            end
            for ii=1:handles.ntemp
                tFrom(ii) = templFromI + halfSweep * (ii-1);
                tTo(ii) = templToI + halfSweep * (ii-1);
                tempTempl = tempTempl + handles.meanLFP (i,tFrom(ii):tTo(ii)); 
                blFromNow(ii) = blFrom + halfSweep * (ii-1);
                blToNow(ii) = blTo + halfSweep * (ii-1);
                respFromNow(ii) = respFrom + halfSweep * (ii-1);
                respToNow(ii) = respTo + halfSweep * (ii-1);

                tFromPW(ii) = PWfromI + halfSweep * (ii-1);
                tToPW(ii) = PWtoI + halfSweep * (ii-1);
                tempTemplConv = tempTemplConv + handles.meanLFP (i,tFromPW(ii):tToPW(ii));
            end
            
            % now average the mean sweeps
            template = tempTempl / handles.ntemp;
            templateConv = tempTemplConv / handles.ntemp;
            % smooth the template by Savitzky-Golay filtering
            template = sgolayfilt(template,3,wSgolay); %previously framelength was 101
            templateConv = sgolayfilt(templateConv,3,wSgolay);
            
            % align the template to baseline=0
            blDelta = blTo - templFromI;    % number of point in baseline
            offset = mean(template(1:blDelta));     % OK
            template = template - offset;
            templateConv = templateConv - offset;
            % compute the amplitude of the template.
            [templatePeak peak_i] = max(abs(template));      % suggestion: what about normalizing the template to 1?
            handles.templatePeakSigned(i) = template(peak_i);
            templateConv = templateConv / templatePeak;
            if templatePeak ~= max(template), templatePeak = -templatePeak;
            end
            template = template  / templatePeak;
            %handles.templ = template;
            handles.templ(i,:) = template; %gab_2018/09/19 -> make handles.templ a matrix adding the "channel" dimension.
            handles.templateConv(i,:) = templateConv; %gab+enr_2018/10/09
            
            % now the template must be fitted to the corresponding segments
            % of each trial. Each fit returns the offset and the linear
            % scaling factor.
            EPamplitude = [];
            for k=trialBegins:handles.nTrials
                for ii=1:handles.ntemp
                    handles.EPcnt = handles.EPcnt + 1;
                    % first extract the LFP segment to fit with the template
                    segment(1:templateLen) = handles.workLFP(i,k,tFrom(ii):tTo(ii));
                    % Offset the segment to align the baseline to 0
                    offset = mean(segment(1:blDelta));
                    segment = segment - offset;
                    %handles.segment(i,k,ii,:) = segment; %GAB&ENR 2018/10/03
                    %X = [ones(length(segment),1) segment];
                    %fitNow = template'\segment      
                    x = fminsearch(@gixres,1);
                    localEPamp = x ; % after the normalisation the following disappears! * templatePeak;
                    
                    handles.EPamplitude(i,k,ii) = x; % this is used to plot the fitted template
                    handles.EPoffset(i,k,ii) = offset;
                    % computation in the frequency domain by using the
                    % power selected (index: BP4power+1) GAB 06/04/2018
                    bltemp = [];
                    rstemp = [];
                    
                    bltemp = handles.bandPassed_LFP(handles.BP4power+1,i,k,blFromNow(ii):blToNow(ii));
                    bl = squeeze(rms (bltemp));  % do we really need to squeeze?
                 
                    %bltemp = handles.EIrat (i,k,blFromSPG + halfSweepSPG*(ii-1):blToSPG + halfSweepSPG*(ii-1));
                    %bl = squeeze(mean (bltemp,3));
                    rstemp = handles.bandPassed_LFP(handles.BP4power+1,i,k,respFromNow(ii):respToNow(ii));
                    rs = squeeze(rms (rstemp));
                    
                    %rstemp = handles.EIrat (i,k,rsFromSPG + halfSweepSPG*(ii-1):rsToSPG + halfSweepSPG*(ii-1));
                    %rs = squeeze(mean (rstemp,3));
                    FDresp = rs-bl;
                    % extract the band passed data
                    tempPW = handles.bandPassed_LFP(handles.BP4power+1,i,k,tFromPW(ii):tToPW(ii));
                    % multiply it by the fitted template. NO! we should
                    % pass it through the template but not the fitted
                    % template, otherwise the correlation between power and
                    % amplitude is artificially injected in the metric.
                    
                    newTempPW =  templateConv .* squeeze(tempPW)';
                    % and now compute the spectral power
                    templateSPW = rms(newTempPW);
                    
                    % stuff everything in the output table!
                    handles.respSummary(handles.EPcnt,1) = cellstr(handles.file_in);  % 1: file name
                    handles.respSummary(handles.EPcnt,2) = num2cell(i);               % 2: channel number
                    handles.respSummary(handles.EPcnt,3) = num2cell(k);               % 3: trial
                    handles.respSummary(handles.EPcnt,4) = num2cell(ii);              % 4: repeat

                    handles.respSummary(handles.EPcnt,5) = num2cell(localEPamp);      % 5: Peak response
                    handles.respSummary(handles.EPcnt,6) = num2cell(0);               % 6: time to peak
                    handles.respSummary(handles.EPcnt,7) = num2cell(bl);              % 7: Baseline mean BP1 power
                    handles.respSummary(handles.EPcnt,8) = num2cell(rs);              % 8: Response mean BP1 power
                    handles.respSummary(handles.EPcnt,9) = num2cell(FDresp);          % 9: Delta BP1 power
                    handles.respSummary(handles.EPcnt,10) = num2cell(templateSPW);    % 10: Delta BP1 power                
                    handles.respSummary(handles.EPcnt,11) = cellstr(handles.expID);   % 11: Patient name-notes        
                end
            end    
            % compute metrics of the mean responses for each stim repeat
            for ii=1:handles.ntemp
                bltemp = [];
                rstemp = [];
                bltemp = handles.meanLFP (i,blFrom + halfSweep * (ii-1):blTo + halfSweep * (ii-1));
                bl = mean (bltemp,2);
                % the baseline must be subtracted before computing the
                % absolute value
                rstemp = abs(handles.meanLFP (i,respFrom + halfSweep * (ii-1):respTo + halfSweep * (ii-1)) - bl);
                % a bit of filtering before computation of the max
                rstemp = sgolayfilt(rstemp,3,floor(wSgolay/2)+1); %previously framelength
                [resp, imax] = max (rstemp,[],2);
                tpeak = handles.sp*double(respFrom+imax);
                resp = resp * sign(handles.meanLFP (i,respFrom+imax) - bl);
                %resp = handles.meanLFP (i,respFrom+imax) - bl;

                % compute the response in the frequency domain    
                bltemp = [];
                rstemp = [];
                bltemp = handles.meanEIrat (i,blFromSPG + halfSweepSPG*(ii-1):blToSPG + halfSweepSPG*(ii-1));
                bl = mean (bltemp,2);

                rstemp = handles.meanEIrat (i,rsFromSPG + halfSweepSPG*(ii-1):rsToSPG + halfSweepSPG*(ii-1));
                rs = mean (rstemp,2);
                FDresp = rs-bl;

                handles.EPcnt = handles.EPcnt + 1;
                handles.respSummary(handles.EPcnt,1) = cellstr(handles.file_in);  % 1: file name
                handles.respSummary(handles.EPcnt,2) = num2cell(i);               % 2: channel number
                handles.respSummary(handles.EPcnt,3) = num2cell(0);               % 3: trial
                handles.respSummary(handles.EPcnt,4) = num2cell(ii);              % 4: repeat

                handles.respSummary(handles.EPcnt,5) = num2cell(resp);            % 5: Peak response
                handles.respSummary(handles.EPcnt,6) = num2cell(tpeak);           % 6: time to peak
                handles.respSummary(handles.EPcnt,7) = num2cell(bl);              % 7: Baseline mean BP1 power
                handles.respSummary(handles.EPcnt,8) = num2cell(rs);              % 8: Response mean BP1 power
                handles.respSummary(handles.EPcnt,9) = num2cell(FDresp);          % 9: Delta BP1 power
%                handles.respSummary(handles.EPcnt,10) = num2cell(TemplateSPW);    % 10: Delta BP1 power                
                handles.respSummary(handles.EPcnt,11) = cellstr(handles.expID);   % 11: Patient name-notes        
            end
        end        
    end
end

set(handles.zebraSaysData.outEvochedResponses,'Data',handles.respSummary);
zebraguidata(hObject, handles);

    function gixout = gixres(b)
        % this nested function computes the residue of the difference
        % between the template and the given sweep.
        gixout = sum((segment(blDelta:end)-b*template(blDelta:end)).^2);
    end
end


function computeMeanData (hObject, handles)
% Compute means of all the data representations on all trials
% April 25, 2017. A checkbox is added to skip the first trial
% March 2018. Modified to allow conditional computation of BPed data

leakageFlag = get (handles.displayLeakCk,'Value');
skipFlag = get (handles.skip1trialCk,'value');
if skipFlag, trialBegins = 2;
else
    trialBegins = 1;
end

ntemp = 1; 
handles.meanSpgDeleaked = [];
handles.meanSpg = [];
for (i=1:handles.nCh)
    tmp = [];        % initialize tmp
    % create temporary matrix for average computation
    tmp (1:handles.nTrials-trialBegins+1,1:handles.dtaLen) = handles.workLFP(i,trialBegins:handles.nTrials,1:handles.dtaLen);
    handles.meanLFP (i,1:handles.dtaLen) = mean (tmp,1);
    % compute mean of band-passed data
    for (k=1:5)
        if handles.BPcomputed(k)
            tmp (1:handles.nTrials-trialBegins+1,1:handles.dtaLen) = handles.bandPassed_LFP(k,i,trialBegins:handles.nTrials,1:handles.dtaLen);
            handles.meanBP (k,i,1:handles.dtaLen) = mean (tmp,1);
        end
    end
    % compute mean spectrogram
    tmp = [];
    if leakageFlag      %GAB: add plotting of deleaked mean spg.
        tmp (1:handles.nTrials-trialBegins+1,1:handles.freqN,1:handles.spgl) = handles.spgDeleaked(i,trialBegins:handles.nTrials,1:handles.freqN,1:handles.spgl);
        handles.meanSpgDeleaked (i,1:handles.freqN,1:handles.spgl) = mean (tmp,1);
    end
    tmp (1:handles.nTrials-trialBegins+1,1:handles.freqN,1:handles.spgl) = handles.spg(i,trialBegins:handles.nTrials,1:handles.freqN,1:handles.spgl);
    handles.meanSpg (i,1:handles.freqN,1:handles.spgl) = mean (tmp,1);
    
    % compute mean E/I index
    tmp = [];
    tmp(1:handles.nTrials-trialBegins+1,1:handles.spgl) = handles.EIrat(i,trialBegins:handles.nTrials,1:handles.spgl);
    handles.meanEIrat (i,1:handles.spgl) = mean (tmp,1);    
    % compute mean 'gamma' power
    tmp = [];
    tmp(1:handles.nTrials-trialBegins+1,1:handles.spgl) = handles.spgPlot(i,trialBegins:handles.nTrials,1:handles.spgl);
    handles.meanSpgPlot (i,1:handles.spgl) = mean (tmp,1);
    
    % the computation of the mean power spectra has been moved to the
    % ZebraSpectre.m file    
end
zebraguidata(hObject, handles);
end