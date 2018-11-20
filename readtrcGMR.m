% [DATAOUT]=readtrcGMR(SETTIINGS)  -   Reads Micromed System 98 *.trc EEG File 
%                                   for import into EEGLAB. This file reads 
%                                   Micromed System Plus EEG *.trc files with 
%                                   header of type 4
%
% USAGE:
%   >> [DATAOUT]=readtrcGMR(SETTINGS)
%
%
% INPUT
%   SETTINGS is a struct holding the parameters for reading the .TRC file
%   SETTINGS has the following fields:
%
%       SETTINGS.filename :                 Name of file to be imported
%       SETTINGS.loadevents.state :         'yes' for loading event triggers
%                                           'no' for not
%       SETTINGS.loadevents.type :          'marker' for event triggers inserted 
%                                           on 'MKR channel
%                                           'eegchan' for triggers inserted on 
%                                           EEG channels
%                                           'none' or
%                                           'both'
%       SETTINGS.loadevents.dig_ch1:        number of name of eegchan marker channel 1
%       SETTINGS.loadevents.dig_ch1_label:  label to give events on eegchan marker channel 1
%       SETTINGS.loadevents.dig_ch2:        number of name of eegchan marker channel 2
%       SETTINGS.loadevents.dig_ch2_label:  label to give events on eegchan marker channel 2
%       SETTINGS.chan_adjust_status:        1 for adjusting amp of channels 0 for not
%       SETTINGS.chans                      channels to load, [ ] for all
%       (default)
%       SETTINGS.chan_adjust                channels to adjust
%
%   Alternant method: enter 11 input arguments each corresponding to one of the above
%   fields.  If only one arg is used,  it must be the struct above.  If more, there
%   must be 11 inputs in the order above i.e. OUT=readtrc(filename,eventstate....etc).
%
% OUTPUT
%   TRC with same fields as EEGLAB EEG file structure. (see eeg_checkset.m)
%   The following fields are use: 
%       TRC.data 
%       TRC.filename
%       TRC.filepath
%       TRC.srate
%       TRC.setname
%       TRC.pnts
%       TRC.nbchan
%       TRC.trials
%       TRC.xmin
%       TRC.ref
%
% Author: Rami K. Niazy
% Copyright (c) University of Oxford.
%
%   see also pop_readtrc.m    eegplugin_trcimport.m

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004-2005 University of Oxford
% Author: Rami K. Niazy
%         rami@fmrib.ox.ac.uk
%
% This program can only be provided from Micromed s.r.l.  You may not
% give away or edit this program. 

% April 19, 2005
% Version 1.1
% Reads triggers from 'MRK' and EEG channels together
% Assigns value of 'MRK' trigger as event 'type' in EEGLAB
% Reads only selected channels into memory

% Dec 23, 2004
% Fixed bug for No MRK case

% Sep 29, 2004
% Fixed Trigger offset

% Sep 22, 2004
% Fixed Trigger offset

% Sep 21, 2004
% Allow usage of ':' in specifying
% channels to exclude or adjust

% Aug 12, 2004
% Added command line usage
% fixed com out for EEGLAB

% Aug 10, 2004
% Adjusted Licence
% Added Chan exclude
% Fixed scaling

% Aug 4, 2004
% Fixed scaling of data to uV
% Fixed Analog Trig size detection

% Date: Julyl 30, 2004
% Initial Setup

% Date: January 14, 2009
% Version 1.2
% Write in comments field of EEG structure some informations for export TRC plugin
% Author: Orsato Raffaele
%         raffaele.orsato@micromed-it.com  

% Date: January 20, 2009
% Close file before return

% Date: January 28, 2009
% Put all informations for the trc export under field comment

% Date: August 2017
% Modified for usage within Zebra Explore by GMR

function [TRC]=readtrcGMR(varargin)

% Definition of the eeg_emptyset structure. August 2017.
 
eeg_emptyset.setname = '';
eeg_emptyset.filename = '';
eeg_emptyset.filepath = '';
eeg_emptyset.subject = '';
eeg_emptyset.group = '';
eeg_emptyset.condition = '';
eeg_emptyset.session = [];
eeg_emptyset.comments = '';
eeg_emptyset.nbchan = 0;
eeg_emptyset.trials = 0;
eeg_emptyset.pnts = 0;
eeg_emptyset.srate = 1;
eeg_emptyset.xmin = 0;
eeg_emptyset.xmax = 0;
eeg_emptyset.times = [];
eeg_emptyset.data = [];
eeg_emptyset.icaact = [];
eeg_emptyset.icawinv = [];
eeg_emptyset.icasphere = [];
eeg_emptyset.icaweights = [];
eeg_emptyset.icachansind = [];
eeg_emptyset.chanlocs = [];
eeg_emptyset.urchanlocs = [];
eeg_emptyset.chaninfo = [];
eeg_emptyset.ref = [];
eeg_emptyset.event = [];
eeg_emptyset.urevent = [];
eeg_emptyset.eventdescription = {};
eeg_emptyset.epoch = [];
eeg_emptyset.epochdescription = {};
eeg_emptyset.reject = [];
eeg_emptyset.stats = [];
eeg_emptyset.specdata = [];
eeg_emptyset.specicaact = [];
eeg_emptyset.splinefile = '';
eeg_emptyset.icasplinefile = '';
eeg_emptyset.dipfit = [];
eeg_emptyset.history = '';
eeg_emptyset.saved = 'no';
eeg_emptyset.etc = [];

% Checking the input arguments
if nargin<1
    error('Not enough input arguments. Please see help file for usage');
elseif nargin==1
    if isstruct(varargin{1})
        PARAM=varargin{1};
    else
        error('Single input usage must be a Structure.  Please see help file.');
    end
elseif nargin < 10
    error('Not enough input arguments. Please see help file for usage');
elseif nargin > 10
    error('Too many input arguments.  Please see help file for usage');
else
    PARAM.filename=varargin{1};
    PARAM.loadevents.state=varargin{2};
    PARAM.loadevents.type=varargin{3};
    PARAM.loadevents.dig_ch1=varargin{4};
    PARAM.loadevents.dig_ch1_label=varargin{5};
    PARAM.loadevents.dig_ch2=varargin{6};
    PARAM.loadevents.dig_ch2_label=varargin{7};
    PARAM.chan_adjust_status=varargin{8};
    PARAM.chan_adjust=varargin{9};
    PARAM.chans=varargin{10};
end
Trigs1=[];
Trigs2=[];



% ---------------- Opening File------------------
trcfile=PARAM.filename;
fid=fopen(trcfile,'r');
if fid==-1
    error('Can''t open *.trc file')
end
TRC=eeg_emptyset;           % eeg_emptyset is explicitly defined above. 

%------------------reading patient & recording info----------
fseek(fid,64,-1);
surname=char(fread(fid,22,'char'))';        
name=char(fread(fid,20,'char'))';           

fseek(fid,128,-1);
day=fread(fid,1,'char');                    
if length(num2str(day))<2
    day=['0' num2str(day)];
else
    day=num2str(day);
end
month=fread(fid,1,'char');
switch month
case 1 
    month='JAN';
case 2 
    month='FEB';
case 3 
    month='MAR';
case 4 
    month='APR';
case 5 
    month='MAY';
case 6 
    month='JUN';
case 7 
    month='JUL';
case 8 
    month='AUG';
case 9 
    month='SEP';
case 10 
    month='OCT';
case 11 
    month='NOV';
case 12 
    month='DEC';            
end
year=num2str(fread(fid,1,'char')+1900);     

%------------------ Reading Header Info ---------

fseek(fid,175,-1);
Header_Type=fread(fid,1,'char');
if Header_Type ~= 4
    error('*.trc file is not Micromed System98 Header type 4')
end

fseek(fid,138,-1);
Data_Start_Offset=fread(fid,1,'uint32');
Num_Chan=fread(fid,1,'uint16');
Multiplexer=fread(fid,1,'uint16');
Rate_Min=fread(fid,1,'uint16');
Bytes=fread(fid,1,'uint16');
fseek(fid,176+8,-1);
Code_Area=fread(fid,1,'uint32');
Code_Area_Length=fread(fid,1,'uint32');
fseek(fid,192+8,-1);
Electrode_Area=fread(fid,1,'uint32');
Electrode_Area_Length=fread(fid,1,'uint32');

fseek(fid,400+8,-1);
Trigger_Area=fread(fid,1,'uint32');
Tigger_Area_Length=fread(fid,1,'uint32');


%----------------- Allocate Memory and Determine Data Type----------

fseek(fid,Data_Start_Offset,-1);

fprintf('Allocating memory...\n');
switch Bytes
case 1   
    bstring='uint8';
case 2
    bstring='uint16';
case 4
    bstring='uint32';
end
       
tracetmp=fread(fid,bstring,(Num_Chan-1)*Bytes)';
traceL=length(tracetmp);
clear tracetmp;

if isempty(PARAM.chans)         % PARAM.chans must be a string containing the channels separated by spaces
    chans=1:Num_Chan;
else
    chans=eval([ '[' PARAM.chans ']' ]);
end

chansL=length(chans);
tracedata=zeros(chansL,traceL);
m=traceL;

%------------------ Reading Code Info -------------
fseek(fid,Code_Area,-1);
code=fread(fid,Num_Chan,'uint16');


for c=1:Num_Chan
    electrode(c).chan_record=code(c);
    fseek(fid,Electrode_Area+code(c)*128,-1);
    fseek(fid,2,0);
    if c <10 
        electrode(c).positive_input=...
            [num2str(c),' -',char(fread(fid,6,'char'))'];
        electrode(c).negative_input=...
            [num2str(c),' -',char(fread(fid,6,'char'))'];
    else
        electrode(c).positive_input=...
            [num2str(c),'-',char(fread(fid,6,'char'))'];
        electrode(c).negative_input=...
            [num2str(c),'-',char(fread(fid,6,'char'))'];
    end
    electrode(c).logical_min=fread(fid,1,'int32');
    electrode(c).logical_max=fread(fid,1,'int32');
    electrode(c).logical_ground=fread(fid,1,'int32');
    electrode(c).physical_min=fread(fid,1,'int32');
    electrode(c).physical_max=fread(fid,1,'int32');
    
    electrode(c).measurement_unit=fread(fid,1,'int16');
    switch electrode(c).measurement_unit
    case -1
        electrode(c).measurement_unit=1e-9; 
    case 0
        electrode(c).measurement_unit=1e-6;
    case 1
        electrode(c).measurement_unit=1e-3;
    case 2
        electrode(c).measurement_unit=1;
    case 100
        electrode(c).measurement_unit='percent';
    case 101
        electrode(c).measurement_unit='bpm';
    case 102
        electrode(c).measurement_unit='Adim';
    otherwise
        warning('Unknown measurement unit. uV assumed.');
        electrode(c).measurement_unit=10e-6;
    end
    fseek(fid,8,0);
    electrode(c).rate_coef=fread(fid,1,'uint16'); 
end


%---------------- Read & Prep Trigger Area Data ----------
fseek(fid,Trigger_Area,-1);
for l=1:Tigger_Area_Length/6
    trigger(1,l)=fread(fid,1,'uint32');
    trigger(2,l)=fread(fid,1,'uint16');
end


first_trigger=trigger(1,1);
tl=length(trigger);
NoTrig=0;
for tr=1:tl
    if ((trigger(1,tr) <= m) & (trigger(1,tr) >= first_trigger))
        NoTrig=NoTrig+1;
    end
end

if NoTrig > 0
   	trigger=trigger(:,1:NoTrig);
else
	trigger=[];
	first_trigger=[];
end

    

%---------------Reading Other Event Data   -------------
TRC=eeg_emptyset;
switch  PARAM.loadevents.state
case 'no';
case 'yes';
    fprintf('Extracting events...\n');
    switch PARAM.loadevents.type
    case 'marker'
        if ~isempty(trigger)            % xche rigger is empty???
            [triggerR,triggerC]=size(trigger)
            for E=1:triggerC
                TRC.event(end+1).type=num2str(trigger(2,E));
                TRC.event(end).latency=trigger(1,E)+1;
            end
        else
            warndlg('No marker triggers to import on ''MRK'' channel',...
                'Import .TRC Warning!');
        end
    case 'eegchan'
        %------find 1st trigger channel-------
        dig_ch=[];
        if str2num(PARAM.loadevents.dig_ch1)
            dig_ch(1)=str2num(PARAM.loadevents.dig_ch1);
        else
            for C=1:Num_Chan
                if ~isempty(findstr(lower(PARAM.loadevents.dig_ch1),...
                        lower(electrode(C).positive_input)))
                    dig_ch(1)=C;
                    break;
                end
            end
        end
        
        %------check for 2nd trigger channel-------
        if ~isempty(PARAM.loadevents.dig_ch2)
            if str2num(PARAM.loadevents.dig_ch2)
                dig_ch(2)=str2num(PARAM.loadevents.dig_ch2);
            else
                for C=1:Num_Chan
                    if ~isempty(findstr(lower(PARAM.loadevents.dig_ch2),...
                            lower(electrode(C).positive_input)))
                        dig_ch(2)=C;
                        break;
                    end
                end
            end
        end
        
                
        %--------read trigs-------------------------
  
        fseek(fid,Data_Start_Offset+(dig_ch(1)-1)*Bytes,-1);
        tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
        if ischar(electrode(dig_ch(1)).measurement_unit)==0
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min)*...
                electrode(dig_ch(1)).measurement_unit;
        else
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min);
        end
        trigval=min(tracedata(1,:));    % rtacedata is the trigger signal conditioned in such a way that the 
                                        % triggers are seen as negative peaks
        Trigs1=find(tracedata(1,:)==trigval);       % these are the indexes of the trigger events
                                                   
        if length(dig_ch)>1
            fseek(fid,Data_Start_Offset+(dig_ch(2)-1)*Bytes,-1);
            tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
            if ischar(electrode(dig_ch(2)).measurement_unit)==0
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min)*...
                    electrode(dig_ch(2)).measurement_unit;
            else
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min);
            end
            Trigs2=find(tracedata(1,:)==trigval);
        end
        
      
        %--------write trigs----------------------------------
        for E=1:length(Trigs1)
            TRC.event(end+1).type=PARAM.loadevents.dig_ch1_label;
            TRC.event(end).latency=Trigs1(E);
        end
            
        if ~isempty(Trigs2)
            for E=1:length(Trigs2)
                TRC.event(end+1).type=PARAM.loadevents.dig_ch2_label;
                TRC.event(end).latency=Trigs2(E);
            end
          
        end
        
    case 'both'
        
        %-----------Marker trigs----------------------------
        if ~isempty(trigger)
            [triggerR,triggerC]=size(trigger);
            for E=1:triggerC
                TRC.event(end+1).type=num2str(trigger(2,E));
                TRC.event(end).latency=trigger(1,E)+1;
            end
        else
            warndlg('No marker triggers to import on ''MRK'' channel',...
                'Import .TRC Warning!');
        end
        
         %------find 1st trigger channel-------
        dig_ch=[];
        if str2num(PARAM.loadevents.dig_ch1)
            dig_ch(1)=str2num(PARAM.loadevents.dig_ch1);
        else
            for C=1:Num_Chan
                if ~isempty(findstr(lower(PARAM.loadevents.dig_ch1),...
                        lower(electrode(C).positive_input)))
                    dig_ch(1)=C;
                    break;
                end
            end
        end
        
        %------check for 2nd trigger channel-------
        if ~isempty(PARAM.loadevents.dig_ch2)
            if str2num(PARAM.loadevents.dig_ch2)
                dig_ch(2)=str2num(PARAM.loadevents.dig_ch2);
            else
                for C=1:Num_Chan
                    if ~isempty(findstr(lower(PARAM.loadevents.dig_ch2),...
                            lower(electrode(C).positive_input)))
                        dig_ch(2)=C;
                        break;
                    end
                end
            end
        end
        
                
        %--------read trigs-------------------------
       
        fseek(fid,Data_Start_Offset+(dig_ch(1)-1)*Bytes,-1);
        tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
        if ischar(electrode(dig_ch(1)).measurement_unit)==0
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min)*...
                electrode(dig_ch(1)).measurement_unit;
        else
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min);
        end
        trigval=min(tracedata(1,:));
        Trigs1=find(tracedata(1,:)==trigval);
        

        
        if length(dig_ch)>1
            fseek(fid,Data_Start_Offset+(dig_ch(2)-1)*Bytes,-1);
            tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
            if ischar(electrode(dig_ch(2)).measurement_unit)==0
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min)*...
                    electrode(dig_ch(2)).measurement_unit;
            else
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min);
            end
            Trigs2=find(tracedata(1,:)==trigval);
        end
        
        %------------write trigs----------------------------------
        
        for E=1:length(Trigs1)
            TRC.event(end+1).type=PARAM.loadevents.dig_ch1_label;
            TRC.event(end).latency=Trigs1(E);
        end
            
        if ~isempty(Trigs2)
            for E=1:length(Trigs2)
                TRC.event(end+1).type=PARAM.loadevents.dig_ch2_label;
                TRC.event(end).latency=Trigs2(E);
            end
          
        end
        
    end   
end


%------------------Reading Data-------------------

fprintf('Reading data...\n');
for c=1:chansL
    fseek(fid,Data_Start_Offset+(chans(c)-1)*Bytes,-1);
    tracedata(c,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
    if ischar(electrode(chans(c)).measurement_unit)==0
        tracedata(c,:)=-((tracedata(c,:)-electrode(chans(c)).logical_ground)/...
            (electrode(chans(c)).logical_max-...
            electrode(chans(c)).logical_min+1))*...
            (electrode(chans(c)).physical_max-...
            electrode(chans(c)).physical_min)*...
            electrode(chans(c)).measurement_unit;
    else
        tracedata(c,:)=-((tracedata(c,:)-electrode(chans(c)).logical_ground)/...
            (electrode(chans(c)).logical_max-...
            electrode(chans(c)).logical_min+1))*...
            (electrode(chans(c)).physical_max-...
            electrode(chans(c)).physical_min);
    end
    
    if ~isempty(Trigs1)
        if chans(c)==dig_ch(1)
            tracedata(c,Trigs1)=...
                (tracedata(c,(Trigs1+1))+tracedata(c,(Trigs1-1)))/2;
        end
    end

    if ~isempty(Trigs2)
        if chans(c)==dig_ch(2)
            tracedata(c,Trigs2)=...
                (tracedata(c,(Trigs2+1))+tracedata(c,(Trigs2-1)))/2;
        end
    end   
    
end


% -----------Reading Fs-------------------------

mean_fs=mean(cat(1,electrode.rate_coef));
switch mean_fs
case 1
    fs=1*Rate_Min;
case 2
    fs=2*Rate_Min;
case 3
    fs=3*Rate_Min;
case 4
    fs=4*Rate_Min;
case 5
    fs=5*Rate_Min;
otherwise
    warning('Unsupported Sampling Frequency');
end


%----------Prep output-------------------------

fprintf('Preparing output...\n');
if PARAM.chan_adjust_status==1
	if length(tracedata)>(fs*60)
        avgvar=mean(var(tracedata(1:Num_Chan-3,1:fs*60)'));
	else
        avgvar=mean(var(tracedata(1:Num_Chan-3,:)'));
	end
    
    ch_adj_t=['[' PARAM.chan_adjust ']'];
    ch_adj=eval(ch_adj_t);
    
    for ch=1:length(ch_adj)
        tracedata(ch_adj(ch),:)=tracedata(ch_adj(ch),:)*avgvar/var(tracedata(ch_adj(ch),:));
    end
end


TRC.data=-tracedata*1e6; % scale to uV and change polarity for EEGLAB
sp=findstr(surname,'  ');
if sp >=1
    TRC.setname=[surname(1:(sp-1)) ', ' name(1) '. ' year month day ' .TRC File'];
else
    TRC.setname=[surname ', ' name(1) '. ' year month day ' .TRC File'];
end

%-----------------------Prepare filename for copy in Export routine--------------------
k = findstr(trcfile,'.');
filecopy = [trcfile(1:k-1),'_copy','.TRC'];
%-----------------------write if adjust was done or not--------------------------------
if PARAM.chan_adjust_status==1
    fileop ='Adjust done'; 
else
    fileop ='Adjust undone';     
end
%-----------------------write the channels read--------------------------------
if isempty(PARAM.chans) 
    channels=['1:',int2str(length(chans))];
else
    channels=PARAM.chans;
end
TRC.filename=trcfile;
%-----------------Used for export to TRC plugin--------------------------------
% TRC.comments=['micromed','|',trcfile,'|',filecopy,'|',fileop,'|',channels];
%------------------------------------------------------------------------------
TRC.subject = [name surname];
TRC.comments = [day ' ' month ' ' year];
TRC.pnts=length(tracedata);
TRC.nbchan=chansL;
TRC.trials=1;
TRC.srate=fs;
TRC.xmin=0;
TRC.xmax=(TRC.pnts-1)/fs;
TRC.ref='common';
fclose(fid);


return;
