function data = zebraguidata(h, data_in)

% CHECKING
% Same checks performed by guidata
narginchk(1,2);
if (nargin == 2)
  error(nargoutchk(0, 0, nargout, 'struct'));
end

fig = [];
if isscalar(h) &&  ishghandle(h)
  fig = getParentFigure(h);
end
if isempty(fig)
  error(message('MATLAB:guidata:InvalidInput'));
end


% GET PART
% The get part is the same as guidata.

if nargin == 1 % GET
    data = guidata(h);

% SET PART
% The set part is different from guidata
else % (nargin == 2) SET
    % Here we have to remove the data. Before guidata!

    if strcmp(fig.Tag,'ZebraMainFig')
        mainfigcheck = 1; % We are in main figure
    else
        mainfigcheck = 0; % We are in non-main figures
    end
    
% Cleaning data_in from all the possible nested parentHandles (only the bad ones)
    structlvl = 1;
    data_in = parentHandles_search(data_in, mainfigcheck, structlvl);
    
    guidata(h, data_in)
end


function mystruct = parentHandles_search(mystruct, mainfigcheck, structlvl)
% structlvl is 1 at the base level (handles structure), when you dig deeper in the
% structures, it gets higher (2 is for handles.autoData structure).
% Structure levels after 3 are not needed.

% If nothing is found, the input list is returned

% This conditional block is executed early, before reaching the structure
% level where the remotion of parentHandles will be performed.
if structlvl < (3 - mainfigcheck)
    structlist = fieldnames(mystruct);
    if ~isempty(structlist)
        for i=1:length(structlist)
            
            tmpfield = mystruct.(structlist{i});
            
% If a field is a structure, assign the value of the function parentHandles_search to the temporary structure
            if isstruct(tmpfield)
                tmpfield = parentHandles_search(tmpfield,mainfigcheck,structlvl+1);
                
% Assign the value in tmpfield to the field structlist(i) of mystruct
% 'caller' means that the variable will be in the function workspace
                mystruct.(structlist{i}) = tmpfield;
            end
        end
    end
    
% This executes when the correct struct level for the remotion of
% parentHandles has been reached.
% This level is 2 for the main figure and is 3 for the other figures
elseif structlvl == (3 - mainfigcheck)
    if isfield(mystruct, 'parentHandles')
        mystruct = rmfield(mystruct, 'parentHandles');
    end
end


function fig = getParentFigure(fig)
% if the object is a figure or figure descendent, return the
% figure.  Otherwise return [].
while ~isempty(fig) && ~strcmp('figure', get(fig,'Type'))
  fig = get(fig,'Parent');
end
