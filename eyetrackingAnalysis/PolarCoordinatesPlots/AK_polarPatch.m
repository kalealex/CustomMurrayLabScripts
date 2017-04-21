function [ hPatch ] = AK_polarPatch( varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% CHECK INPUT

% check for axes handle
if ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes') &&...
        isa(varargin{2},'double') && isa(varargin{3},'double') && isa(varargin{4},'double') && isa(varargin{5},'double') &&...
        isa(varargin{6},'double') && isa(varargin{7},'double') &&...
        length(varargin{2}) == length(varargin{3}) && length(varargin{2}) == length(varargin{4}) && length(varargin{2}) == length(varargin{5}) && length(varargin{2}) == length(varargin{7})
    ah = varargin{1};
    distM = varargin{2};
    angleM = varargin{3};
    distSD = varargin{4};
    angleSD = varargin{5};
    colors = varargin{6};
    colorIdx = varargin{7};
    varargin(1:7) = [];
    newAx = false;  
elseif isa(varargin{1},'double') && isa(varargin{2},'double') && isa(varargin{3},'double') && isa(varargin{4},'double') &&...
        isa(varargin{5},'double') && isa(varargin{6},'double') &&...
        length(varargin{1}) == length(varargin{2}) && length(varargin{1}) == length(varargin{3}) && length(varargin{1}) == length(varargin{4}) && length(varargin{1}) == length(varargin{6})
    ah = gca;
    distM = varargin{1};
    angleM = varargin{2};
    distSD = varargin{3};
    angleSD = varargin{4};
    colors = varargin{5};
    colorIdx = varargin{6};
    varargin(1:6) = [];
    % if the axes have children, it's not new (important for adjusting
    % limits below)
    newAx = isempty(get(ah,'Children'));
else 
    error(['Input should be as follows: '...
        'axis handle [optional],'...
        'distance means (vector),'...
        'angle means (vector),'...
        'distance SDs (vector),'...
        'angle SDs (vector),'...
        'colors (matrix w/ dimensions length(vector) by RGB),'...
        'color index (vector)']); % message
end

%% convert and plot data

% remove NaNs from data
badData = ~isfinite(distM) | ~isfinite(angleM) | ~isfinite(distSD) | ~isfinite(angleSD);
distM(badData) = [];
angleM(badData) = [];
distSD(badData) = [];
angleSD(badData) = [];
colorIdx(badData) = [];

if ~isempty(distM) && ~isempty(angleM) && ~isempty(distSD) && ~isempty(angleSD)
    % convert polar units to x,y
    [Xdata,Ydata] = pol2cart(angleM',distM'); % will throw error if there are not an equal number of conditions represented in each set
    
    % scale patch width by window size
    xl = xlim; yl = ylim;
    xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);
    K = .005*mean([xrange yrange]); % arbitrary constant multiplier for angleSD representation
    
    % scale patch color saturation
    saturation = .30;
    
    % generate patch coordinates
    patchX = nan(length(colorIdx),4);
    patchY = nan(length(colorIdx),4);
    patchC = nan(length(colorIdx),3);
    for iD = 1:length(colorIdx)
        % assign x and y coordinates for patch relative to origin
        clear Xorigin Yorigin
        Xorigin = [-distSD(iD) distSD(iD) distSD(iD) -distSD(iD)];
        Yorigin = [K*tan(angleSD(iD)) K*tan(angleSD(iD)) -K*tan(angleSD(iD)) -K*tan(angleSD(iD))];
        % rotate coordinates at angleM to match lines in compass plot; then
        % translate to xy-coordinates of vector average
        [patchX(iD,:), patchY(iD,:)] = AK_rotateTranslate2D(Xorigin,Yorigin,angleM(iD),Xdata(iD),Ydata(iD));
        % color scaling
        clear tempColor
        tempColor = rgb2hsv(colors(colorIdx(iD),:));
        tempColor(2) = tempColor(2)*saturation;
        patchC(iD,:) = hsv2rgb(tempColor);
    end
    
    % set axis limits:
    maxX = max(max(patchX))+.05*max(max(patchX)); maxY = max(max(patchY))+.05*max(max(patchY)); % find max values and use to create axis limits
    axis([-maxX maxX -maxY maxY])
    
    for iData = 1:length(colorIdx)
        % plot patches
        hPatch = patch(patchX(iData,:),patchY(iData,:),patchC(iData,:),'EdgeColor','w');
        % plot vector averages if new axis
        if newAx 
            hPl(iData) = scatter(ah,Xdata(iData),Ydata(iData),'Marker','o');
            set(hPl(iData),'Color',colors(colorIdx(iData),:)) % color by condition
        end
    end
    
end

end

