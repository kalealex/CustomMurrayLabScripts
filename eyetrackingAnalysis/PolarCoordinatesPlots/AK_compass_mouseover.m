function [ ph ] = AK_compass_mouseover( varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% CHECK INPUT

% check for axes handle
if ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes') &&...
        isa(varargin{2},'double') && isa(varargin{3},'double') && isa(varargin{4},'double') && isa(varargin{5},'double') &&...
        ~iscell(varargin{6}) &&...
        length(varargin{2}) == length(varargin{3}) && length(varargin{2}) == length(varargin{5})
    ah = varargin{1};
    dist = varargin{2};
    angle = varargin{3};
    colors = varargin{4};
    colorIdx = varargin{5};
    mouseover_labels = [];
    varargin(1:5) = [];
    newAx = false;
elseif ~iscell(varargin{1}) && length(varargin{1}) == 1 && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes') &&...
        isa(varargin{2},'double') && isa(varargin{3},'double') && isa(varargin{4},'double') && isa(varargin{5},'double') &&...
        iscell(varargin{6}) && isa(varargin{7},'double') &&...
        length(varargin{2}) == length(varargin{3}) && length(varargin{2}) == length(varargin{5}) && length(varargin{2}) == length(varargin{7})
    ah = varargin{1};
    dist = varargin{2};
    angle = varargin{3};
    colors = varargin{4};
    colorIdx = varargin{5};
    mouseover_labels = varargin{6};
    labelIdx = varargin{7};
    varargin(1:7) = [];
    newAx = false;    
elseif isa(varargin{1},'double') && isa(varargin{2},'double') && isa(varargin{3},'double') && isa(varargin{4},'double') &&...
        iscell(varargin{5}) && isa(varargin{6},'double') &&...
        length(varargin{1}) == length(varargin{2}) && length(varargin{1}) == length(varargin{4}) && length(varargin{1}) == length(varargin{6})
    ah = gca;
    dist = varargin{1};
    angle = varargin{2};
    colors = varargin{3};
    colorIdx = varargin{4};
    mouseover_labels = varargin{5};
    labelIdx = varargin{6};
    varargin(1:6) = [];
    % if the axes have children, it's not new (important for adjusting
    % limits below)
    newAx = isempty(get(ah,'Children'));
else 
    ah = gca;
    dist = varargin{1};
    angle = varargin{2};
    colors = varargin{3};
    colorIdx = varargin{4};
    mouseover_labels = [];
    varargin(1:4) = [];
    % if the axes have children, it's not new (important for adjusting
    % limits below)
    newAx = isempty(get(ah,'Children'));
end

%% convert and plot data

% remove NaNs from data
badData = ~isfinite(dist) | ~isfinite(angle);
dist(badData) = [];
angle(badData) = [];
colorIdx(badData) = [];
if ~isempty(mouseover_labels)
    labelIdx(badData) = [];
end 

if ~isempty(dist) && ~isempty(angle)
    % convert polar units to x,y
    [Xdata,Ydata] = pol2cart(angle',dist'); % will throw error if there are not an equal number of conditions represented in each set
    
    % set as input for line function
    lineOrigin = zeros(size(Xdata));
    lineX = horzcat(lineOrigin,Xdata);
    lineY = horzcat(lineOrigin,Ydata);
    
    % set axis limits:
    if newAx
        maxX = max(Xdata)+.05*max(Xdata); maxY = max(Ydata)+.05*max(Ydata); % find max values and use to create axis limits
        axis([-maxX maxX -maxY maxY])
    end
    
    % plot data using compass command
%     ph = compass(ah,Xdata,Ydata); % plot
    for iData = 1:length(colorIdx)
        ph(iData) = line(lineX(iData,:),lineY(iData,:),'Marker','o');
        set(ph(iData),'Color',colors(colorIdx(iData),:)) % color by condition
        
%         text(Xdata(iData),Ydata(iData),mouseover_labels{labelIdx(iData)},... % constant labels
%         'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
%         'hittest','off');
    end
    
end


%apply mouse motion function
if ~isempty(mouseover_labels)
    set(gcf,'windowbuttonmotionfcn',{@mousemove,ph,Xdata,Ydata,mouseover_labels,labelIdx});
end
 

function mousemove(src,ev,L,Xdata,Ydata,labels,labelIdx)
 
%since this is a figure callback, the first input is the figure handle:
f=handle(src);
 
%like all callbacks, the second input, ev, isn't used. 
 
%determine which object is below the cursor:
obj=hittest(f); %<-- the important line in this demo
 
if any(any(obj==L)) %if over any plot... 
    % change pointer when hovering over plot
    set(f,'Pointer','crosshair');
    
    % determine which plot
    Lindex = find(obj==L);
%     if mod(Lindex,length(L(:,1))) ~= 0
%         Lrow = mod(Lindex,length(L(:,1)));
%     else
%         Lrow = length(L(:,1));
%     end
%     Lcol = ceil(Lindex/length(L(:,1)));
    
    %get cursor coordinates in its axes:
    a=get(L(Lindex),'parent');
    point=get(a,'currentpoint');
    xclick=point(1,1,1);
    yclick=point(1,2,1);
 
    %determine which point we're over:
    idx=findclosestpoint2D(xclick,yclick,L(Lindex));
 
    %make a "tool tip" that displays data point
    xl = xlim; yl = ylim;
    xrange = xl(2)-xl(1); yrange = yl(2)-yl(1);
    xoffset=xrange/1000;
    yoffset=yrange/1000;
 
    delete(findobj(f,'tag','mytooltip')); %delete last tool tip
    text(Xdata(Lindex)+xoffset,Ydata(Lindex)+yoffset,labels{labelIdx(Lindex)},...
        'backgroundcolor',[1 1 .8],'tag','mytooltip','edgecolor',[0 0 0],...
        'hittest','off');
else
    % change pointer back and delete plot
    set(gcf,'Pointer','arrow');
    delete(findobj(f,'tag','mytooltip')); %delete last tool tip
 
end
 
 
function index=findclosestpoint2D(xclick,yclick,datasource)
%this function checks which point in the plotted line "datasource"
%is closest to the point specified by xclick/yclick. It's kind of 
%complicated, but this isn't really what this demo is about...
 
xdata=get(datasource,'xdata');
ydata=get(datasource,'ydata');
 
activegraph=get(datasource,'parent');
 
pos=getpixelposition(activegraph);
xlim=get(activegraph,'xlim');
ylim=get(activegraph,'ylim');
 
%make conversion factors, units to pixels:
xconvert=(xlim(2)-xlim(1))/pos(3);
yconvert=(ylim(2)-ylim(1))/pos(4);
 
Xclick=(xclick-xlim(1))/xconvert;
Yclick=(yclick-ylim(1))/yconvert;
 
Xdata=(xdata-xlim(1))/xconvert;
Ydata=(ydata-ylim(1))/yconvert;
 
Xdiff=Xdata-Xclick;
Ydiff=Ydata-Yclick;
 
distnce=sqrt(Xdiff.^2+Ydiff.^2);
 
index=distnce==min(distnce);
 
index=index(:); %make sure it's a column.
 
if sum(index)>1
    thispoint=find(distnce==min(distnce),1);
    index=false(size(distnce));
    index(thispoint)=true;
end


