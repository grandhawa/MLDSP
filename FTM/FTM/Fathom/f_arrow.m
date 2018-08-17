function hh=arrow5(varargin)
% - draw an arrow
%   ARROW(tail,head) draw an arrow between the tail (x,y) and head (x,y)
%
% ARROW(tail,head,'Property1',PropertyValue1,'Property2',PropertyValue2,...)
%   draw an arrow, and sets the values of the specified properties of the arrow.
%
%   ARROW('Property1',PropertyValue1,'Property2',PropertyValue2,...)
%   draw an arrow, and set the values of the specified properties of the arrow.
%   The 'head' and 'tail' properties must be specified.
%
%   ARROW(H,'Property1',PropertyValue1,'Property2',PropertyValue2,...)
%   sets the values of the specified properties of an existing arrow.
%
%   H = ARROW(...) returns the handle to the line object used as the arrow.
%
%   ARROW  redraws all arrows (rescaling if required). This may affect the
%   the values of gca and gcf.
%
%   Additional properties supported are:
%     'size'   length of the arrow head, in inches (default 0.125)
%     'direction'   direction of the arrow, 'forward', 'backward', or 'double'
%     'angle'  angle of the arrowhead
%     'tail'   define [x,y] of tail of arrow
%     'head'   define [x,y] of head of arrow
%
%   NOTES: The arrowhead is scaled to be drawn correctly when printed. It may
%     appear skewed in the figure window. Also, any of the following may cause
%     existing arrows to be rendered improperly:
%       - the axes 'position' is changed
%       - the axes 'xlim' or 'ylim' is changed
%       - the figure 'PaperPosition' is changed
%     In general, it is recommended that arrows be drawn last.
%     Use ARROW with no arguments to redraw all arrows
%
% Posted to news://comp.soft-sys.matlab, March-2003
% by AJ \"no z\" Johnson (aj.jozhnson@lmco.com)


h = [];
head = [];
tail = [];

if nargin == 0
   redraw_all_arrows; % sub-fuction speified below
   return
end

if mod(nargin,2) == 0
   % Even number of arguments
   if isreal(varargin{1}) & isreal(varargin{2})
      head = varargin{1};
      tail = varargin{2};
      varargin(1:2) = [];
   end
else
   % Odd number of arguments
   if isreal(varargin{1})
      % First arg is is: hopefully handle(s) to existing arrow(s)
      h = varargin{1};
      tag = get(h,'tag');
      if ~strcmp(tag,'arrow')
         error('Handle must an existing arrow')
      end
      varargin(1) = [];
      varargin = [get(h,'userdata'),varargin];
   else
      error('Incorrect number of input arguments')
   end
end


keeplist = [];
discardlist = [];
headsize = 0.125;
ang = 15; % Angle of arrow
direction = 'forward';
for i=1:2:length(varargin)
   if ~isstr(varargin{i})
      error('Invalid property name');
   else
      switch lower(varargin{i})
         case 'direction', direction = varargin{i+1};
         case 'size', headsize = varargin{i+1};
         case 'angle', ang = varargin{i+1};
         case 'head', head = varargin{i+1};
         case 'tail', tail = varargin{i+1};
         otherwise, keeplist = [keeplist,i,i+1];
      end
   end
end

varargin = varargin(keeplist);


x1 = head(1);
x2 = tail(1);
y1 = head(2);
y2 = tail(2);
%  x=[head(1),tail(1)]
%  y=[head(2),tail(2)]
if length(h) == 0
   h = plot([x1 x2],[y1 y2]); % Plot the basic arrow line
   a = gca; % Parent axes handle
   f = gcf; % Parent figure handle
else
   a = get(h,'parent');
   f = get(a,'parent');
end

xlim = get(a,'xlim');
ylim = get(a,'ylim');
xd = diff(xlim);
yd = diff(ylim);

ppos_inch = get(f,'PaperPosition'); % Figure position on paper, in inches (left, bottom, width, height)
apos_norm = get(a,'Position'); % Axes position on figure, normalized
apos_inch = [
   ppos_inch(1)+apos_norm(1)*ppos_inch(3),
   ppos_inch(2)+apos_norm(2)*ppos_inch(4),
   apos_norm(3)*ppos_inch(3),
   apos_norm(4)*ppos_inch(4)];
xscale = xd/apos_inch(3);
yscale = yd/apos_inch(4);





len = sqrt(sum((head-tail).^2)); % Length of arrow

baseline = (tail-head)./len; % unit vector long baseline
baseline = baseline(1)/xscale+j*baseline(2)/yscale;
baseline = baseline ./ abs(baseline);
orient = headsize.*exp(j*ang*pi/180);

p1=orient.*baseline;
p2=conj(orient).*baseline;

xp1=real(p1)*xscale;
yp1=imag(p1)*yscale;
xp2=real(p2)*xscale;
yp2=imag(p2)*yscale;

if strcmp(direction, 'backward') | strcmp(direction, 'double')
   x = [x1,x1+xp1,x1,x1+xp2,x1];
   y = [y1,y1+yp1,y1,y1+yp2,y1];
else
   x=[x1];
   y=[y1];
end

if strcmp(direction, 'forward') | strcmp(direction, 'double')
   x = [x, x2,x2-xp1,x2,x2-xp2,x2];
   y = [y, y2,y2-yp1,y2,y2-yp2,y2];
else
   x = [x, x2];
   y = [y, y2];
end

userdata = {'tail', tail, 'head', head, 'direction', direction, 'size', headsize, 'angle', ang};
set(h,'tag','arrow','userdata',userdata,'xdata',x','ydata',y,varargin{:});

if nargout > 0
   hh = h;
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_all_arrows(varargin)

for f=get(0,'children')' % Loop over all figures
   for a=get(f,'children')' % Loop over all axes
      h_list = get(a,'children');
      for h=h_list' % Loop over all types
         Type = get(h,'Type');
         if isstr(Type) & strcmp(Type,'line')
            if strcmp(get(h,'Tag'),'arrow') & length(get(h,'userdata')) > 6
               args = [{h},varargin];
               f_arrow(args{:}); % Redraw arrow
            end
         end
      end
   end
end

