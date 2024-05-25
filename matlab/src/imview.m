function fg = imview(u,varargin)
%% 
% Usage: fg = imview(u,Name,Value)
%
% Input(s)/Output(s):
%
%   u  : (matrix of double) image graylevels
% 
%   fg : output figure handler
% 
% Optional Name-Value pair arguments:
% 
%   ['scale',z] : (scalar, default z = 1) upscale factor (must have integer
%                 value), the displayed image will be kron(u, ones(z)),
%                 which corresponds to a nearest-neighbor upscaling of u
%                 with factor z along both axes
%
%   ['black',b] : (scalar, default b = min(u(:))) graylevel value to
%                 display in black (or in the first color of the colormap 
%                 if not set to gray) 
% 
%   ['white',w] : (scalar, default w = max(u(:))) graylevel value to
%                 display in white (or in the last color of the colormap if
%                 not set to gray)
%
%   ['colormap',cmap] : (default, cmap = gray(256)) colormap to use for
%                       the display
%
% Interactive commands (keyboard): 
% 
%  'h' : list interactive commands
%  'r' : redraw the figure (for some reasons, MATLAB might draw the figure
%        with incorrect dimensions, this phenomenom happens randomly)
%  'q' : quit (close the figure)
% 
% Description: display and image into a Matlab Figure with tight borders
% AND precise control of the interpolation (by default, one pixel of the
% screen = one pixel of the image). 

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('u');
p.addParameter('scale',1);
p.addParameter('black',min(u(:)));
p.addParameter('white',max(u(:)));
p.addParameter('colormap',gray(256));
parse(p,u,varargin{:});
z = p.Results.scale;
m = p.Results.black;
M = p.Results.white; 
cmap = p.Results.colormap; 

%% consistency checks
% input u (2D matrix of double real numbers)
if(~isreal(u) || numel(size(u)) ~= 2)
    help imview;
    error('input ''u'' must be a 2D array of double real numbers');
end
% input z (real scalar number > 0, no decimal part)
if(~isreal(z) || ~isscalar(z) || z ~= floor(z) || z <= 0)
    help imview;
    error('optional input ''z'' must be a positive real scalar number, without decimal part (z == floor(z))');
end
% input b (real scalar number)
if(~isreal(m) || ~isscalar(m))
    help imview;
    error('optional input ''b'' must be a real scalar number');
end
% input w (real scalar number)
if(~isreal(M) || ~isscalar(M))
    help imview;
    error('optional input ''w'' must be a real scalar number');
end
% input cmap (2D array with 3 columns, containing real numbers in [0,1])
if(~isreal(cmap) || numel(size(cmap)) ~= 2 || size(cmap,2) ~= 3 || min(cmap(:)) < 0 || max(cmap(:)) > 1)
    help imview;
    error('optional input ''cmap'' must be a 2D array with three columns and values in [0,1]');
end

%% local functions (event handlers)

    function callback_keyboard_pressed(fg,evt,u,m,M) % keyboard event handler function
        switch evt.Key

            case 'h' % help (display interactive commands)
                fprintf("Interactive commands (keyboard)\n\n");
                fprintf("  'h' : list interactive commands\n"); 
                fprintf("  'r' : redraw the figure\n");
                fprintf("  'q' : quit (close the figure)\n"); 
                fprintf("\n"); 
                
            case 'q' % close figure
                close(fg);
                
            case 'r' % redraw
             disp('redraw figure');
             drawimage(u,m,M); 
        end
    end

    function fg = drawimage(u,m,M)
        gr = groot();
        u = double(u);
        [ny,nx] = size(u);
        fg = gcf();
        image(255*(u-m)/(M-m),'CDataMapping','direct','XData',0:nx-1,'YData',0:ny-1);
        axis 'image';
        set(fg,'Colormap',cmap,'KeyPressFcn',{@callback_keyboard_pressed,u,m,M});
        set(gca(),'XLim',[0,nx-1],'YLim',[0,ny-1],'Visible','off','Units','pixels','ActivePositionProperty','position','Position',[1,1,nx,ny]);
        pos = get(fg,'Position');
        pos(1:2) = [min(100,gr.ScreenSize(3)-nx),min(gr.ScreenSize(4)-ny-100,gr.ScreenSize(4)-ny)];
        set(fg,'Position',pos);
        set(fg,'Visible','on');
        pause(1e-6);
        pos(3:4) = [nx,ny];
        set(fg,'Position',pos);
    end

%% CORE OF THE MODULE
fg = drawimage(kron(u,ones(z)),m,M); 

end
