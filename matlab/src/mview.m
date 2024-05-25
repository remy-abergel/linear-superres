function fg = mview(u,varargin)
%% 
% Usage: fg = mview(u,Name,Value)
%
% Input(s)/Output(s):
%
%   u  : (hypmatrix of double) image sequence
% 
%   fg : output figure handler
% 
% Optional Name-Value pair arguments:
%
%   ['scale',z] : (scalar, default z = 1) upscale factor (must have integer
%                 value), the displayed frames will be kron(u(:,:,id),
%                 ones(z)), for id = 1..size(u,3), this corresponds to
%                 frame-by-frame nearest-neighbor upscaling with factor z
%                 along both axes
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
%  'f' : go forward (show next frame)
%  'b' : go backward (show previous frame)
%  'r' : redraw the figure (for some reasons, MATLAB might draw the figure
%        with incorrect dimensions, this phenomenom happens randomly)
%  'q' : quit (close the figure)
% 
% Description : frame-by-frame movie displayer

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

%% consistency checks (TODO)

%% local functions (event handlers)
    function callback_keyboard_pressed(fg,evt,u,m,M) % keyboard event handler function
        
        im_hdl = findobj('Tag',sprintf('mview_%d',fg.Number)); 
        
        switch evt.Key
            
            case 'h' % help (display interactive commands)
                fprintf("Interactive commands (keyboard)\n\n");
                fprintf("  'h' : list interactive commands\n"); 
                fprintf("  'f' : go forward (show next frame)\n"); 
                fprintf("  'b' : go backward (show previous frame)\n"); 
                fprintf("  'r' : redraw the figure\n");
                fprintf("  'q' : quit (close the figure)\n"); 
                fprintf("\n"); 
               
            case 'q' % quit
                close(fg);
                
            case 'f' % forward 
                id = 1+mod(fg.UserData.frame_id,fg.UserData.nb_frames);
                im_hdl.CData = fg.UserData.contrast_slope*im_hdl.UserData(:,:,id) + fg.UserData.contrast_offset; 
                fg.UserData.frame_id = id; 
                fg.Name = sprintf("%s (frame #%d/%d)",fg.UserData.name,id,fg.UserData.nb_frames); 
                
            case 'b' % backward
                id = 1+mod(fg.UserData.frame_id-2,fg.UserData.nb_frames);
                im_hdl.CData = fg.UserData.contrast_slope*im_hdl.UserData(:,:,id) + fg.UserData.contrast_offset; 
                fg.UserData.frame_id = id; 
                fg.Name = sprintf("%s (frame #%d/%d)",fg.UserData.name,id,fg.UserData.nb_frames); 
                
            case 'r' % redraw
                disp('redraw figure');
                drawmovie(u,m,M); 
        end
    end

    function fg = drawmovie(u,m,M)
        gr = groot();
        u = double(u);
        [ny,nx,nimages] = size(u);
        fg = gcf();
        a = 255./(M-m); b = -m*a;
        if(isnan(a)); a = 1; b = 0; end
        image(a*u(:,:,1)+b,'CDataMapping','direct','XData',0:nx-1,'YData',0:ny-1,'UserData',u,'Tag',sprintf('mview_%d',fg.Number));
        axis 'image';
        fg.UserData = struct('contrast_slope',a,'contrast_offset',b,'nb_frames',nimages,'frame_id',1,'name',fg.Name);
        fg.Name = sprintf("%s (frame #1/%d)",fg.Name,nimages);
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
disp("mview (frame-by-frame movie displayer): press 'f' to go forward, 'b' to go backward and 'h' to display the other interactive commands.");
if z == 1
    fg = drawmovie(u,m,M); 
else
    [ny,nx,nim] = size(u);
    v = reshape(kron(reshape(u,[ny,nx*nim]),ones(z)),[z*ny,z*nx,nim]);
    fg = drawmovie(v,m,M);
end
end
