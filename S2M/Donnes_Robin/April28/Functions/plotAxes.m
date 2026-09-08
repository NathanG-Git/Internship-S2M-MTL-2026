function h = plotAxes(varargin)
% Show system of axes (a 4x4 column matrix of rototranslation RT or 3x3 matrix of rotation)
% where RT = [Xx Yx Zx Tx
%             Xy Yy Zy Ty
%             Xz Yz Zz Tz
%             0  0  0  1 ] 
% If a 3x3 is sent, T = [0 0 0]'.
%
% SYNTAX
%  h = plotAxes(RT) plots a new system of axes, if a graph already exists,
%    it turns on hold on and plot on it
%  h = plotAxes(RT, h) updates to a new RT
%  h = plotAxes(h, RT) updates to a new RT
%  h = plotAxes(h, options) changes options ('param', value), keep same RT values 
%  h = plotAxes(RT, options) creates a new plot and changes options ('param', value) see Options section.
%  h = plotAxes(RT, h, options) updates to a new RT and changes options ('param', value)
%  h = plotAxes(h, RT, options) updates to a new RT and changes options ('param', value)
%  h = plotAxes() create new plot, shows nothing
%  h = plotAxes(options) create new plot options ('param', value), shows nothing
%
% OPTIONS - options are ('optionName', value), [] are default values 
%   'length'      = [0.05];        ==> length of axes (default optimize for
%                                      data expressed in meters.
%   'LineWidth'   = [3];           ==> width of axes
%   'Color'       = ['rgb'];       ==> Color of x, y and z axes respectively,
%                                      must b three letters
%   'arrow'       = [false], true; ==> Show arrows instead of axes
%   'MaxHeadSize' = [8];           ==> If arrow is true, this changes the
%                                      size of the arrows' head
%
% EXAMPLE
%   RT = eye(4);                        % Let define a new RT
%   h = plotAxes(RT);                   % Create a new plot and get the handle
%   RT(1:3,4) = 0.2;                    % Let moves this RT
%   plotAxes(h, RT, 'color', 'rrr');    % Update its position and change axes color to red 
%   RT2 = eye(4);                       % Let create a new RT
%   h2 = plotAxes('color', 'bbb');      % Prepare a new plot with blue axes
%   plotAxes(RT2, h2);                  % Plot it
%
% DISCLAIMER and COPYRIGHT
%   This function was written by 
%                             Pariterre
%                             laboratoire S2M
%                             Université de Montréal
%   Version 1 : April 13th, 2016


    % Arguments identification
    changeRT = true; % If this flag is true, RT is redrawn
    
    % if the first argument is the RT matrix
    if ~isempty(varargin) && isnumeric(varargin{1})
        % Check the second argument if is the handler, if
        if length(varargin) > 1 && isstruct(varargin{2})
            % Get the handler and set the parser
            h = varargin{2};
            parseur = 3:2:length(varargin);
        else
            % Create new handler and set the parser
            h = setDefault();
            parseur = 2:2:length(varargin);
        end
        % Replace RT with the new one
        RT = varargin{1};
        
    % if the first argument is the handler
    elseif ~isempty(varargin) && isstruct(varargin{1})
        h = varargin{1};
        % Check the second argument if is the handler, if
        if length(varargin) > 1 && isnumeric(varargin{2})
            % Get the new RT and set the parser
            RT = varargin{2};
            parseur = 3:2:length(varargin);
        else
            % Create new handler and set the parser
            changeRT = false;
            parseur = 2:2:length(varargin);
        end
        
    % If nothing was sent or just options
    else 
        % Prepare a new structure with default values and prepare parser
        RT = nan(4);      % Empty RT
        h = setDefault();
        parseur = 1:2:length(varargin);
    end

    % Parsing arguments 
    somethingChanged = false;
    for i = parseur
        if ~ischar(varargin{i})
            error('(''Paramètre'',''Valeur'' doit être respecté');
        end
        if length(varargin) < i+1
            error('La valeur pour %s n''est pas présente', varargin{i});
        end

        if strcmpi(varargin{i}, 'arrow')
            h.Arrow = varargin{i+1};
            somethingChanged = true;
        elseif strcmpi(varargin{i}, 'length')
            h.Length = varargin{i+1};
            somethingChanged = true;
        elseif strcmpi(varargin{i}, 'linewidth')
            h.LineWidth = varargin{i+1};
            somethingChanged = true;
        elseif strcmpi(varargin{i}, 'color')
            h.Color = varargin{i+1};
            if length(h.Color) == 1
                h.Color = [h.Color h.Color h.Color];
            end
            somethingChanged = true;
        elseif strcmpi(varargin{i}, 'maheadsize')
            h.MaxHeadSize = varargin{i+1};
            somethingChanged = true;
        else
            error('%s n''est pas une propriété', varargin{i});
        end
    end
    
    if isempty(h.hs)
        hold on
        if h.Arrow
            h.hs(1) = quiver3(nan,nan,nan,nan,nan,nan);
            h.hs(2) = quiver3(nan,nan,nan,nan,nan,nan);
            h.hs(3) = quiver3(nan,nan,nan,nan,nan,nan);
        else
            h.hs(1) = plot3(nan,nan,nan);
            h.hs(2) = plot3(nan,nan,nan);
            h.hs(3) = plot3(nan,nan,nan);
        end
        somethingChanged = true;
    end

    if somethingChanged
        if h.Arrow
            set(h.hs(1), 'LineWidth', h.LineWidth, 'color', h.Color(1), 'maxheadsize', h.MaxHeadSize);
            set(h.hs(2), 'LineWidth', h.LineWidth, 'color', h.Color(2), 'maxheadsize', h.MaxHeadSize);
            set(h.hs(3), 'LineWidth', h.LineWidth, 'color', h.Color(3), 'maxheadsize', h.MaxHeadSize);
        else
            set(h.hs(1), 'LineWidth', h.LineWidth, 'color', h.Color(1));
            set(h.hs(2), 'LineWidth', h.LineWidth, 'color', h.Color(2));
            set(h.hs(3), 'LineWidth', h.LineWidth, 'color', h.Color(3));
        end
    end
    

    if changeRT
        % Check RT size
        if size(RT, 1) == 3 || size(RT, 2) == 3
            RT(1:3,4) = 0;
            RT(4,:) = [0 0 0 1];
        elseif size(RT, 1) ~= 4 || size(RT, 2) ~= 4
            error('RT must be a 4x4 homogeneous matrix or a 3x3 rotation matrix');
        end
        
        % Convert RT into axes
        X = RT * [0 0 0 1;h.Length 0      0      1]';
        Y = RT * [0 0 0 1;0      h.Length 0      1]';
        Z = RT * [0 0 0 1;0      0      h.Length 1]';

        % Plot them
        if h.Arrow
            set(h.hs(1), 'Xdata', X(1,1), 'Ydata', X(2,1), 'Zdata', X(3,1), 'Udata', X(1,2)-X(1,1), 'Vdata', X(2,2)-X(2,1), 'Wdata', X(3,2)-X(3,1));
            set(h.hs(2), 'Xdata', Y(1,1), 'Ydata', Y(2,1), 'Zdata', Y(3,1), 'Udata', Y(1,2)-Y(1,1), 'Vdata', Y(2,2)-Y(2,1), 'Wdata', Y(3,2)-Y(3,1));
            set(h.hs(3), 'Xdata', Z(1,1), 'Ydata', Z(2,1), 'Zdata', Z(3,1), 'Udata', Z(1,2)-Z(1,1), 'Vdata', Z(2,2)-Z(2,1), 'Wdata', Z(3,2)-Z(3,1));
        else
            set(h.hs(1), 'Xdata', X(1,:), 'Ydata', X(2,:), 'Zdata', X(3,:));
            set(h.hs(2), 'Xdata', Y(1,:), 'Ydata', Y(2,:), 'Zdata', Y(3,:));
            set(h.hs(3), 'Xdata', Z(1,:), 'Ydata', Z(2,:), 'Zdata', Z(3,:));
        end
    end

    if nargout == 0
        clear h
    end
end

function h = setDefault()
    % Set default param
    h.hs = [];          % Handler on plots
    h.Length = 0.05;    % Length of the lines
    h.LineWidth = 3;    % Width of the lines
    h.Color = 'rgb';    % Colors of line (X,Y,Z respectively)
    h.Arrow = false;    % Is this an arrow
    h.MaxHeadSize = 8;  % Head size of the arrow
end