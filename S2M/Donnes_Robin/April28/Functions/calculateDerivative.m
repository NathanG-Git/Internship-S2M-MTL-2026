function derivativeData = calculateDerivative(data, frameRate)
    % Validate input data
    if ~ismatrix(data) || isempty(data)
        error('Input data must be a non-empty 2D matrix.');
    end

    % Validate frame rate
    if ~isscalar(frameRate) || frameRate <= 0
        error('Frame rate must be a positive scalar value.');
    end

    % Number of frames and axes in the data
    [nFrames, nAxes] = size(data);

    % Validate if there are enough frames to calculate the derivative
    if nFrames < 2
        error('Data must contain at least two time frames to calculate the derivative.');
    end

    % Initialize the derivative data matrix
    derivativeData = zeros(nFrames, nAxes);

    % Calculate the time interval for the central difference
    timeInterval = 2 / frameRate;  % The interval covers two frames

    % Loop over each axis
    for iFr = 2:nFrames-1
        for iAx = 1:nAxes
            % Central difference for interior points
            derivativeData(iFr, iAx) = (data(iFr + 1, iAx) - data(iFr - 1, iAx)) / timeInterval;
        end
    end

    % Forward difference for the first data point
    derivativeData(1, :) = (data(2, :) - data(1, :)) / (1 / frameRate);

    % Backward difference for the last data point
    derivativeData(end, :) = (data(end, :) - data(end - 1, :)) / (1 / frameRate);
end



