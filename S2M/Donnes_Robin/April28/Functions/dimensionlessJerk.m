function Dim_Jerk = dimensionlessJerk(Jerk_values, Velocity_values, Range_of_motion, time_length)
% This function calculates the dimensionless jerk using the integrated square root method

% Notes:
% Jerk_values input matrix must be of length based on time and columns based on the number of joints
% Velocity_values input matrix must be of the corresponding data set to the jerk values. If it is not, the calculation will be incorrect. 
% Range_of_motion input matrix must be of the corresponding data set to the jerk and velocity values. If they are not, the calculation will be incorrect. 
% Since this calculation method controls for the amount of time, the number of frames/data points is essential for this calculation

% Input sizes validation
if length(Range_of_motion) ~= length(time_length)
    error('ROM and time_length arrays must have matching sizes.');
end

frames    = length(time_length);
numValues = size(Jerk_values, 2);
Dim_Jerk  = zeros(frames, numValues);

for iF = 1:frames
    for iJ = 1:numValues
        jerk = Jerk_values(time_length(iF,1):time_length(iF,2), iJ);
        Integrated_Jerk = trapz(jerk.^2);
        
        Vel_peak_squared = max(Velocity_values(time_length(iF,1):time_length(iF,2), iJ))^2;

        if Vel_peak_squared == 0
            warning('Velocity peak squared is zero. Skipping current joint computation.');
            continue;
        end

        Ang_Disp_squared = Range_of_motion(iF)^2;

        frequency = length(time_length(iF,1):time_length(iF,2));
        Denominator = (frequency^5) / (Vel_peak_squared * Ang_Disp_squared);

        Dim_Jerk(iF, iJ) = sqrt(0.5 * Denominator * Integrated_Jerk);
    end
end

end
