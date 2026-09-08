function result = custom_integral(j, theta, t)
    % Parameters:
    %   j        : Array representing j(t); in rad/s^3
    %   theta    : Array representing theta(t); in radians
    %   t        : Time points corresponding to j and theta; in seconds
    
    % Returns:
    %   result   : Result of the integral

    % Compute the 3rd derivative of theta with respect to t
   % d3_theta = diff(theta, 3) ./ diff(t, 3);

    % Adjust size of j and t arrays to match size of d3_theta
    j = j(1:end-3);
    t = t(1:end-3);
    theta = theta(1:end-3);

    % Define the integrand
    integrand = (j.^2) .* (t.^5 ./ theta.^2);

    % Compute the integral using the trapezoidal rule for discrete data
    integral_val = trapz(t, integrand);

    % Multiply by the factor outside the integral
    result = sqrt(1/2) * sqrt(integral_val);
end
