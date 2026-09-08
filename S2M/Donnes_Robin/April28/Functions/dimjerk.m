function result = dimjerk(j, Vmean, t)

    % Compute the integral
    integralValue = trapz(t, j.^2);

    D = t(end) - t(1);

    % Compute the final result
    result = (D^3 * integralValue) / (Vmean^2);

end
