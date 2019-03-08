function rh = fun_rel_hum(theta, x, pres)
    % to calculate the relative humidity from the absolute humidity
    % 1 ... temperature in °C
    % 2 ... absolute humidity
    % 3 ... pressure

    ps = exp(23.462-(3978.205 ./ (233.349+theta)) ); %% for theta>0 based on Bauphysik presentations
    rh = x./(x+0.622) .* (pres./ps);
end