

I_neighbor=find(firings(1:end,2)==repr+1);                     % index of a first nighbor, once fired
II_neighbor=find(firings(1:end,2)==repr+2);                    % index of a second neighbor, once fired

wave_delay=(firings(II_neighbor(1)) - firings(I_neighbor(1)))*dt;         % wave delay in steps

wave_speed=1/wave_delay;                                            % speed in arbitrary units

% just to save the current speed
% wave_speed=wave_speed;
% wave_delay_8=wave_delay;