function [g] = greens_matrix(source_type, k, source_xyz, mic_xyz)
% Package   : libmatlabpackage.sph
% Function  : greens_matrix
% Creation  : 14-03-23
% Author    : Lachlan Birnie
% Version   : -
%
% Description:
%
%   Returns a matrix of greens_function terms as elements.
%   Matrix shape is [mics by sources (by k)]

    % import xyz2rtp
    % import rtp2xyz

    source_type = lower(source_type);

    switch source_type

        case {'pw', 'planewave', 'plane_wave', 'plane-wave', ...
              'pw-', 'planewave-', 'plane_wave-', 'plane-wave-', ...
              'pw+', 'planewave+', 'plane_wave+', 'plane-wave+'}

            % Incident angle of planewaves.
            [~, pw_theta, pw_phi] = xyz2rtp(source_xyz);
            pw_ky = rtp2xyz(k, pw_theta, pw_phi);

            % Distance term, d = ky dot x.
            pw_d = mic_xyz * pw_ky.';  % [mic by source]

            % Far-field greens.
            switch source_type
                case {'pw+', 'planewave+', 'plane_wave+', 'plane-wave+'}
                    % Incoming (inward / sink) planewave.
                    g = exp((+1i) .* pw_d);
                otherwise
                    % Outgoing (outward) planewave.
                    g = exp((-1i) .* pw_d);
            end

        % TODO: case {'pointsource'}
        case {'ps', 'pointsource', 'point_source', 'point-source', ...
              'ps+', 'pointsource+', 'point_source+', 'point-source+', ...
              'ps-', 'pointsource-', 'point_source-', 'point-source-'}
          
            % Distance of source to receiver ||x-y||.
            pw_d = sqrt((source_xyz(:,1).' - mic_xyz(:,1)).^2 ...
                    +(source_xyz(:,2).' - mic_xyz(:,2)).^2 ...
                    +(source_xyz(:,3).' - mic_xyz(:,3)).^2 ...
                    );  % [Q,L]
            
            % Point source Green's function.
            switch source_type
                case {'ps', 'pointsource', 'point_source', 'point-source', ...
                      'ps+', 'pointsource+', 'point_source+', 'point-source+'}
                    g = exp((+1i) .* k .* pw_d) ./ (4 .* pi .* pw_d);  % [Q,L]
                otherwise
                    g = exp((-1i) .* k .* pw_d) ./ (4 .* pi .* pw_d);  % [Q,L]
            end
            


        % TODO: case {'mixedwave'}
            
        otherwise
            temp = dpstack;
            error('%s(): failed to determine source_type:%s', temp(1).name, source_type);
    end

end