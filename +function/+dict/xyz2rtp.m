function [r, t, p] = xyz2rtp(x, y, z)
% Package   : libmatlabpackage
% Function  : xyz2rtp
% Creation  : 27-01-21
% Author    : Lachlan Birnie
% Version   : 02-03-23
%
% Description
%
%   Convert cartesian coordinates to spherical coordinates:
%   (radius, theta, phi)
%   (radius, elevation, azimuth)
%
% Inputs
%
%   xyz     : [Q by 3] (x,y,z) position of each point (Q points total).
%   - or -    (input type depends on number of inputs).
%   x,y,z   : [Q by 1] vector of each points x,y,z position.
%
% Outputs
%
%   rtp     : [Q by 3] (r,t,p) position of each point (Q points total).
%             (radius, elevation 0 to pi, azimuth 0 to 2pi).
%   - or -    (output type depends on number of outputs).
%   r,t,p   : [Q by 1] vecotr of each points r, t, and p position.

    if (nargin == 1)
        % Split matrix input into vectors.
        temp = x;
        x = temp(:,1);
        y = temp(:,2);
        z = temp(:,3);
    end
    
    % Make sure positions are columns.
    if size(x,2) > size(x,1), x = x.'; end
    if size(y,2) > size(y,1), y = y.'; end
    if size(z,2) > size(z,1), z = z.'; end
    
    r = sqrt(x.^2 + y.^2 + z.^2);
    t = (pi/2) - atan2(z, sqrt(x.^2 + y.^2));
    p = wrapTo2Pi(atan2(y, x));
    
    if (nargout == 1)
        % Combine outputs into matrix.
        r = [r, t, p];
    end
    
end