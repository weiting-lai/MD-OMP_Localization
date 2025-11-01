function [vs_xyz, vs_count, vsx, vsy, vsz] = generate_grid(xspec, yspec, zspec)
% Function  : generate_grid
% Creation  : 24-07-25
% Author    : Weiting Lai
% Version   : -
% -------------------------------------------------------------
% GENERATE_GRID  Generate 3D virtual source points.
%
% [vs_xyz, vs_count, vsx, vsy, vsz] = pos.generate_grid(xspec, yspec, zspec)
%
% Input:
%   xspec, yspec, zspec - Each can be one of the following:
%       1) A 3-element vector [min max N] → generates linspace(min, max, N)
%       2) A direct coordinate vector, e.g., linspace(-1.4, 1.4, 30)
%
% Output:
%   vs_xyz   : (N × 3) array of [x y z] coordinates
%   vs_count : Total number of virtual source points
%   vsx, vsy, vsz : Coordinate vectors used for each axis
%
% Example:
%   [points, count] = pos.generate_grid([-1.4 1.4 30], [-2.0 2.0 45], [-0.2 0.2 5]);

    vsx = parse_spec(xspec, 'xspec');
    vsy = parse_spec(yspec, 'yspec');
    vsz = parse_spec(zspec, 'zspec');

    % Use ndgrid for efficient 3D grid generation
    [X, Y, Z] = ndgrid(vsx, vsy, vsz);
    vs_xyz = [X(:), Y(:), Z(:)];
    vs_count = size(vs_xyz, 1);
end

% ---- Local helper -------------------------------------------------------
function v = parse_spec(spec, name)
%PARSE_SPEC  Helper to interpret input as either [min max N] or a direct vector.

    validateattributes(spec, {'numeric'}, {'vector','nonempty'}, mfilename, name);
    spec = spec(:).';  % Ensure row vector
    if numel(spec) == 3
        n = max(1, round(spec(3)));
        v = linspace(spec(1), spec(2), n);
    else
        v = spec; % Use directly as coordinate vector
    end
end
