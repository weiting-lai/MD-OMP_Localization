function plot_room(mic_xyz, true_src, est_src, room_dims, cfg)
% Function  : plot_room
% Creation  : 24-07-25
% Author    : Weiting Lai
% Version   : -
% -------------------------------------------------------------
% PLOT_ROOM  3D room view: microphones (red), sources (blue), estimate (green).
%
% plt.plot_room(mic_xyz, true_src, est_src, room_dims, cfg)
%
% Inputs
%   mic_xyz   : M×3 microphone positions
%   true_src  : S×3 ground-truth source positions
%   est_src   : K×3 estimated positions (for the selected method); [] to skip
%   room_dims : [L W H] in meters (x, y, z)
%   cfg       : plotting config (all optional fields)
%                 .save_path = '' | 'fig_room.png'
%                 .view      = [-35 20]     (az, el)
%                 .legend_loc= 'best'

   % ---- defaults --------------------------------------------------------
    if nargin < 5 || isempty(cfg), cfg = struct; end
    if ~isfield(cfg,'save_path'), cfg.save_path = ''; end
    if ~isfield(cfg,'view'),      cfg.view      = [-35 20]; end
    if ~isfield(cfg,'legend_loc'),cfg.legend_loc= 'best'; end

    L = room_dims(1); W = room_dims(2); H = room_dims(3);

    % ---- room patch ------------------------------------------------------
    V = [ 0 0 0;  L 0 0;  L W 0;  0 W 0; ...
          0 0 H;  L 0 H;  L W H;  0 W H ];
    F = [ 1 2 3 4;  5 6 7 8;  1 2 6 5;  2 3 7 6;  3 4 8 7;  4 1 5 8 ];

    figure('Color','w', 'Position', [100 100 1000 800]);
    patch('Vertices',V,'Faces',F,'FaceColor','none','EdgeColor','k','LineWidth',1.5);
    hold on;

    % ---- microphones (red) & true sources (blue) -------------------------
    h_mic = plot3(mic_xyz(:,1), mic_xyz(:,2), mic_xyz(:,3), ...
        'k*', 'MarkerSize',10, 'MarkerFaceColor','k');
    h_src = plot3(true_src(:,1), true_src(:,2), true_src(:,3), ...
        'bsquare', 'MarkerSize',10, 'MarkerFaceColor','b');

    % ---- estimated points (selected method) ------------------------------
    if ~isempty(est_src)
        h_est = plot3(est_src(:,1), est_src(:,2), est_src(:,3), ...
            'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
        legend([h_mic, h_src, h_est], {'Microphones','Sources','Estimated'}, ...
               'Location', cfg.legend_loc);
    else
        legend([h_mic, h_src], {'Microphones','Sources'}, 'Location', cfg.legend_loc);
    end

    xlabel('$x$ [m]','Interpreter','latex');
    ylabel('$y$ [m]','Interpreter','latex');
    zlabel('$z$ [m]','Interpreter','latex');

    set(gca,'FontName','Times New Roman','FontSize',28);
    xlim([0 L]); ylim([0 W]); zlim([0 H]);
    axis equal;
    axis([0 L 0 W 0 H]);
    grid on;
    view(cfg.view);
    hold off;

    % ---- save (optional) -------------------------------------------------
    if ~isempty(cfg.save_path)
        exportgraphics(gca, cfg.save_path, 'Resolution', 300);
    end
end