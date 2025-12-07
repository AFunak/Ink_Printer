%% Continuous Inkjet Printer (CIJ) Deflection Simulation

%% -------- Physical Parameters --------
Vx = 2;                 % horizontal velocity (m/s)
Vmax = 200;              % maximum voltage applied (V)
W  = 1e-3;               % capacitor plate spacing (m)
L1 = (5e-3);           % capacitor length (m)
L2 = 1.25e-2;            % distance from capacitor to paper (m)

% Droplet Diameter
droplet_diameter_um = 84; % Droplet diameter in micrometers (µm)
droplet_diameter    = droplet_diameter_um * 1e-6; % Droplet diameter in meters (m)


%% Droplet Mass (m) and Charge (q) Calculation
% NOTE: We assume the fluid is water (or similar) with density rho_ink.
% The original parameters (m_orig, q_orig) are used to derive these constants
% for a default diameter (D_orig).

m_orig = 2.7e-7;             % original droplet mass (kg)
q_orig = (-1.9e-10);           % original droplet charge (C)
D_orig = 50.1e-6;            % Implied diameter based on original mass (assuming rho_ink ~1000 kg/m^3)

% 1. Recalculate Mass (m)
% Mass is proportional to Volume (D^3).
m = m_orig * (droplet_diameter / D_orig)^3;

% 2. Recalculate Charge (q)
% We assume charge q is proportional to D^2 for a "smaller charge density"
% effect, which means larger droplets require proportionally more voltage.
q = q_orig * (droplet_diameter / D_orig)^2;

fprintf('--- Droplet Parameters ---\n');
fprintf('Droplet Diameter: %.2f µm\n', droplet_diameter_um);
fprintf('Droplet Mass (m): %.2e kg\n', m);
fprintf('Droplet Charge (q): %.2e C\n', q);
fprintf('--------------------------\n');


%% -------- Timing and Geometry --------
dpi = 100;              % dots per inch
inch = 0.0254;          % meters in 1 inch
line_height = 0.01;    % 10 mm tall letter
num_dots = round(dpi * (line_height / inch));

T = L1 / Vx;            % time inside capacitor
t_flight = L2 / Vx;     % time from capacitor to paper
x_paper = L1 + L2;      % horizontal location of paper

%% -------- Droplet motion with no voltage --------
t_no_voltage = x_paper / Vx;
fprintf('Time for droplet to reach paper with no voltage: %.6f s\n', t_no_voltage);

%% -------- Vertical target positions --------
% Calculate the required vertical position on the paper
y_target = linspace(-line_height/2, line_height/2, num_dots);

%% -------- Compute Voltages (CORE PHYSICS) --------
% 1. Find the necessary vertical velocity (Vy) at the capacitor exit (L1)
%    y_target = Vy * t_flight  =>  Vy = y_target / t_flight
Vy = y_target / t_flight;

% 2. Find the necessary constant vertical acceleration (ay) in the capacitor
%    Vy = ay * T              =>  ay = Vy / T
ay = Vy / T;

% 3. Find the required Electric Field (E) based on F=ma and F=qE (F_y = m*ay = q*E)
%    E = (m * ay) / q
E  = (m .* ay) ./ q;

% 4. Find the required Voltage (V) based on E=V/W (V = E * W)
V  = E .* W;

% Clamp voltages to ±Vmax
% This scales the required voltages to fit within the maximum allowed voltage.
% Note: If max(abs(V)) is > Vmax, the print height will be reduced.
scale_factor = Vmax / max(abs(V));
V = V * scale_factor;

fprintf('Maximum required voltage (unclamped): %.1f V\n', max(abs(E)).*W);
fprintf('Actual maximum applied voltage: %.1f V\n', max(abs(V)));


%% -------- Compute time stamps --------
% Assuming droplets are released at equal time intervals to hit the paper sequentially
time_stamps = (0:num_dots-1) * (T + t_flight);
total_time = time_stamps(end);
fprintf('Total time to print %d dots: %.6f s\n', num_dots, total_time);


%% -------- Plot staircase voltage --------
figure('Name','Voltage Staircase','Color','w','Position',[100 100 800 400]);
stairs(time_stamps, V, 'LineWidth', 2, 'Color', 'b');
hold on;
plot(time_stamps, V, 'ro', 'MarkerFaceColor','r');
xlabel('Time (s)'); ylabel('Voltage V(t) [V]');
title(sprintf('Staircase Voltage Profile for Letter "I" (D=%.2fµm)', droplet_diameter_um));
grid on;
ylim([-Vmax-20, Vmax+20]);
legend('Voltage Profile','Droplet Times','Location','best');
drawnow;

%% -------- Animation Settings --------
speed_scale = 10;
Vx_slow = Vx / speed_scale;
pause_time = 0.00019;
dt = 0.00029;
animation_enabled = true;

    %% -------- Animation Setup --------
    if animation_enabled
    
    % --- VISUAL MARKER SIZE ADJUSTMENT ---
    % Use a fixed, small marker size (4) for visual clarity in the animation,
    % regardless of the physical droplet diameter (84 um).
    visual_marker_size = 4; 
    
    fig_anim = figure('Name','Droplet Animation','Color','w','Position',[900 100 800 600]);
    hold on; axis equal;
    xlabel('x (m)'); ylabel('y (m)');
    xlim([-0.002, L1 + L2 + 0.003]);
    ylim([-line_height, line_height]*2);  % scale to show paper and capacitor
    title(sprintf('Droplet Traveling Through Capacitor - Letter "I" (D=%.2fµm)', droplet_diameter_um));

    %% Draw capacitor
    x_offset = 0.00025;
    cap_length_scale = 0.5;
    L1_draw = L1 * cap_length_scale;

    cap_x = [-0.001 + x_offset, -0.001 + L1_draw + x_offset, ...
             -0.001 + L1_draw + x_offset, -0.001 + x_offset];
    cap_y = [-W/2, -W/2, W/2, W/2];

    fill(cap_x, cap_y, [0.8 0.8 0.8], 'FaceAlpha',0.3, 'EdgeColor','k','LineWidth',1.5);
    line([-0.001 + x_offset, -0.001 + L1_draw + x_offset], [W/2 W/2], 'Color','r','LineWidth',2);
    line([-0.001 + x_offset, -0.001 + L1_draw + x_offset], [-W/2 -W/2], 'Color','b','LineWidth',2);
    text(-0.001 + L1_draw/2 + x_offset, W/2 + 0.001, '+ Plate','Color','r','HorizontalAlignment','center');
    text(-0.001 + L1_draw/2 + x_offset, -W/2 - 0.001, '– Plate','Color','b','HorizontalAlignment','center');

    %% Draw paper
    paper_x = -0.001 + L1 + L2;
    paper_width = 0.002;
    paper_height = 0.05;  % slightly larger for visibility
    px = [paper_x, paper_x+paper_width, paper_x+paper_width, paper_x];
    py = [-paper_height/2, -paper_height/2, paper_height/2, paper_height/2];
    fill(px, py, [1,1,0.8], 'EdgeColor','k','FaceAlpha',0.3,'LineWidth',1.5);
    text(paper_x + paper_width/2, paper_height/2 + 0.001, 'Paper', ...
         'HorizontalAlignment','center','FontSize',10,'Color',[0.3 0.3 0.3]);

    %% Pre-allocate final dot positions
    dot_positions = zeros(num_dots, 2);

    %% -------- Droplet Animation --------
    fprintf('Animating %d droplets...\n', num_dots);

    for k = 1:num_dots
    % Initialize droplet
    % Use the fixed, small size for the moving marker
    droplet = plot(-0.001, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', visual_marker_size);
    trajectory = plot(NaN, NaN, 'k-', 'LineWidth', 0.5);

    % -------- Initialize position and velocity --------
    x = -0.001 - 0.00175;   % Start 1.75 mm farther back
    y = 0;
    Vy_inst = 0;

    % Determine acceleration based on voltage V(k), mass (m), and charge (q)
    % E = V/W, F = qE, a = F/m = qV / (mW)
    ay_now = (q * V(k)) / (m * W);

    hit_capacitor = false;

    while x < (L1 + L2)
        % -------- Acceleration only inside capacitor --------
        a_now = (x >= 0 && x <= L1) * ay_now;

        % Update velocity & position
        Vy_inst = Vy_inst + a_now*dt;
        y = y + Vy_inst*dt;
        x = x + Vx_slow*dt;

        % Check for capacitor collision
        if (x >= 0 && x <= (L1-.0031)) && (y >= W/2 || y <= -W/2)
            hit_capacitor = true;
            break;
        end

        % Update droplet marker
        if ishandle(droplet)
            set(droplet,'XData',x,'YData',y);
        end
        drawnow limitrate;
        pause(pause_time);
    end

    % Record final position if droplet didn't hit plates
    if ~hit_capacitor
        dot_positions(k,:) = [x, y];
        % Use a small size for the impact point
        plot(x, y, 'bo', 'MarkerFaceColor','b','MarkerSize', visual_marker_size * 0.75);
    end

    % Remove moving droplet
    if ishandle(droplet); delete(droplet); end
    if ishandle(trajectory); delete(trajectory); end
end

    % Plot all final dots (Use a slightly larger marker for the final printed pattern)
    plot(dot_positions(:,1), dot_positions(:,2), 'b.', 'MarkerSize', visual_marker_size * 1.5);
    title(sprintf('Completed Letter "I" (D=%.2fµm)', droplet_diameter_um));
else
    fprintf('Animation skipped. Voltage plot complete.\n');
end

%% -------- Final Summary --------
fprintf('\n=== SIMULATION RESULTS ===\n');
fprintf('Droplet Diameter: %.2f µm\n', droplet_diameter_um);
fprintf('Droplet Mass (m): %.2e kg\n', m);
fprintf('Droplet Charge (q): %.2e C\n', q);
fprintf('Number of droplets: %d\n', num_dots);
fprintf('Total print time: %.6f s\n', total_time);
fprintf('Voltage range: %.1f V to %.1f V\n', min(V), max(V));
fprintf('Vertical span: %.1f mm\n', (max(y_target)-min(y_target))*1000);
disp('Simulation complete!');