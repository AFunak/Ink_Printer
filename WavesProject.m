%% Inkjet Printer Simulation - Letter "I" 5mm, ±200V
clear; clc; close all;

%% -------- Physical Parameters --------
Vx = 20;                 % horizontal velocity (m/s)
m  = 2.7e-7;             % droplet mass (kg)
q  = -1.9e-10;             % droplet charge (C)
W  = 1e-3;              % capacitor plate spacing (m)
L1 = 5e-3;              % capacitor length (m) .5 mm
L2 = 1.25e-2;              % distance from capacitor to paper (m) 1.25mm
Vmax = 200;             % maximum voltage applied (V)

dpi = 100;              % dots per inch
inch = 0.0254;          % meters in 1 inch
line_height = 0.01;    % 5 mm tall letter
num_dots = round(dpi * (line_height / inch));

T = L1 / Vx;            % time inside capacitor
t_flight = L2 / Vx;     % time from capacitor to paper
x_paper = L1 + L2;      % horizontal location of paper

%% -------- Droplet motion with no voltage --------
t_no_voltage = x_paper / Vx;
fprintf('Time for droplet to reach paper with no voltage: %.6f s\n', t_no_voltage);

%% -------- Vertical target positions --------
y_target = linspace(-line_height/2, line_height/2, num_dots);

%% -------- Compute Voltages --------
Vy = y_target / t_flight;      % vertical velocity needed at exit
ay = Vy / T;                    % acceleration in capacitor
E  = (m .* ay) ./ q;            % required electric field
V  = E .* W;                     % voltage required

% Clamp voltages to ±Vmax
V = V * Vmax / max(abs(V));

%% -------- Compute time stamps --------
time_stamps = (0:num_dots-1) * (T + t_flight);
total_time = time_stamps(end);
fprintf('Total time to draw letter "I": %.6f s\n', total_time);

%% -------- Plot staircase voltage --------
figure('Name','Voltage Staircase','Color','w','Position',[100 100 800 400]);
stairs(time_stamps, V, 'LineWidth', 2, 'Color', 'b');
hold on;
plot(time_stamps, V, 'ro', 'MarkerFaceColor','r');
xlabel('Time (s)'); ylabel('Voltage V(t) [V]');
title('Staircase Voltage Profile for Letter "I"');
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
    fig_anim = figure('Name','Droplet Animation','Color','w','Position',[900 100 800 600]);
    hold on; axis equal;
    xlabel('x (m)'); ylabel('y (m)');
    xlim([-0.002, L1 + L2 + 0.003]);
    ylim([-line_height, line_height]*2);  % scale to show paper and capacitor
    title('Droplet Traveling Through Capacitor - Letter "I"');

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
    droplet = plot(-0.001, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize',6);
    trajectory = plot(NaN, NaN, 'k-', 'LineWidth', 0.5);

    % -------- Initialize position and velocity --------
    x = -0.001 - 0.00175;   % Start 1.75 mm farther back
    y = 0;
    Vy_inst = 0;

    % Determine acceleration direction based on voltage sign
    if V(k) >= 0
        ay_now = abs(q * V(k)/W) / m;   % positive voltage: top to bottom
    else
        ay_now = -abs(q * V(k)/W) / m;  % negative voltage: bottom to top
    end

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
        plot(x, y, 'bo', 'MarkerFaceColor','b','MarkerSize',4);
    end

    % Remove moving droplet
    if ishandle(droplet); delete(droplet); end
    if ishandle(trajectory); delete(trajectory); end
end

    % Plot all final dots
    plot(dot_positions(:,1), dot_positions(:,2), 'b.', 'MarkerSize',10);
    title('Completed Letter "I"');
else
    fprintf('Animation skipped. Voltage plot complete.\n');
end

%% -------- Final Summary --------
fprintf('\n=== SIMULATION RESULTS ===\n');
fprintf('Number of droplets: %d\n', num_dots);
fprintf('Total print time: %.6f s\n', total_time);
fprintf('Voltage range: %.1f V to %.1f V\n', min(V), max(V));
fprintf('Vertical span: %.1f mm\n', (max(y_target)-min(y_target))*1000);
disp('Simulation complete!');