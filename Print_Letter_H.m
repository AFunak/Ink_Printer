%% CIJ Deflection Simulation - Straight "H" Legs, Wider Paper, Purple Filled Region
clearvars; close all; clc;

%% -------- Physical Parameters --------
Vx = 20;                 
Vmax = 200;              
W  = 5e-3;               
L1 = 5e-3;               
L2 = 1.25e-2;            

%% Droplet properties
droplet_diameter_um = 84;               
droplet_diameter = droplet_diameter_um * 1e-6; 

m_orig = 2.7e-7;             
q_orig = -1.9e-10;       
D_orig = 50.1e-6;            

m = m_orig * (droplet_diameter / D_orig)^3;
q = q_orig * (droplet_diameter / D_orig)^2;  

fprintf('--- Droplet Parameters ---\n');
fprintf('Droplet Diameter: %.2f µm\n', droplet_diameter_um);
fprintf('Droplet Mass (m): %.2e kg\n', m);
fprintf('Droplet Charge (q): %.2e C\n', q);
fprintf('--------------------------\n\n');

%% -------- Print geometry & timing --------
dpi = 300;              
inch = 0.0254;          
line_height = 0.01;     
line_width  = 0.008;    

num_dots_total = round(dpi * (line_height / inch));
num_dots = floor(num_dots_total / 3) * 3;   
N_s = num_dots / 3;

T_flight_cap = L1 / Vx;      
T_flight_coast = L2 / Vx;    
x_paper_plane = L1 + L2; 


%% -------- Vertical and Horizontal Targets for "H" --------
y_max = line_height/2;

% Left vertical leg
N_leg = N_s;
y1 = linspace(-y_max, y_max, N_leg);
x1 = -line_width/2 * ones(1,N_leg);

% Middle horizontal crossbar
N_cross = N_s;
y2 = zeros(1,N_cross);
x2 = linspace(-line_width/2, line_width/2, N_cross);

% Right vertical leg
y3 = linspace(-y_max, y_max, N_leg);
x3 = line_width/2 * ones(1,N_leg);

% Combine full arrays
y_target = [y1, y2, y3];
x_target = [x1, x2, x3];

% Absolute horizontal positions on paper
paper_x_positions = x_paper_plane + x_target;
num_dots = length(y_target); 

%% -------- Compute required ay and voltages --------
denominator = (0.5 * T_flight_cap^2 + T_flight_cap * T_flight_coast);
ay = y_target / denominator;  

E_required = (m .* ay) ./ q;  
V_required = E_required * W;  

% Clamp voltages to ±Vmax
Vmax_actual = max(abs(V_required));
V = V_required;
if Vmax_actual > Vmax
    scale_factor = Vmax / Vmax_actual;
    V = V * scale_factor;
    fprintf('WARNING: Required voltage (%.1f V) exceeds Vmax (%.1f V). Scaling applied (%.2f%%)\n',...
        Vmax_actual, Vmax, scale_factor*100);
end
fprintf('Maximum required voltage (unclamped): %.1f V\n', Vmax_actual);
fprintf('Actual maximum applied voltage: %.1f V\n\n', max(abs(V)));

%% -------- Plot voltage staircase --------
time_stamps_release = (0:num_dots-1)/10000;  
figure('Name','Voltage Staircase','Color','w','Position',[100 100 600 300]);
stairs(time_stamps_release*1000, V, 'LineWidth', 2);
hold on;
plot(time_stamps_release*1000, V, 'ko', 'MarkerFaceColor','k');
xlabel('Droplet Release Time (ms)');
ylabel('Voltage V(t) [V]');
title('Staircase Voltage Profile (horizontal plates)');
grid on;
ylim([-Vmax-20, Vmax+20]);
%% -------- Horizontal (side) capacitors voltage staircase --------
% Define 3-step voltage: left, middle, right (example ±Vmax/2, 0, ±Vmax/2)
N_s = num_dots / 3;  % number of dots per step
Vx_side = [ones(1,N_s)*(-Vmax/2), ones(1,N_s)*(0), ones(1,N_s)*(Vmax/2)];

% Ensure length matches num_dots
Vx_side = Vx_side(1:num_dots);

% Time stamps (same as droplet release times)
time_stamps_release = (0:num_dots-1) * (1/10000);  % if f_drop = 10 kHz

% Plot staircase
figure('Name','Side Capacitors Voltage','Color','w','Position',[200 150 800 400]);
stairs(time_stamps_release*1000, Vx_side, 'LineWidth', 2);
hold on;
plot(time_stamps_release*1000, Vx_side, 'ko', 'MarkerFaceColor','k');
xlabel('Droplet Release Time (ms)');
ylabel('Voltage V(t) [V]');
title('Voltage Across Side Capacitors (3-step Staircase)');
grid on;
ylim([-Vmax-20, Vmax+20]);
legend('Applied Voltage','Droplet Times','Location','best');

%% -------- Simulation & Animation settings --------
dt = 0.00029;           
pause_time = 0.0005;    
animation_enabled = true;

%% -------- Animate droplets --------
if animation_enabled
    visual_marker_size = 3;   
    fig = figure('Name','Droplet Animation','Color','w','Position',[900 100 700 500]);
    hold on; axis equal;
    xlabel('x (m)'); ylabel('y (m)');
    xlim([-0.003, x_paper_plane + 0.006]);
    ylim([-line_height*1.0, line_height*1.0]);
    title(sprintf('Droplet Through Capacitor - Letter "H" (D=%.2f µm)', droplet_diameter_um));

    cap_draw_start = -0.001 + 0.00025;    
    cap_draw_end   = cap_draw_start + L1;

    % Horizontal capacitor plates (top & bottom)
    cap_y_top = [W/2, W/2, W/2+0.00025, W/2+0.00025];
    cap_y_bot = [-W/2-0.00025, -W/2-0.00025, -W/2, -W/2];
    cap_x_draw = [cap_draw_start, cap_draw_end, cap_draw_end, cap_draw_start];
    fill(cap_x_draw, cap_y_top, [0.9 0.9 0.9], 'FaceAlpha',0.4, 'EdgeColor','k');
    fill(cap_x_draw, cap_y_bot, [0.9 0.9 0.9], 'FaceAlpha',0.4, 'EdgeColor','k');
    line([cap_draw_start, cap_draw_end], [W/2, W/2], 'Color','r','LineWidth',2);
    line([cap_draw_start, cap_draw_end], [-W/2, -W/2], 'Color','b','LineWidth',2);

    % Purple region fully filling the area between the top and bottom plates
    purple_color = [0.6 0 0.8];  % purple
    purple_alpha = 0.3;          % semi-transparent
    fill([cap_draw_start, cap_draw_end, cap_draw_end, cap_draw_start], ...
         [-W/2, -W/2, W/2, W/2], ...
         purple_color, 'FaceAlpha', purple_alpha, 'EdgeColor','none');

    % Draw wider paper
    paper_center_vis = cap_draw_start + x_paper_plane; 
    paper_width = 0.012;   
    paper_height = 0.05;
    px = [paper_center_vis - paper_width/2, paper_center_vis + paper_width/2, ...
          paper_center_vis + paper_width/2, paper_center_vis - paper_width/2];
    py = [-paper_height/2, -paper_height/2, paper_height/2, paper_height/2];
    fill(px, py, [1 1 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'k');
    text(paper_center_vis, paper_height/2 + 0.002, 'PAPER', 'HorizontalAlignment','center','FontWeight','bold');

    dot_positions = NaN(num_dots, 2);   
    hit_flags = false(num_dots,1);

    fprintf('Animating %d droplets...\n', num_dots);

    for k = 1:num_dots
        droplet_vis = plot(cap_draw_start, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', visual_marker_size);

        x_phys = 0;
        y_phys = 0;
        Vy_inst = 0;
        ay_now = ay(k);
        hit_cap = false;

        x_phys_end = x_paper_plane + 0.002;
        while x_phys < x_phys_end
            in_cap = (x_phys >= 0) && (x_phys <= L1);
            a_now = in_cap * ay_now;

            Vy_inst = Vy_inst + a_now * dt;    
            y_phys = y_phys + Vy_inst * dt;    
            x_phys = x_phys + Vx * dt;         

            if in_cap && (abs(y_phys) >= W/2)
                hit_cap = true;
                break;
            end

            if ishandle(droplet_vis)
                set(droplet_vis, 'XData', x_phys + cap_draw_start, 'YData', y_phys);
            end

            drawnow limitrate;
            pause(pause_time);
        end

        if ~hit_cap
            x_land = paper_x_positions(k);   
            y_land = y_phys;
            dot_positions(k, :) = [x_land, y_land];
            plot(x_land + cap_draw_start, y_land, 'k.', 'MarkerSize', visual_marker_size * 1.5);
        else
            hit_flags(k) = true;
            plot(x_phys + cap_draw_start, y_phys, 'r.', 'MarkerSize', visual_marker_size * 1.5);
        end

        if ishandle(droplet_vis); delete(droplet_vis); end
    end

    % Final printed pattern
    figure('Name','Final Printed Pattern','Color','w','Position',[100 550 500 500]);
    hold on; axis equal; grid on;

    valid = ~isnan(dot_positions(:,1));
    x_mm = (dot_positions(valid,1) - x_paper_plane) * 1000;  
    y_mm = dot_positions(valid,2) * 1000;

    scale = 0.6; 
    plot(x_mm*scale, y_mm*scale, 'k.', 'MarkerSize', visual_marker_size*2);

    num_hits = sum(hit_flags);
    fprintf('Droplets that collided with capacitor plates: %d / %d\n', num_hits, num_dots);

    xlabel('Horizontal Position on Paper (mm)');
    ylabel('Vertical Deflection (mm)');
    title(sprintf('Final Printed "H" Pattern (D=%.2f µm)', droplet_diameter_um));
    xlim([-paper_width/2*1000*scale, paper_width/2*1000*scale]);
    ylim([-line_height*0.8, line_height*0.8]*1000*scale);
end

%% -------- Final summary --------
fprintf('\n=== SIMULATION RESULTS ===\n');
fprintf('Droplet Diameter: %.2f µm\n', droplet_diameter_um);
fprintf('Droplet Mass (m): %.2e kg\n', m);
fprintf('Droplet Charge (q): %.2e C\n', q);
fprintf('Number of droplets: %d\n', num_dots);
fprintf('Voltage range (applied): %.1f V to %.1f V\n', min(V), max(V));
fprintf('Vertical target span: %.1f mm\n', (max(y_target)-min(y_target))*1000);
disp('Simulation complete!');
