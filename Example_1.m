
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/9/2021
% Description:  This script generates Figures 1a and 1b.
%               Forward time reachable sets are predicted using MM.

clc; clear all; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Predict Reachable Sets using d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intervals Defining Initial Set
X0 = [-1/2 , 1/2; ....
      -1/2 , 1/2];
% Check to make sure X0 is a valid rectangle
if X0(1, 2) < X0(1, 1) || X0(2, 2) < X0(2, 1)
    print('Error 1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .002;   % Timestep for Simulation
T  = 1;      % Prediction Time-Horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0_Boundary = makeRectangle(X0);
Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
holder = Phi;
T_size = size(0:dt:T, 2);
xu = zeros(2, T_size + 1);  xu(:, 1) = X0(:, 1);
xo = zeros(2, T_size + 1);  xo(:, 1) = X0(:, 2);
 
% Compute Time = 1 Second Reachable Set of System
% Compute MM approximation of Reachable Set
for t = 1:T_size
    % reachable set computation
    holder2 = zeros(2, Phi_size);
    for i = 1:size(holder, 2)
            x = holder(:, i);
            x_next = x + dt*dxdt(x);
            holder2(:, i) = x_next;
    end
    holder = holder2;
    Phi = [Phi, holder2];
    
    % MM approximation
    xu(:, t + 1) = xu(:, t) + dt*d(xu(:, t), xo(:, t));
    xo(:, t + 1) = xo(:, t) + dt*d(xo(:, t), xu(:, t));
end
x_T = holder;
xu_T = xu(:, t + 1);
xo_T = xo(:, t + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca;
axis([-1, 5, -1, 3])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r' , ...
            'LineWidth', 1.25, ...
            'FaceAlpha', .7, ...
            'HandleVisibility', 'off');
% Plot MM Overapproximation of RF(1, X0)
SF = makeRectangle([xu_T, xo_T]);
patch(SF(1, :), SF(2, :), 'g', ...
                          'FaceAlpha', .1, ...
                          'LineWidth', 1.25, ...
                          'HandleVisibility', 'off');
% Plot Time = 1 Reachable Set RF(1, X0)
patch(x_T(1, :), x_T(2, :), 'g', ...
                            'FaceAlpha', .9, ...
                            'LineWidth', 1.25, ...
                            'HandleVisibility', 'off');
scatter(xu_T(1, 1), xu_T(2, 1), 'k', 'filled', 'HandleVisibility', 'off');
scatter(xo_T(1, 1), xo_T(2, 1), 'k', 'filled', 'HandleVisibility', 'off');


Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;
hold on; grid on;
ax = gca;
axis([-1, 5, -1, 5])
xticks([-1, 0, 1, 2, 3, 4, 5])
yticks([-1, 0, 1, 2, 3, 4, 5])
xlabel('$x_1$','Interpreter','latex')
ylabel('$\widehat{x}_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')

% Plot Cone Coresponding to X0
patch([X0(1, 1), X0(1, 1), X0(1, 2)], ...
      [X0(1, 2), X0(1, 1), X0(1, 2)], 'r', ...
            'LineWidth', 1, ...
            'FaceAlpha', .7, ...
            'HandleVisibility', 'off');
% Plot Cone Coresponding to Phi^e
patch([xu_T(1, 1), xu_T(1, 1), xo_T(1, 1)], ...
      [xo_T(1, 1), xu_T(1, 1), xo_T(1, 1)], 'g', ...
            'FaceAlpha', .2, ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');
% Plot trajectory of embedding system that yeilds approximation    
plot(xu(1, :), xo(1, :), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
scatter(X0(1, 1), X0(1, 2),60,  'b', 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');
scatter(xu_T(1, 1), xo_T(1, 1),60, 'b','filled', ...
                'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off');

% Plot Diagnol of Embedding Space
plot([-1, 5], [-1, 5], 'k', 'LineWidth', 1.25)

Leg = legend();
set(Leg,'visible','off')

grid on;
ax.Layer = 'top';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = d(x, xh)
    if x(2, 1) >= 0 && x(2, 1) >= - xh(2, 1)  
        out(1, 1) = x(2, 1)^2 + 2;   
    elseif xh(2, 1) <= 0 && x(2, 1) <= -xh(2, 1)
        out(1, 1) = xh(2, 1)^2 +  2;  
    elseif (x(2, 1) <= 0) && (xh(2, 1) >= 0)
        out(1, 1) = 2;
    end

    out(2, 1) = x(1); 
end

function out = dprime(x, xh)
    out = [xh(2)^2 + 2 + 10*(x(2) - xh(2)); ...
           x(1)]; 
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/10 :X0(1, 2), ...
                            X0(2, 1): d(2)/10 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function out = dxdt(x)
    out = [ x(2, 1)^2 + 2;...
            x(1, 1) ];
end