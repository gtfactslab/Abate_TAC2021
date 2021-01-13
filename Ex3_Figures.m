
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 1/13/2021
% Description:  This script generates Figure 3.
%               Forward and backward-time reachable sets are approximated
%               via the MM property.

clc; clear all;

% Initial Set
X0 = [ -.25 ,.25; ....
       -.25 , 0.25 ];
   
% Disturbance Bound
W = [0, .25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .005;   % Timestep for simulation
T  = 1;      % Simulation time-horizon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0_Boundary = makeRectangle(X0);
Phi0 = X0_Boundary;
Phi_size = size(Phi0, 2);

Phi = Phi0;
nPhi = Phi0;

holder = Phi;
nholder = Phi;

T_size = size(0:dt:T, 2);
xu = zeros(2, T_size + 1);  xu(:, 1) = X0(:, 1);
xo = zeros(2, T_size + 1);  xo(:, 1) = X0(:, 2);

nxu = zeros(2, T_size + 1);  nxu(:, 1) = X0(:, 1);
nxo = zeros(2, T_size + 1);  nxo(:, 1) = X0(:, 2);
 
for t = 1:T_size
    % print current time
    current_time = dt*(t-1)
    
    % forward-time reachable set computation
    holder2 = [];
    for i = 1:size(holder, 2)
            x = holder(:, i);
            for w = W(1):(W(2) - W(1))/4:W(2)
                x_next = x + dt*F(x, w);
                holder2 = [holder2, x_next];
            end
    end
    k = boundary(holder2(1, :)',holder2(2, :)',0.02);
    holder = holder2(:, k);
    Phi = [Phi, holder];
    
    % forward-time MM approximation
    xu(:, t + 1) = xu(:, t) + dt*d(xu(:, t), W(1), xo(:, t), W(2));
    xo(:, t + 1) = xo(:, t) + dt*d(xo(:, t), W(2), xu(:, t), W(1));
    
    % backward-time reachable set computation
    nholder2 = [];
    for i = 1:size(nholder, 2)
            nx = nholder(:, i);
            for w = W(1):.25:W(2)
                nx_next = nx - dt*F(nx, w);
                nholder2 = [nholder2, nx_next];
            end
    end
    k = boundary(nholder2(1, :)', nholder2(2, :)',0.02);
    nholder = nholder2(:, k);
    nPhi = [nPhi, nholder];
    
    % backward-time MM approximation
    nxu(:, t + 1) = nxu(:, t) + dt*nd(nxu(:, t), W(1), nxo(:, t), W(2));
    nxo(:, t + 1) = nxo(:, t) + dt*nd(nxo(:, t), W(2),  nxu(:, t), W(1));
    
end
x_T = holder;
xu_T = xu(:, t + 1);
xo_T = xo(:, t + 1);

nx_T = nholder;
nxu_T = nxu(:, t + 1);
nxo_T = nxo(:, t + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1, 3.5, -2.5, 2.5]);
xticks([-1 0 1 2 3])
yticks([-2 -1 0 1 2])

% Plot Initial Set X0
patch(X0_Boundary(1, :), X0_Boundary(2, :), 'r', 'LineWidth', 1.25);

% Plot forward time reachable set and MM approximation
SF = makeRectangle([xu_T, xo_T]);
patch(SF(1, :), SF(2, :), 'g', 'FaceAlpha', .1, 'LineWidth', 1.25);
patch(x_T(1, :), x_T(2, :), 'g', 'FaceAlpha', .7, 'LineWidth', 1.25);

scatter([xu_T(1, 1), xo_T(1, 1)], [xu_T(2, 1), xo_T(2, 1)], 'k', 'filled');

% Plot backward time reachable set and MM approximation
NSF = makeRectangle([nxu_T, nxo_T]);
patch(NSF(1, :), NSF(2, :), 'b', 'FaceAlpha', .05, 'LineWidth', 1.25);
patch(nx_T(1, :), nx_T(2, :), 'b', 'FaceAlpha', .6, 'LineWidth', 1.25);
scatter(nxu_T(1, 1), nxu_T(2, 1), 'k', 'filled');
scatter(nxo_T(1, 1), nxo_T(2, 1), 'k', 'filled');

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% system dynamics
function out = F(x, w)
    out = [  x(1, 1)*x(2, 1) + x(2, 1) + w;...
            x(1, 1) + 1];
end

% tight decomposition function for forward-time dynamics
function out = d(x, w, xhat, what)
    if x(1, 1) >= -1
        out(1, 1) = x(1, 1)*x(2, 1) + x(2, 1) + w;
    else
        out(1, 1) = x(1, 1)*xhat(2, 1) + xhat(2, 1) + w;
    end
    out(2, 1) = x(1, 1) + 1;
end

% tight decomposition function for backard-time dynamics
function out = nd(x, w, xhat, what)
    if x(1, 1) <= -1
        out(1, 1) = -x(1, 1)*x(2, 1) - x(2, 1) - what;
    else
        out(1, 1) = -x(1, 1)*xhat(2, 1) - xhat(2, 1) - what;
    end
    out(2, 1) = - xhat(1, 1) - 1;
end

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end