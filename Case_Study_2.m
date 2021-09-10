
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/9/2021
% Description:  This script generates Figure 4.
%               A robustly forward invariant and attractive set for a 
%               4-dimensional polynomial system is identified from a limit 
%               cycle in the embedding system.  3 reachable sets are
%               computed and plotted to demonstrate the attractiveness of
%               the invariant set.

clc; clear all;

% Disturbance bound
global W
W = [-.1, .1]; % Disturbance Bound

% Initial Rectangle
X0 = [1, 1.5;
      1, 1.5;
      1, 1;
      0, 0]
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on
axis([-.6, 1.6, -.5, 1.6])
xlabel('$x_1$','Interpreter','latex')
xticks([-.5, 0, .5, 1, 1.5])
ylabel('$x_2$','Interpreter','latex')
yticks([-.5, 0, .5, 1, 1.5])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Invariant Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Case_Study_2_Data.mat') % Get set from previous simulation
patch(outer(1, :), outer(2, :), 'b', 'FaceAlpha', .3, 'LineWidth', 1.3, 'HandleVisibility', 'off');  % Plot
patch(inner(1, :), inner(2, :), 'w', 'LineWidth', 1.3, 'HandleVisibility', 'off');                   % Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get initial points for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = .005;
[a1, a2] = meshgrid(X0(1, 1):dx:X0(1, 2), ...
                    X0(2, 1):dx:X0(2, 2));
a1 = a1(:);
a2 = a2(:);

k = boundary(a1, a2);
a1 = a1(k);
a2 = a2(k);
Initial_Points = zeros(4, size(a1, 1));
for i = 1:1:size(a1, 1)
    Initial_Points(:, i) = [a1(i); a2(i); X0(3, 1); X0(4, 1)];
end
clear a1 a2 k

% plot initial set
patch(Initial_Points(1, :), Initial_Points(2, :), 'w')
patch(Initial_Points(1, :), Initial_Points(2, :), 'r', 'LineWidth', 1.3, 'FaceAlpha', .7)
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct Simulation (To Get True Reachable Set)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = .001;
T_fin = 14;
dw = .05;

holder = Initial_Points;
holder2 = [];

Embedding_State = X0(:);
corner_holder = [];

for t = 0:dt:T_fin
    t
    
    if t ~= 0
        holder = Next_Set;
    end
    
    if t == .05
        dw = .2;
    end
    
    % propegate reachable set forward in time
    for i = 1:size(holder, 2)
        x_now = holder(:, i);
        for w1 = W(1, 1) : dw : W(1, 2)
            x_next = x_now + dt*F(x_now, w1);
            holder2 = [holder2, x_next];
        end
    end
    %k = boundary(holder2(1, :)',holder2(2, :)', fid);
    [k,av] = convhull(holder2([1,2], :)');
    Next_Set = holder2(:, k);
    holder2 = [];

    % propegate embedding state forward in time
    Embedding_State = Embedding_State + dt*E(Embedding_State);

    if t == 10 || t == 12 || t == 14 
        % plot specific reachable sets
        rect = makeRectangle([Embedding_State(1:2), Embedding_State(5:6)]);
        patch(rect(1, :), rect(2, :), 'w', 'HandleVisibility', 'off');
        patch(rect(1, :), rect(2, :), 'g', 'FaceAlpha', .2, 'HandleVisibility', 'off');
        patch(Next_Set(1, 1:2:end), Next_Set(2, 1:2:end), 'g', 'HandleVisibility', 'off');
        drawnow
    end
end

drawnow
Leg = legend();
set(Leg,'visible','off')



% true dynamics
function out = F(x, w)
    out = [-2*x(1) + x(2)*(1+x(1)) + x(3) + w(1); ...
           -x(2) + (1 - x(2))*x(1) + 0.1; ...
           -x(4); ...
           x(3)];
end

% embedding function
function out = E(a)
    global W
    x = a(1:4);
    xhat = a(5:8);
    
    out = [decomp(x,W(:, 1), xhat, W(:, 2)); ...
           decomp(xhat, W(:, 2), x, W(:, 1))];
end

% decomposition function
function out = decomp(x, w, xhat, what)
    if x(1) >= -1
        out(1, 1) = -2*x(1) + x(2)*(1+x(1)) + x(3) + w(1);
    else
        out(1, 1) = -2*x(1) + xhat(2)*(1+x(1)) + x(3) + w(1);
    end
    
    if (1 - x(2)) >= 0 
        out(2, 1) = -x(2) + (1 - x(2))*x(1) + 0.1;
    else
        out(2, 1) = -x(2) + (1 - x(2))*xhat(1) + 0.1;
    end
    out(3:4, 1) = [-xhat(4); x(3)];
end

% plotting tool
function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

    