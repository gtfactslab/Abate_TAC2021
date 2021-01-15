
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 1/15/2021
% Description:  This script generates Figures 4a and 4b.
%               A robustly forward invariant and gloally attractive region 
%               is computed from an equilibrium in the embeddign space.  
%               This procedure is repeated using the backward-time dynamics 
%               to compute a second forward invariant region. 
%               We also employ a transformed system to compute forward
%               ivnariant regions aswell.


clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global W
% Bound on disturbance
W = [-3/4, 3/4; ...
     -3/4, 3/4];
 
global T
% Transformation matrix for parallelotope computation
T = [3,-1; ...
     1, 3];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import smallest attractive region which was computed externally 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('CS1_Data.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute forward invariant regions using MM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yeq = fsolve(@E, [-5; -5; 5; 5]);       % get forward time eq
yeq_back = fsolve(@EG, [-1; -1; 1; 1]); % get backward time eq

% create invariant regions
yeq = reshape(yeq, 2, 2)    % Invariant and attractive region
[a, b] = meshgrid(yeq(1, :), flip(yeq(2, :)));
X_EQ = makeRectangle(yeq);

yeq_back = reshape(yeq_back, 2, 2)  % Invariant region
[a, b] = meshgrid(yeq_back(1, :), flip(yeq_back(2, :)));
X_EQB = [a(:)'; b(1, 1), b(2, 1), b(2, 2), b(1, 2)];

if prod(yeq(:, 1) >= yeq(:, 2)) ||...
   prod(yeq_back(:, 1) >= yeq_back(:, 2))
    error('Equilibrium does not have correct ordering');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1.25, 1.25, -1.25, 1.25]);
xticks([-1 0 1])
yticks([-1 0 1])
Leg = legend();
set(Leg,'visible','off');
ax.Layer = 'top';

patch(X_EQ(1, :), X_EQ(2, :), 'w', ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', 1);
patch(X_EQ(1, :), X_EQ(2, :), 'g', ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', .5);

patch(outer_REACH(1, :), outer_REACH(2, :), 'w', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', 1);
patch(outer_REACH(1, :), outer_REACH(2, :), 'b',  ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', .5); % outer ring

patch(inner_REACH(1, :), inner_REACH(2, :), 'w', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', 1);
patch(inner_REACH(1, :), inner_REACH(2, :), 'g', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', .5); % inner ring

patch(X_EQB(1, :), X_EQB(2, :), 'w', ...
                                'HandleVisibility', 'off', ...
                                'FaceAlpha', 1);
patch(X_EQB(1, :), X_EQB(2, :), 'r', ...
                                'HandleVisibility', 'off', ...
                                'FaceAlpha', .5);

scatter(yeq(1, :), yeq(2, :), 'k', 'filled', ...
                              'HandleVisibility', 'off');
scatter(yeq_back(1, :), yeq_back(2, :), 'k', 'filled', ...
                                        'HandleVisibility', 'off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute forward invariant regions using MM with Transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimoptions('fsolve','OptimalityTolerance', 10^-10);
yeq = fsolve(@E2, 2*[-5; -5; 5; 5], options);       % get forward time eq
yeq_back = fsolve(@EG2, [-1; -1; 1; 1]); % get backward time eq

% create invariant regions
yeq = reshape(yeq, 2, 2)    % Invariant and attractive region
[a, b] = meshgrid(yeq(1, :), flip(yeq(2, :)));
X_EQ = T*makeRectangle(yeq);

yeq = T*yeq;

yeq_back = reshape(yeq_back, 2, 2)  % Invariant region
[a, b] = meshgrid(yeq_back(1, :), flip(yeq_back(2, :)));
X_EQB = T*makeRectangle(yeq_back);

yeq_back = T*yeq_back;

if prod(yeq(:, 1) >= yeq(:, 2)) ||...
   prod(yeq_back(:, 1) >= yeq_back(:, 2))
    error('Equilibrium does not have correct ordering');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;
hold on; grid on;
ax = gca;
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1.5, 1.5, -1.5, 1.5]);
xticks([-1 0 1])
yticks([-1 0 1])
Leg = legend();
set(Leg,'visible','off');
ax.Layer = 'top';

patch(X_EQ(1, :), X_EQ(2, :), 'w', ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', 1);
patch(X_EQ(1, :), X_EQ(2, :), 'g', ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', .5);

patch(outer_REACH(1, :), outer_REACH(2, :), 'w', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', 1);
patch(outer_REACH(1, :), outer_REACH(2, :), 'b',  ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', .5); % outer ring

patch(inner_REACH(1, :), inner_REACH(2, :), 'w', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', 1);
patch(inner_REACH(1, :), inner_REACH(2, :), 'g', ...
                                            'HandleVisibility', 'off', ...
                                            'FaceAlpha', .5); % inner ring
patch(X_EQB(1, :), X_EQB(2, :), 'w', ...
                                'HandleVisibility', 'off', ...
                                'FaceAlpha', 1);
patch(X_EQB(1, :), X_EQB(2, :), 'r', ...
                                'HandleVisibility', 'off', ...
                                'FaceAlpha', .5);

scatter(yeq(1, :), yeq(2, :), 'k', 'filled', ...
                              'HandleVisibility', 'off');
scatter(yeq_back(1, :), yeq_back(2, :), 'k', 'filled', ...
                                        'HandleVisibility', 'off');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end


%%%%%%%%%%%%%%%%%%%%
% without transformation
%%%%%%%%%%%%%%%%%%%%
% dynamics 
function out = F(x, w)
    out = [-x(2) + x(1)*(4 - 4*x(1)^2 - x(2)^2) + w(1); ...
            x(1) + x(2)*(4 - x(1)^2 - 4*x(2)^2) + w(2)];
end 

% forward-time decomposition function
function out = d(y, w, yhat, what)
    
    if prod([y; w] <= [yhat; what])
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) [1, 0] * F([y(1); x], [0; 0]);
        fun_y2 = @(x) [0, 1] * F([x; y(2)], [0; 0]);
        fun_w1 = @(v) [1, 0] * F([0; 0], [v; 0]);
        fun_w2 = @(v) [0, 1] * F([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w(1), what(1), options);
        out_w2 = fminbnd(fun_w2, w(2), what(2), options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
    elseif prod([y; w] >= [yhat; what])
        fun_y1 = @(x) -[1, 0]* F([y(1); x], [0; 0]);
        fun_y2 = @(x) -[0, 1]* F([x; y(2)], [0; 0]);
        fun_w1 = @(v) -[1, 0]* F([0; 0], [v; 0]);
        fun_w2 = @(v) -[0, 1]* F([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what(1), w(1), options);
        out_w2 = fminbnd(fun_w2, what(2), w(2), options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
    else
        out = [inf; inf];
    end
end

% forward-time embedding function
function out = E(xm)
    global W
    x = xm(1:2); 
    xhat = xm(3:4);
    
    out = [d(x, W(:, 1), xhat, W(:, 2)); ...
           d(xhat, W(:, 2), x, W(:, 1))];
end

% backward-time decomposition function
function out = dG(y, w, yhat, what)
    
    if prod([y; w] <= [yhat; what])
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) -[1, 0] * F([y(1); x], [0; 0]);
        fun_y2 = @(x) -[0, 1] * F([x; y(2)], [0; 0]);
        fun_w1 = @(v) -[1, 0] * F([0; 0], [v; 0]);
        fun_w2 = @(v) -[0, 1] * F([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w(1), what(1), options);
        out_w2 = fminbnd(fun_w2, w(2), what(2), options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
    elseif prod([y; w] >= [yhat; what])
        fun_y1 = @(x) [1, 0]* F([y(1); x], [0; 0]);
        fun_y2 = @(x) [0, 1]* F([x; y(2)], [0; 0]);
        fun_w1 = @(v) [1, 0]* F([0; 0], [v; 0]);
        fun_w2 = @(v) [0, 1]* F([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what(1), w(1), options);
        out_w2 = fminbnd(fun_w2, what(2), w(2), options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
        
    else
        out = [inf; inf];
    end
end

% backward-time embedding function
function out = EG(xm)
    global W
    x = xm(1:2); xhat = xm(3:4);
    out = [dG(x, W(:, 1), xhat, W(:, 2)); ...
           dG(xhat, W(:, 2), x, W(:, 1))];
end



%%%%%%%%%%%%%%%%%%%%
% with transformation
%%%%%%%%%%%%%%%%%%%%

function out = F2(x, w)
    global T
    y = T*x;
    
    out = [-y(2) + y(1)*(4 - 4*y(1)^2 - y(2)^2) + w(1); ...
            y(1) + y(2)*(4 - y(1)^2 - 4*y(2)^2) + w(2)];
        
    out = T\out; 
end 

% forward-time decomposition function with transformation
function out = d2(y, w, yhat, what)
    
    if prod([y; w] <= [yhat; what])
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) [1, 0] * F2([y(1); x], [0; 0]);
        fun_y2 = @(x) [0, 1] * F2([x; y(2)], [0; 0]);
        fun_w1 = @(v) [1, 0] * F2([0; 0], [v; 0]);
        fun_w2 = @(v) [0, 1] * F2([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w(1), what(1), options);
        out_w2 = fminbnd(fun_w2, w(2), what(2), options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
    elseif prod([y; w] >= [yhat; what])
        fun_y1 = @(x) -[1, 0]* F2([y(1); x], [0; 0]);
        fun_y2 = @(x) -[0, 1]* F2([x; y(2)], [0; 0]);
        fun_w1 = @(v) -[1, 0]* F2([0; 0], [v; 0]);
        fun_w2 = @(v) -[0, 1]* F2([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what(1), w(1), options);
        out_w2 = fminbnd(fun_w2, what(2), w(2), options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
    else
        out = [inf; inf];
    end
end

% forward-time embedding function with transformation
function out = E2(xm)
    global W
    x = xm(1:2); 
    xhat = xm(3:4);
    
    out = [d2(x, W(:, 1), xhat, W(:, 2)); ...
           d2(xhat, W(:, 2), x, W(:, 1))];
end

% backward-time decomposition function with transformation
function out = dG2(y, w, yhat, what)
    
    if prod([y; w] <= [yhat; what])
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) -[1, 0] * F2([y(1); x], [0; 0]);
        fun_y2 = @(x) -[0, 1] * F2([x; y(2)], [0; 0]);
        fun_w1 = @(v) -[1, 0] * F2([0; 0], [v; 0]);
        fun_w2 = @(v) -[0, 1] * F2([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w(1), what(1), options);
        out_w2 = fminbnd(fun_w2, w(2), what(2), options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
    elseif prod([y; w] >= [yhat; what])
        fun_y1 = @(x) [1, 0]* F2([y(1); x], [0; 0]);
        fun_y2 = @(x) [0, 1]* F2([x; y(2)], [0; 0]);
        fun_w1 = @(v) [1, 0]* F2([0; 0], [v; 0]);
        fun_w2 = @(v) [0, 1]* F2([0; 0], [0; v]);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what(1), w(1), options);
        out_w2 = fminbnd(fun_w2, what(2), w(2), options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
        
    else
        out = [inf; inf];
    end
end

% backward-time embedding function with transformation
function out = EG2(xm)
    global W
    x = xm(1:2); xhat = xm(3:4);
    out = [dG2(x, W(:, 1), xhat, W(:, 2)); ...
           dG2(xhat, W(:, 2), x, W(:, 1))];
end






