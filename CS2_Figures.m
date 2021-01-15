
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 1/15/2021
% Description:  This script generates Figures 5.
%               A decomposition function for x^+ = sin(x) is computed and
%               plotted on [-10, 10]x[-10, 10].  A forward invariant region
%               for the system in Case Study 2 is also computed.


clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global W
W = [-2, 2]; % Didurbance Bound


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Hyperrectangular RFI Set Using MM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute an equilibrium in embedding space.
% This point defines a robustly forward invairant set for the original
% dynamics.
thing = fsolve(@thng, -ones(4, 1));
XE = reshape(thing, 2, 2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot decomposition function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = .4;
X_range = [-10, 10];

[X, Xh] = meshgrid(X_range(1):dx:X_range(2), X_range(1):dx:X_range(2));
num = size(X, 1);
               
points = [X(:)'; Xh(:)'];

real = [];
% get diagnol
for x = X_range(1):dx:X_range(2)
    out = sin(x);
    real = [real, [x; x; out]];
end

holder = [];
for i = 1:size(points, 2)
    xnow = points(:, i);
    out2 = delta(xnow(1), xnow(2));
    
    holder = [holder, [xnow; out; out2]];    
end
X = reshape(holder(1, :), [num, num]);
Y = reshape(holder(2, :), [num, num]);
Z1 = reshape(holder(3, :), [num, num]);
Z2 = reshape(holder(4, :), [num, num]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
ax = gca;
xlabel('$y$','Interpreter','latex')
ylabel('$\widehat{y}$','Interpreter','latex')
zlabel('$\delta(y,\, y)$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis(10*[-1, 1, -1, 1]);
Leg = legend();
set(Leg,'visible','off');

view([-1, .5, 1])
%surf(X, Y, Z2, 'FaceColor', 'g')
s = surf(X, Y, Z2, 'EdgeColor', 'none', 'FaceColor', 'interp');
plot3(real(1, :),real(2, :),real(3, :), 'r', 'LineWidth', 2);

%%Extract X,Y and Z data from surface plot
x=s.XData;
y=s.YData;
z=s.ZData;
%%Create vectors out of surface's XData and YData
x=x(1,:);
y=y(:,1);
%%Divide the lengths by the number of lines needed
xnumlines = 10; % 10 lines
ynumlines = 10; % 10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);
%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k');
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dynamics
function out = F(x, w)
    out = [sin(x(2)); 
           1/2*x(1) + w];
end

% function for computing equilibrium in discrete-time
function out = thng(a)
    out = E(a(1:2), a(3:4)) - a;
end

% embedding function
function out = E(x, xhat)
    global W
    out = [ d(x, W(:, 1), xhat, W(:, 2)); 
            d(xhat, W(:, 2), x, W(:, 1)) ];
end

% decomposition function
function out = d(x, w, xhat, what)
    out = [ delta(x(2), xhat(2)); 
           1/2*x(1) + w];
end

% decomposition function for sin(x)
function out = delta(y, yhat)
    if abs(y - yhat) >= 2*pi
        out = sign(y - yhat);
    elseif cos(y) <= 0 && cos(yhat) <= 0 && abs(y - yhat) <= pi
        out = sin(yhat);
    elseif cos(y) <= 0 && cos(yhat) >= 0
        out = sign(y - yhat);
    elseif cos(y) >= 0 && cos(yhat) >= 0 && abs(y - yhat) <= pi
        out = sin(y);
    elseif cos(y)*cos(yhat) >= 0 && abs(y - yhat) >= pi
        out = sign(y - yhat);
    elseif yhat >= y && cos(y) >= 0 && cos(yhat) <= 0
        out = min([sin(y), sin(yhat)]);
    elseif y >= yhat && cos(y) >= 0 && cos(yhat) <= 0
        out = max([sin(y), sin(yhat)]);
    end
end
