
% Title: Computing Robustly Forward Invariant Sets for Mixed-Monotone
%        Systems
% Submitted to: Transactions on Automatic Control (TAC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/9/2021
% Description:  This script generates Figure 3.
%               A forward invariant set is computed for a nonlinear
%               quadrotor model using the theory of mixed-monotonicity.  A
%               linear transformation of the dynamics is formed, which is
%               mixed monotone, and an equilibrium in the embedding space
%               identifies a forward invariant parallelotope for the
%               initial dynamics.

clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%
% Define Dynamics
%%%%%%%%%%%%%%%%%%%%%%%
g = -9.8; % Gravitaional Constant

% Linearisation Matixes at the origin
A = [0, 0, 1, 0, 0; ...
     0, 0, 0, 1, 0; ...
     0, 0, 0, 0, g; ...
     0, 0, 0, 0, 0; ...
     0, 0, 0, 0, 0];

B = [0, 0; ...
     0, 0; ...
     0, 0; ...
     1, 0; ...
     0, 1];

K = [0, -40, 0, -20,   0; ...
     3,   0, 6,   0, -25];

%%%%%%%%%%%%%%%%%%%%%%%
% Get Transformation From Linearisation Of xdot = F(x) At The Origin x = 0
%%%%%%%%%%%%%%%%%%%%%%%

J_cl = A+B*K;

[T, ~] = eig(J_cl);

%%%%%%%%%%%%%%%%%%%%%%%
% Hypothesize Invariant Set Geometry
%%%%%%%%%%%%%%%%%%%%%%%

xe = [0; 0; 0; 1; 3.6075];

% confirm set invariant
if SE(zeros(10,1), E(-xe, xe))
    disp('Set Is Invariant');
else
    disp('Set Is Not Invariant');
end

%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Invariant Set
%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on; grid on;
drawnow;
ax = gca;
xlabel('$x_2$','Interpreter','latex')
ylabel('$x_4$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
Leg = legend();
set(Leg,'visible','off');
ax.Layer = 'top';

EE = [-xe; xe];

% get corners of parellelogram in original space
corn = T([2, 4], [4, 5])*[EE(4), EE(9), EE(9), EE(4), EE(4); ...
                          EE(5), EE(5), EE(10), EE(10), EE(5)];

% Plot Invariant Parallelogram
patch(corn(1, :), corn(2, :), [0, 1, 0], ...
                              'HandleVisibility', 'off', ...
                              'FaceAlpha', .4);

% Plot collums of T
plot(-2*[0, T(2, 4)], -2*[0, T(4, 4)], 'k', 'LineWidth', 2)
plot(1*[0, T(2, 5)], 1*[0, T(4, 5)], 'k', 'LineWidth', 2)

function out = F(x)
z1 = x(1);
z2 = x(2);
z3 = x(3);
z4 = x(4);
z5 = x(5);

out(1, 1) = 0.1306*z2 - 25.3218*z1 - 0.5178*z3 + 0.3212*sin(0.9163*z1 - 0.0415*z2 + 0.1645*z3)*(17.7179*z4 + 2.0604*z5 + 9.8000);
out(2, 1) = 10.6888*z1 - 1.1900*z2 + 1.9189*z3 - 1.1904*sin(0.9163*z1 - 0.0415*z2 + 0.1645*z3)*(17.7179*z4 + 2.0604*z5 + 9.8000);
out(3, 1) = 18.7611*z1 - 0.8497*z2 + 1.5118*z3 - 2.0894*sin(0.9163*z1 - 0.0415*z2 + 0.1645*z3)*(17.7179*z4 + 2.0604*z5 + 9.8000);
out(4, 1) = 2.5820*z4 + 2.3639*z5 - 1.1473*cos(0.9163*z1 - 0.0415*z2 + 0.1645*z3)*(17.7179*z4 + 2.0604*z5 + 9.8000) + 11.2437;
out(5, 1) = 0.1592*cos(0.9163*z1 - 0.0415*z2 + 0.1645*z3)*(17.7179*z4 + 2.0604*z5 + 9.8000) - 2.5820*z5 - 2.8202*z4 - 1.5599; 
end

% Southeast Order
function out = SE(a, b)
    n = size(a, 1)/2;

    if prod(a(1:n) <= b(1:n)) && prod(a(n+1:2*n) >= b(n+1:2*n))
        out = 1;
    else
        out = 0;
    end
end

% Embedding Function
function out = E(x, xh)
    out = [d(x, xh); d(xh, x)];
end

% Decomposition Function
function out = d(x, xh)
    out = [d1(x, xh);...
           d2(x, xh);...
           d3(x, xh);...
           d4(x, xh);...
           d5(x, xh)];
end

function out = d1(x, xh)
z1 = x(1);
z2 = x(2);
z3 = x(3);
z4 = x(4);
z5 = x(5);
zh1 = xh(1);
zh2 = xh(2);
zh3 = xh(3);
zh4 = xh(4);
zh5 = xh(5);

A1 = [z1, z1];
A2 = [z2, zh2];
A3 = [z3, zh3];
A4 = [z4, zh4];
A5 = [z5, zh5];

q1 = [];
q2 = [];
q3 = [];
q4 = [];


for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
%
q1 = [q1, 0.3212*17.7179*A4(i4)*(sin(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) + 0.0415*A2(i2) + 0.1645*A3(i3))];
q2 = [q2, - 0.3212*17.7179*A4(i4)*(0.0415*A2(i2) + 0.1645*A3(i3))];
%
q3 = [q3, 0.3212*2.0604*A5(i5)*(sin(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) +  0.0415*A2(i2) + 0.1645*A3(i3))];
q4 = [q4, - 0.3212*2.0604*A5(i5)*( 0.0415*A2(i2) + 0.1645*A3(i3))];

                end
            end
        end
    end
end

if prod(x <= xh)
    out(1, 1) = min(q1) + min(q2) + ...
                min(q3) + min(q4);
elseif prod(xh <= x)
    out(1, 1) = max(q1) + max(q2) + ...
                max(q3) + max(q4);
else 
    error
end

out(1, 1) = out(1, 1) + 0.3212*9.8000*sin(0.9163*z1 - 0.0415*z2 + 0.1645*z3) ...
                      + 0.0415*0.3212*9.8000*(z2 - zh2) ...
                      + 0.1645*0.3212*9.8000*(z3 - zh3);



out(1, 1) = out(1, 1) + 0.1306*z2 - 25.3218*z1 - 0.5178*zh3;

end

function out = d2(x, xh)
z1 = x(1);
z2 = x(2);
z3 = x(3);
z4 = x(4);
z5 = x(5);
zh1 = xh(1);
zh2 = xh(2);
zh3 = xh(3);
zh4 = xh(4);
zh5 = xh(5);

A1 = [z1, zh1];
A2 = [z2, z2];
A3 = [z3, zh3];
A4 = [z4, zh4];
A5 = [z5, zh5];

q1 = [];
q2 = [];
q3 = [];
q4 = [];


for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
%
q1 = [q1, - 1.1904*17.7179*A4(i4)*(sin(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) - 0.9163*A1(i1) - 0.1645*A3(i3))];
q2 = [q2, - 1.1904*17.7179*A4(i4)*(0.9163*A1(i1) + 0.1645*A3(i3))];

%
q3 = [q3, - 1.1904*2.0604*A5(i5)*(sin(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) - 0.9163*A1(i1) - 0.1645*A3(i3))];
q4 = [q4, - 1.1904*2.0604*A5(i5)*(0.9163*A1(i1) + 0.1645*A3(i3))];

                end
            end
        end
    end
end

if prod(x <= xh)
    out(1, 1) = min(q1) + min(q2) + ...
                min(q3) + min(q4);
elseif prod(xh <= x)
    out(1, 1) = max(q1) + max(q2) + ...
                max(q3) + max(q4);
end

out(1, 1) = out(1, 1) - 1.1904*9.8000*sin(0.9163*z1 - 0.0415*z2 + 0.1645*z3) ...
                      + 0.9163*1.1904*9.8000*(z1 - zh1) ...
                      + 0.1645*1.1904*9.8000*(z3 - zh3);

out(1, 1) = out(1, 1) + 10.6888*z1 - 1.1900*z2 + 1.9189*z3;


end

function out = d3(x, xh)
z1 = x(1);
z2 = x(2);
z3 = x(3);
z4 = x(4);
z5 = x(5);
zh1 = xh(1);
zh2 = xh(2);
zh3 = xh(3);
zh4 = xh(4);
zh5 = xh(5);

A1 = [z1, zh1];
A2 = [z2, zh2];
A3 = [z3, z3];
A4 = [z4, zh4];
A5 = [z5, zh5];

q1 = [];
q2 = [];
q3 = [];
q4 = [];


for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
%
q1 = [q1, - 2.0894*17.7179*A4(i4)*(sin(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) + 0.9163*A1(i1) - 0.0415*A2(i2))];
q2 = [q2,  2.0894*17.7179*A4(i4)*(0.9163*A1(i1) - 0.0415*A2(i2))];
%
q3 = [q3, - 2.0894*2.0604*A5(i5)*(sin(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) + 0.9163*A1(i1) - 0.0415*A2(i2))];
q4 = [q4,  2.0894*2.0604*A5(i5)*(0.9163*A1(i1) - 0.0415*A2(i2))];

                end
            end
        end
    end
end

if prod(x <= xh)
    out(1, 1) = min(q1) + min(q2) + ...
                min(q3) + min(q4);
elseif prod(xh <= x)
    out(1, 1) = max(q1) + max(q2) + ...
                max(q3) + max(q4);
end

out(1, 1) = out(1, 1) - 2.0894*9.8000*sin(0.9163*z1 - 0.0415*z2 + 0.1645*z3) ...
                      + 0.9163*2.0894*9.8000*(z1 - zh1) ...
                      + 0.0415*2.0894*9.8000*(z2 - zh2);

out(1, 1) = out(1, 1) + 18.7611*z1 - 0.8497*zh2 + 1.5118*z3;

end

function out = d4(x, xh)
z1 = x(1);
z2 = x(2);
z3 = x(3);
z4 = x(4);
z5 = x(5);
zh1 = xh(1);
zh2 = xh(2);
zh3 = xh(3);
zh4 = xh(4);
zh5 = xh(5);

A1 = [z1, zh1];
A2 = [z2, zh2];
A3 = [z3, zh3];
A4 = [z4, z4];
A5 = [z5, zh5];

q1 = [];
q2 = [];
q3 = [];
q4 = [];


for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
%
q1 = [q1, - 1.1473*17.7179*A4(i4)*(cos(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) - 0.9163*A1(i1) + 0.0415*A2(i2) - 0.1645*A3(i3))];
q2 = [q2, - 1.1473*17.7179*A4(i4)*(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3))];
%
q3 = [q3, - 1.1473*2.0604*A5(i5)*(cos(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) - 0.9163*A1(i1) + 0.0415*A2(i2) - 0.1645*A3(i3))];
q4 = [q4, - 1.1473*2.0604*A5(i5)*(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3))];

                end
            end
        end
    end
end

if prod(x <= xh)
    out(1, 1) = min(q1) + min(q2) + ...
                min(q3) + min(q4);
elseif prod(xh <= x)
    out(1, 1) = max(q1) + max(q2) + ...
                max(q3) + max(q4);
           
else
    error();
end

out(1, 1) = out(1, 1) - 1.1473*9.8000*cos(0.9163*z1 - 0.0415*z2 + 0.1645*z3) ...
                      + 0.9163*1.1473*9.8000*(z1 - zh1) ...
                      + 0.0415*1.1473*9.8000*(z2 - zh2) ...
                      + 0.1645*1.1473*9.8000*(z3 - zh3);

out(1, 1) = out(1, 1) + 2.5820*z4 + 2.3639*z5 + 11.2437;

end

function out = d5(x, xh)
z1 = x(1);
z2 = x(2);
z3 = x(3);
z4 = x(4);
z5 = x(5);
zh1 = xh(1);
zh2 = xh(2);
zh3 = xh(3);
zh4 = xh(4);
zh5 = xh(5);

A1 = [z1, zh1];
A2 = [z2, zh2];
A3 = [z3, zh3];
A4 = [z4, zh4];
A5 = [z5, z5];

q1 = [];
q2 = [];
q3 = [];
q4 = [];


for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
    %
    q1 = [q1, 0.1592*17.7179*A4(i4)*(cos(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) + 0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3))];
    q2 = [q2, - 0.1592*17.7179*A4(i4)*(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3))];
    %
    q3 = [q3, 0.1592*2.0604*A5(i5)*(cos(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3)) + 0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3))];
    q4 = [q4, - 0.1592*2.0604*A5(i5)*(0.9163*A1(i1) - 0.0415*A2(i2) + 0.1645*A3(i3))];
    
                    end
                end
            end
        end
    end
    
    if prod(x <= xh)
        out(1, 1) = min(q1) + min(q2) + min(q3) + min(q4);
    elseif prod(xh <= x)
        out(1, 1) = max(q1) + max(q2) + max(q3) + max(q4);
    end
    
    out(1, 1) = out(1, 1) + 0.1592*9.8000*cos(0.9163*z1 - 0.0415*z2 + 0.1645*z3) ...
                          + 0.9163*0.1592*9.8000*(z1 - zh1) ...
                          + 0.0415*0.1592*9.8000*(z2 - zh2) ...
                          + 0.1645*0.1592*9.8000*(z3 - zh3);
    
    out(1, 1) = out(1, 1) - 2.5820*z5 - 2.8202*zh4 - 1.5599;

end





