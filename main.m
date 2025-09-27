clear; close all; 

Nx = 10;   % No.of cells in x-direction 
Ny = 10;   % No.of cells in y-direction 
Lx = 1;    % Length of the mesh in x-direction 
Ly = 0.5;  % Length of the mesh in y-direction
%% Plotting Mesh
[x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly);

figure(1);
xline(x_lineee);
yline(y_lineee);
hold on;
for c = 1:length(nodesx)
    plot(nodesx(c), nodesy,'r.','MarkerSize',10);
end
hold on;
for r = 1:length(nodesx)
    plot(nodesx, nodesy(r),'r.','MarkerSize',10);
end
hold off;
%%
S = -1.5; q = -5000;
epsilon = 1e-4; % Error Tolerance for convergence 

T = zeros(Ny+2, Nx+2);
num_iter = 100000;
for n = 1:num_iter
    fprintf('Iteration no. %d \n', n)
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly);
    [T, cnvrg_val] = GaussScidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly, S, q);
    residuals(n) = cnvrg_val;
    disp(cnvrg_val) 

    if cnvrg_val < epsilon
        disp('Convergence has been achieved')
        break
    else
        disp('Not yet, try harder')
    end
end
%% Checking mesh independence
% Nx = [10, 20, 40, 80];
% Ny = [10, 20, 40, 80];
% 
% S = -1.5; q = -5000;
% epsilon = 1e-4; % Error Tolerance for convergence 
% num_iter = 100000;
% for l = 1:length(Nx)
%     fprintf('Solving for %dx%d \n', Nx(l), Ny(l))
%     [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx(l), Ny(l), Lx, Ly);
%     T = zeros(Ny(l)+2, Nx(l)+2);
% 
%     for n = 1:num_iter
%         [T, cnvrg_val] = GaussScidelSolver(T, nodesx, nodesy, dx, dy, Nx(l), Ny(l), Lx, Ly, S, q);
%         if cnvrg_val < epsilon
%             break
%         else
%             continue
%         end
%     end
% 
%     [X, Y] = meshgrid(nodesx, nodesy);
% 
%     avgT_LeftSec(l)   = mean(T(X <= Lx/3));
%     avgT_MiddleSec(l) = mean(T(X > Lx/3 & X <= 2*Lx/3));
%     avgT_RightSec(l)  = mean(T(X > 2*Lx/3));
%     disp('Done')
% end
% 
% fprintf('Std Deviation of Avg Temp in left section is %.2f \n', std(avgT_LeftSec))
% fprintf('Std Deviation of Avg Temp in middle section is %.2f \n', std(avgT_MiddleSec))
% fprintf('Std Deviation of Avg Temp in right section is %.2f \n', std(avgT_RightSec))
%%
dTdx = zeros(Ny+2, Nx+2);
dTdy = zeros(Ny+2, Nx+2);

[X_nodes, Y_nodes] = meshgrid(nodesx, nodesy);
for i = 1:Ny+2 
    dTdx(i, 1) = (T(i, 2) - T(i, 1)) / (nodesx(2) - nodesx(1));

    for j = 2:Ny+1
        dTdx(i, j) = (T(i, j+1) - T(i, j-1)) / (nodesx(j+1) - nodesx(j-1));
    end

    dTdx(i, end) = (T(i, end) - T(i, end-1)) / (nodesx(end) - nodesx(end-1));
end


for j = 1:Nx+2
    dTdy(1, j) = (T(2, j) - T(1, j)) / (nodesy(2) - nodesy(1));

    for i = 2:Ny+1
        dTdy(i, j) = (T(i+1, j) - T(i-1, j)) / (nodesy(i+1) - nodesy(i-1));
    end

    dTdy(end, j) = (T(end, j) - T(end-1, j)) / (nodesy(end) - nodesy(end-1));
end

K = 16 .* (nodesy / Ly) + 16;
qx = -K .* dTdx;
qy = -K .* dTdy;
%% Plots
% Convergence Plot
figure;
plot(1:length(residuals), residuals, '-o');
title('Convergence History');
xlabel('Iteration Number');
ylabel('Residual');
grid on;

% Heat-flux Vector plot
figure;
contourf(X_nodes, Y_nodes, T, 20, "LineStyle", 'none');
colorbar;
hold on; 

quiver(X_nodes, Y_nodes, qx, qy, 'k', 'AutoScaleFactor', 2);

hold off;
title('Vector plot of Heat Flux');
xlabel('x');
ylabel('y');
set(gca, 'YDir', 'normal');
%% Gauss Scidel Iterative Solver
function [T, cnvrg_val] = GaussScidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly, S, q)
    % Boundary Conditions for 1
    T(1, :) = 15;
    % Boundary Condition for 3
    T(end, :) = 10;
    for i = 2:Ny+1
        y = nodesy(i) - nodesy(1);
        k = 16*(y/Ly) + 16;

        % Boundary Condition for 2
        T(i, end) = 5 - 5*(y/Ly) + 15*sin(pi*y/Ly)  ;
    
        for j = fliplr(2:Nx+1)
            [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S);

            TW = T(i, j-1);
            TE = T(i, j+1);
            TS = T(i-1, j);
            TN = T(i+1, j);
            
            TP = (aW*TW + aE*TE + aS*TS + aN*TN + Su) / aP;
            T(i, j) = TP;

            % Boundary Condition for 4
            if j == 2
                T(i, 1) = T(i, 2) - (q/k) * (nodesx(2) - nodesx(1));
            end

        end
    end

    cnvrg_val = calcResErr(T, nodesx, nodesy, Nx, Ny, dx, dy, Lx, Ly, S);
end
%% Residual Error Calculation
function cnvrg_val = calcResErr(T, nodesx, nodesy, Nx, Ny, dx, dy, Lx, Ly, S)
    ResErr_Matrx = zeros(Nx+2, Ny+2);
    F = abs(S * Lx * Ly);

    for i = 2:Ny+1
        y = nodesy(i) - nodesy(1);
        k = 16*(y/Ly) + 16;

        for j = 2:Nx+1
            [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S);
            
            TW = T(i, j-1);
            TE = T(i, j+1);
            TS = T(i-1, j);
            TN = T(i+1, j);
            TP = T(i, j);

            R = abs((aW*TW + aE*TE + aS*TS + aN*TN + Su) - aP*TP); % Residue Error
            ResErr_Matrx(i, j) = R/F; 
        end
    end

    cnvrg_val = sum(ResErr_Matrx(:));
end
%% Function to get coefficients for a given cell
function [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S)
    delx = dx(j);
    dely = dy(i);

    distx_back = nodesx(j) - nodesx(j-1);
    distx_front = nodesx(j+1) - nodesx(j);
    disty_down = nodesy(i) - nodesy(i-1);  
    disty_up = nodesy(i+1) - nodesy(i);

    aW = k * dely / distx_back;
    aE = k * dely / distx_front;
    aS = k * delx / disty_down;
    aN = k * delx / disty_up;
    aP = aW + aE + aS + aN;

    Su = S * delx * dely;
end
%% Create Mesh
function [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly)
    start = 0;
    
    strx = 0;
    stry = 0;
    dx = zeros(1,Nx);
    dy = zeros(1,Ny);
    
    %dx(1) = strx;
    %dy(1) = stry;
    
    dx(1) = Lx/Nx;
    dy(1) = Ly/Ny;
    
    x_line = zeros(1,Nx+1);
    y_line = zeros(1,Ny+1);
    nodesx = zeros(1,Nx+2);
    nodesy = zeros(1,Ny+2);
    
    nodesx(1) = start;
    nodesy(1) = start;
    
    x_line(1) = dx(1);
    y_line(1) = dx(1);
    
    for l = 1:Nx
        dx(l+1) = dx(l) - dx(l)*strx; % substitute the minus with plus for stretching in opposite direction
        x_line(l+1) = x_line(l) + dx(l);
    end 
    
    for k = 1:Ny
        dy(k+1) = dy(k) + dy(k)*stry; % substitute the plus with minus for stretching in opposite direction
        y_line(k+1) = y_line(k) + dy(k);
     end
    
    mx = min(x_line(:)); 
    my = min(y_line(:));
    rangex = max(x_line(:)) - mx;
    rangey = max(y_line(:)) - my;
    x_linee = (x_line - mx) / rangex;
    y_linee = (y_line - my) / rangey;
    range_x = Lx - start;
    range_y = Ly - start;
    x_lineee = (x_linee * range_x) + start;
    y_lineee = (y_linee * range_y) + start;
    
    for  g = 1:Nx
        nodesx(g+1) = (x_lineee(g) + x_lineee(g+1)) / 2;
        nodesy(g+1) = (y_lineee(g) + y_lineee(g+1)) / 2;
    end
    
    nodesx(Nx+2) = x_lineee(Nx+1);
    nodesy(Nx+2) = y_lineee(Nx+1);
end




