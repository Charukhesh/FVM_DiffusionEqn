clear; close all;

Lx = 1;    % Length of the mesh in x-direction
Ly = 0.5;  % Length of the mesh in y-direction

fprintf(['Plot the mesh: Section no.1 \nGauss Seidel Solver: Section no.2 \n' ...
    'Prove mesh independence: Section no.3 \nPlot Heat flow (vector plot): Section no.4 \n'])
section = input('What do you want to do? Pls enter the Section No. ');

if section ~= 3
    Nx = input('Enter the no.of cells you want in x-direction: ');
    Ny = input('Enter the no.of cells you want in y-direction: ');
    str_x = input('Enter the stretch(%) you would like in x-direction: ');
    str_y = input('Enter the stretch(%) you would like in y-direction:  ');
end
%% Plotting the mesh (Section-1)
if section == 1
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, str_x, str_y);
    
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
end
%% Solving using Gauss-Seidel Method (Section-2) 
if section == 2
    
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, str_x, str_y);
    epsilon = 1e-4; % Error Tolerance for convergence
    
    T = zeros(Ny+2, Nx+2);
    num_iter = 100000;
    for n = 1:num_iter
        fprintf('Iteration no. %d \n', n)
        [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly);
        residuals(n) = cnvrg_val;
        disp(cnvrg_val)
    
        if cnvrg_val < epsilon
            disp('Convergence has been achieved')
            break
        else
            disp('Not yet, try harder')
        end
    end

    % Convergence Plot
    figure;
    plot(1:length(residuals), residuals, '-o');
    title('Convergence History');
    xlabel('Iteration Number');
    ylabel('Residual');
    grid on;
end
%% Checking mesh independence (Section-3)
if section == 3
    Nx = [10, 20, 40, 80];
    Ny = [10, 20, 40, 80];
    str_x = 0; str_y = 0;
    
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx(l), Ny(l), Lx, Ly, str_x, str_y);
    epsilon = 1e-4; % Error Tolerance for convergence
    num_iter = 100000;
    for l = 1:length(Nx)
        fprintf('Solving for %dx%d \n', Nx(l), Ny(l))
        T = zeros(Ny(l)+2, Nx(l)+2);
    
        for n = 1:num_iter
            [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx(l), Ny(l), Lx, Ly);
            if cnvrg_val < epsilon
                break
            else
                continue
            end
        end
    
        [X, Y] = meshgrid(nodesx, nodesy);
    
        avgT_LeftSec(l)   = mean(T(X <= Lx/3));
        avgT_MiddleSec(l) = mean(T(X > Lx/3 & X <= 2*Lx/3));
        avgT_RightSec(l)  = mean(T(X > 2*Lx/3));
        disp('Done')
    end
    
    fprintf('Std Deviation of Avg Temp in left section is %.2f \n', std(avgT_LeftSec))
    fprintf('Std Deviation of Avg Temp in middle section is %.2f \n', std(avgT_MiddleSec))
    fprintf('Std Deviation of Avg Temp in right section is %.2f \n', std(avgT_RightSec))
end
%% Plotting Heat Flow (Section-4)
if section == 4
    str_x = 25; str_y = 15;
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, str_x, str_y);
    epsilon = 1e-4; % Error Tolerance for convergence
    
    T = zeros(Ny+2, Nx+2);
    num_iter = 100000;
    for n = 1:num_iter
        [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly);   
        if cnvrg_val < epsilon
            break
        else
            continue
        end
    end

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
    
    % Heat-flux Vector plot
    figure;
    contourf(X_nodes, Y_nodes, T, 20, "LineStyle", 'none');
    c = colorbar;
    c.Label.String = 'Temperature';
    hold on;
    
    quiver(X_nodes, Y_nodes, qx, qy, 'k', 'AutoScaleFactor', 1.5);
    
    hold off;
    title('Vector plot of Heat Flux');
    xlabel('x');
    ylabel('y');
    set(gca, 'YDir', 'normal');
end
%% Gauss Scidel Iterative Solver (Case-4)
function [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly)
S = -1.5; q = -5000;
% Boundary Conditions for 1
T(1, :) = 15;
% Boundary Condition for 3
T(end, :) = 10;
for i = 2:Ny+1
    y = nodesy(i) - nodesy(1);
    k = 16*(y/Ly) + 16;

    % Boundary Condition for 2
    T(i, end) = 5 - 5*(y/Ly) + 15*sin(pi*y/Ly);

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
%% Gauss Scidel Iterative Solver (switching Neumann & Dirichlet Conditions)
function [T, cnvrg_val] = GaussSeidelSolvernew(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly)
        S = -1.5; q = -5000;
        % Boundary Conditions for 1
        T(1, :) = 15;
        % Boundary Condition for 4
        T(:, 1) = 20;
        for i = 2:Ny+1
            y = nodesy(i) - nodesy(1);
            k = 16*(y/Ly) + 16;

            % Boundary Condition for 2
            T(i, end) = 5 - 5*(y/Ly) + 15*sin(pi*y/Ly);
            for j = fliplr(2:Nx+1)
                [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S);

                TW = T(i, j-1);
                TE = T(i, j+1);
                TS = T(i-1, j);
                TN = T(i+1, j);

                TP = (aW*TW + aE*TE + aS*TS + aN*TN + Su) / aP;
                T(i, j) = TP;

                % Boundary Condition for 3
                if i == Nx+1
                    T(end, j) = T(end-1, j) - (q/k) * (nodesy(end) - nodesy(end-1));
                end

            end
        end

        cnvrg_val = calcResErr(T, nodesx, nodesy, Nx, Ny, dx, dy, Lx, Ly, S);
    end
%% Gauss Scidel Iterative Solver (Case-3)
    function [T, cnvrg_val] = GaussSeidelSolverCase3(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly)
        S = -1.5; q = 5000;
        % Boundary Conditions for 1
        T(1, :) = 15;
        % Boundary Condition for 3
        T(end, :) = 10;
        for i = 2:Ny+1
            y = nodesy(i) - nodesy(1);
            k = 16*(y/Ly) + 16;

            for j = 2:Nx+1
                [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S);

                TW = T(i, j-1);
                TE = T(i, j+1);
                TS = T(i-1, j);
                TN = T(i+1, j);

                TP = (aW*TW + aE*TE + aS*TS + aN*TN + Su) / aP;
                T(i, j) = TP;

                % Boundary Condition for 2
                if j == Nx+1
                    T(i, end) = T(i, end-1) - (q/k) * (nodesx(end) - nodesx(end-1));
                end

                % Boundary Condition for 4
                if j == 2
                    T(i, 1) = T(i, 2);
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
function [x_lines, y_lines, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, change_x, change_y)
    dx = Lx/Nx; 
    dy = Ly/Ny;
    Cx = ceil(Nx/2);
    Cy = ceil(Ny/2);
    x_lines = 0; y_lines = 0;
    
    for j = 1:Cx
        dx = [dx, (1 - change_x*0.01)*dx(1)];
        dx = [(1 - change_x*0.01)*dx(1), dx];
    end
    for l = 1:Nx
        x_lines(l+1) = x_lines(l) + dx(l);
    end
    x_lines = x_lines*Lx/x_lines(end);
    
    for r = 1:Cy
        dy = [dy, (1 - change_y*0.01)*dy(1)];
        dy = [(1 - change_y*0.01)*dy(1), dy];
    end
    for k = 1:Ny
        y_lines(k+1) = y_lines(k) + dy(k);
    end
    y_lines = y_lines*Ly/y_lines(end);
    
    nodesx = x_lines(1);
    for var1 = 1:length(x_lines)-1
        nodesx = [nodesx, (x_lines(var1+1) + x_lines(var1))*0.5];
    end
    nodesx = [nodesx, x_lines(end)];
    
    nodesy = y_lines(1);
    for var2 = 1:length(y_lines)-1
        nodesy = [nodesy, (y_lines(var2+1) + y_lines(var2))*0.5];
    end
    nodesy = [nodesy, y_lines(end)];
end


