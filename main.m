clear; close all;

Lx = 1;    % Length of the mesh in x-direction
Ly = 0.5;  % Length of the mesh in y-direction

fprintf(['Plot the mesh: Section no.1 \nGauss Seidel Solver: Section no.2 \n' ...
    'Prove mesh independence: Section no.3 \nEffect of different error tolerance: Section no.4 \n' ...
    'Plot Heat flow (vector plot): Section no.5 \nChange BCs: Section no.6 \n'])
section = input('What do you want to do? Pls enter the Section No.');

if section ~= 3
    Nx = input('Enter the no.of cells you want in x-direction: ');
    Ny = input('Enter the no.of cells you want in y-direction: ');
    str_x = input('Enter the stretch(%) you would like in x-direction: ');
    str_y = input('Enter the stretch(%) you would like in y-direction: ');
end

if section ~= 1 && section ~= 4
    epsilon = input('Enter error tolerance you would like to set: ');
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
    
    T = zeros(Ny+2, Nx+2);
    max_iter = 100000;
    for n = 1:max_iter
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
    
    max_iter = 100000;
    for l = 1:length(Nx)
        [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx(l), Ny(l), Lx, Ly, str_x, str_y);
        fprintf('Solving for %dx%d \n', Nx(l), Ny(l))
        T = zeros(Ny(l)+2, Nx(l)+2);
    
        for n = 1:max_iter
            [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx(l), Ny(l), Lx, Ly);
            if cnvrg_val < epsilon
                break
            end
        end
    
        [X, Y] = meshgrid(nodesx, nodesy);
    
        avgT_LeftSec(l)   = mean(T(X <= Lx/3));
        avgT_MiddleSec(l) = mean(T(X > Lx/3 & X <= 2*Lx/3));
        avgT_RightSec(l)  = mean(T(X > 2*Lx/3));
        disp('Done')
    end
    
    fprintf('Avg Temp in left section with variation is %.2f(+-)%.2f \n', mean(avgT_LeftSec), std(avgT_LeftSec))
    fprintf('Avg Temp in middle section with variation is %.2f(+-)%.2f \n', mean(avgT_MiddleSec), std(avgT_MiddleSec))
    fprintf('Avg Temp in right section with variation is %.2f(+-)%.2f \n', mean(avgT_RightSec), std(avgT_RightSec))
end
%% Effect of different error tolerance (Section-4)
if section == 4
    figure
    colors = {'r', 'g', 'b'};              
    markers = {'-o', '-s', '-^'};         

    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, str_x, str_y);
    epsilon = [1e-4, 1e-5, 1e-6];
    T_dict = cell(1, length(epsilon));
    max_iter = 100000;
    for i = 1:length(epsilon)
        T = zeros(Ny+2, Nx+2);
        eps = epsilon(i);
        residuals = [];
        for n = 1:max_iter
            [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly);
            residuals(n) = cnvrg_val;
        
            if cnvrg_val < eps
                disp('Convergence has been achieved')
                fprintf('No.of iterations taken for eps = %.0e is: %d \n', eps, n)
                T_dict{i, 1} = T;
                break
            end
        end
 
         plot(1:length(residuals), residuals, markers{i}, 'Color', colors{i}, ...
             'DisplayName', sprintf('eps = %.0e', eps));
        hold on
    end
    xlabel('Iteration Number');
    ylabel('Residual');
    legend show
    grid on;
    
    err_norm_1 = norm(T_dict{1}(:) - T_dict{3}(:), 2);
    err_norm_2 = norm(T_dict{2}(:) - T_dict{3}(:), 2);   % L2 norm
    fprintf('L2 norm difference (eps=1e-4 vs 1e-6): %.6e\n', err_norm_1);
    fprintf('L2 norm difference (eps=1e-5 vs 1e-6): %.6e\n', err_norm_2);
end
%% Plotting Heat Flow (Section-5)
if section == 5
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, str_x, str_y);
    
    T = zeros(Ny+2, Nx+2);
    max_iter = 100000;
    for n = 1:max_iter
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
    contourf(X_nodes, Y_nodes, T, 30, "LineStyle", 'none'); 
    c = colorbar;
    c.Label.String = 'Temperature';
    c.Label.FontSize = 12;
    hold on;
    
    quiver(X_nodes, Y_nodes, qx, qy, 'k', 'AutoScaleFactor', 2);
    
    hold off;
    ax = gca;
    ax.FontSize = 14;
    xlabel('x', 'FontSize', 20);
    ylabel('y', 'FontSize', 20);
    set(gca, 'YDir', 'normal');
end
%% Changing boundary conditions (Section-6)
if section == 6
    [x_lineee, y_lineee, dx, dy, nodesx, nodesy] = generateMesh(Nx, Ny, Lx, Ly, str_x, str_y);
    
    T_dict = cell(1, 2);
    for a = 1:2
        T = zeros(Ny+2, Nx+2);
        max_iter = 100000;
        for n = 1:max_iter
            if a == 1
                [T, cnvrg_val] = GaussSeidelSolvernew(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly);
                T_dict{a} = T;
            else
                [T, cnvrg_val] = GaussSeidelSolver(T, nodesx, nodesy, dx, dy, Nx, Ny, Lx, Ly);
                T_dict{a} = T;
            end

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
        contourf(X_nodes, Y_nodes, T, 30, "LineStyle", 'none');
        c = colorbar;
        c.Label.String = 'Temperature';
        c.Label.FontSize = 12;
        hold on;
        
        quiver(X_nodes, Y_nodes, qx, qy, 'k', 'AutoScaleFactor', 2);
        
        hold off;
        ax = gca;
        ax.FontSize = 14;
        xlabel('x', 'FontSize', 20);
        ylabel('y', 'FontSize', 20);
        set(gca, 'YDir', 'normal');
    end

    [maxT1, linearIndex1] = max(T_dict{1}(:));
    [maxT2, linearIndex2] = max(T_dict{2}(:));   
    [row1, col1] = ind2sub(size(T_dict{1}), linearIndex1);
    [row2, col2] = ind2sub(size(T_dict{2}), linearIndex2);
    fprintf('Max Temp = %.4f at (row=%d, col=%d)\n', maxT1, Ny+1-row1, col1);
    fprintf('Max Temp = %.4f at (row=%d, col=%d)\n', maxT2, Ny+1-row2, col2);
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
        TP_old = T(i, j-1);
        [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S, TP_old);

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
                TP_old = T(i, j-1);
                [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S, TP_old);

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
%% Residual Error Calculation
function cnvrg_val = calcResErr(T, nodesx, nodesy, Nx, Ny, dx, dy, Lx, Ly, S)
ResErr_Matrx = zeros(Nx+2, Ny+2);
F = abs(S * Lx * Ly);

for i = 2:Ny+1
    y = nodesy(i) - nodesy(1);
    k = 16*(y/Ly) + 16;

    for j = 2:Nx+1
        TP_old = T(i, j-1);
        [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S, TP_old);

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
function [aW, aE, aS, aN, aP, Su] = getCoeffs(i, j, nodesx, nodesy, k, dx, dy, S, TP_old)
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

Su = 0;
Sp = (abs(S) * delx * dely) / TP_old;

aP = aW + aE + aS + aN - Sp;
end
%% Create Mesh
% To stretch/contract from the centre of the mesh in both directions
% Use (+ve) change to contract on both sides from the centre
% Use (-ve) change to stretch on both sides from the centre
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


