clear; close all;

L = 1; H = 0.5; S = -1.5;

Nx = 10; Ny = 10;
del_x = L/Nx; del_y = H/Ny; 

x = L/2; y = H/2;
function T = GaussScidelSolver(x, y, del_x, del_y)
    TW = find_T_val(4, y); % Boundary 4
    TS = find_T_val(1, y); % Boundary 1
    TE = 0; TN = 0; % Initial Guess
    Su = S * del_x * del_y;
    for i = 1:Ny
        for j = 1:Nx
            aW = find_a_val(x-L, y, 1);
            aE = find_a_val(x+L, y, 1);
            aS = find_a_val(x, y-L, 0);
            aN = find_a_val(x, y+L, 0);
            aP = aW + aE + aS + aN;

            TP = (aW*TW + aE*TE + aS*TS + aN*TN + Su) / aP;
            x = x + L;
            TW = TP;
            TE = 0; % Guess
        end
        y = y+ H;
        TS = TP;
        TN = 0; % Guess
    end
    T = [TW, TE, TS, TN];
end

% 0 for 'n' or 's' and 1 for 'e' or 'w'
function a = find_a_val(x, y, side)
    if side == 0
        a = find_k_val(y) * del_x / x;
    elseif side == 1
        a = find_k_val(y) * del_y / x;
    end
end

function k = find_k_val(y)
    k = 16*(y/H) + 16; 
end

function T = find_T_val(bdy_no, y)
    if bdy_no == 1
        T = 15;
    elseif bdy_no == 2
        T = 5 - 5*(y/H) + 15*sin(pi*y/H);
    elseif bdy_no == 3
        T = 10;
    elseif bdy_no == 4
        q = -5000;
        dT_dx = q/find_k_val(y);
    end
end







