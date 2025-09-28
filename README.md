# Computational Fluid Dynamics (CFD)

### Finite Volume Method for 2D Heat Diffusion using Gauss-Seidel Iteration in MATLAB

**Authors**: Charukhesh B.R & Bhadra S

## 📖 Project Overview

This repository provides a MATLAB implementation of the **Finite Volume Method (FVM)** to solve the steady-state 2-dimensional heat diffusion equation. The numerical solution is obtained using the **Gauss-Seidel iterative method** on a structured, non-uniform (stretchable) mesh.

The project explores several key aspects of a CFD simulation, including mesh generation, solver convergence, mesh independence, and the physical interpretation of results under different boundary conditions.

### The Governing Equation

The simulation solves the following steady-state 2D heat diffusion equation with a variable thermal conductivity `k` and a source term `S`:

$$
\frac{\partial}{\partial x}\left(k \frac{\partial T}{\partial x}\right) + \frac{\partial}{\partial y}\left(k \frac{\partial T}{\partial y}\right) + S = 0
$$

### Problem Definition

The physical domain is a rectangle of length `Lx = 1.0` and height `Ly = 0.5`. The thermal conductivity `k` is a function of the vertical position `y`, and a constant source term `S` is applied.

- **Thermal Conductivity:** `k = 16 * (1 + y/Ly)`
- **Source Term:** `S = -1.5`

The simulation primarily focuses on the following set of boundary conditions:

| Boundary | Location | Type      | Condition                                    |
| :------- | :------- | :-------- | :------------------------------------------- |
| Bottom   | `y = 0`  | Dirichlet | `T = 15`                                     |
| Right    | `x = Lx` | Dirichlet | `T = 5(1 - y/Ly) + 15sin(πy/Ly)`               |
| Top      | `y = Ly` | Dirichlet | `T = 10`                                     |
| Left     | `x = 0`  | Neumann   | Heat Flux `q = -5000` (`-k * ∂T/∂x = -5000`) |

## 📁 Repository Structure

The project is contained within a single main script with all necessary functions defined internally.

```
FVM_DiffusionEqn/
│
├── main.m
└── README.md
```

- **main.m**: The main MATLAB script containing all the code. It is interactive and allows the user to select from various analysis modes, including mesh generation, solving, convergence studies, and visualization.

## 🚀 How to Run and Use This Project

1.  **Clone the Repository**
    ```bash
    git clone https://github.com/Charukhesh/FVM_DiffusionEqn.git
    cd FVM_DiffusionEqn
    ```

2.  **Open and Run the File in MATLAB**
    -   Open `main.m` in MATLAB.
    -   Run the script. You will be prompted in the command window to choose an analysis to perform.

### Interactive Analysis Sections

The `main.m` script is divided into sections. When you run the file, you will first see this prompt:

```
Plot the mesh: Section no.1
Gauss Seidel Solver: Section no.2
Prove mesh independence: Section no.3
Effect of different error tolerance: Section no.4
Plot Heat flow (vector plot): Section no.5
Change BCs: Section no.6
What do you want to do? Pls enter the Section No.
```

Below is a guide to each section and the inputs required.

#### **Section 1: Plot the Mesh**
-   **Description**: Generates and visualizes the computational mesh. This is useful for seeing the effect of mesh stretching.
-   **Inputs**:
    -   `Section No.`: **1**
    -   `No. of cells in x-direction (Nx)`: e.g., `40`
    -   `No. of cells in y-direction (Ny)`: e.g., `40`
    -   `Stretch (%) in x-direction`: Positive values contract the mesh towards the center, refining it near the left/right boundaries. e.g., `25`
    -   `Stretch (%) in y-direction`: Positive values contract the mesh, refining it near the top/bottom boundaries. e.g., `15`

---

#### **Section 2: Run the Gauss-Seidel Solver**
-   **Description**: Solves the heat diffusion equation for a single mesh configuration and plots the convergence history (residual vs. iteration number).
-   **Inputs**:
    -   `Section No.`: **2**
    -   `Nx`, `Ny`, `Stretch %`: As described above.
    -   `Error tolerance (epsilon)`: The convergence criterion. The solver stops when the residual drops below this value. e.g., `1e-6`

---

#### **Section 3: Prove Mesh Independence**
-   **Description**: Automatically runs the simulation on a series of increasingly finer grids (`10x10`, `20x20`, `40x40`, `80x80`) to demonstrate that the solution becomes independent of the mesh resolution. It calculates and compares the average temperatures in the left, middle, and right thirds of the domain.
-   **Inputs**:
    -   `Section No.`: **3**
    -   `Error tolerance (epsilon)`: e.g., `1e-6`

---

#### **Section 4: Effect of Different Error Tolerances**
-   **Description**: Compares the convergence behavior and final solutions for different error tolerances (`1e-4`, `1e-5`, `1e-6`). It plots the convergence histories on a single graph and calculates the L2 norm of the difference between the resulting temperature fields.
-   **Inputs**:
    -   `Section No.`: **4**
    -   `Nx`, `Ny`, `Stretch %`: As described above.

---

#### **Section 5: Plot Heat Flow (Vector Plot)**
-   **Description**: Solves the system and then generates a detailed visualization, showing the temperature field as a contour plot with the heat flux vectors (`qx`, `qy`) overlaid. This is the best way to visualize the physical results.
-   **Example Usage**: To generate the heat flow plot for a `40x40` uniform grid with a tolerance of `1e-6`:
    1.  Enter `5` for the section number.
    2.  Enter `40` for `Nx`.
    3.  Enter `40` for `Ny`.
    4.  Enter `0` for x-stretch and `0` for y-stretch.
    5.  Enter `1e-6` for the error tolerance.

---

#### **Section 6: Analyze Changed Boundary Conditions**
-   **Description**: This section demonstrates the model's flexibility. It runs the simulation twice: once with the original boundary conditions, and a second time with the top (Dirichlet) and left (Neumann) boundary conditions swapped. It then generates heat flow plots for both cases to compare the results.
-   **Inputs**:
    -   `Section No.`: **6**
    -   `Nx`, `Ny`, `Stretch %`: As described above.
    -   `Error tolerance (epsilon)`: e.g., `1e-6`

## 📊 Key Analyses and Insights

This project demonstrates several important CFD concepts:
1.  **Mesh Independence**: The study in `Section 3` confirms that as the grid is refined, the solution converges to a stable result. For this problem, the average temperatures vary by less than 1.2 K across different grid densities, indicating that a moderately fine mesh (e.g., 40x40) is sufficient.

2.  **Mesh Stretching for High-Gradient Regions**: The code allows for mesh refinement near boundaries (`Section 1`). This is crucial for accurately capturing the steep temperature gradients caused by the Neumann and complex Dirichlet boundary conditions without needing an excessively fine mesh throughout the domain.

3.  **Solver Convergence**: `Section 4` shows that a stricter error tolerance (e.g., `1e-6` vs. `1e-4`) leads to a more accurate solution at the cost of more iterations. However, the analysis reveals that the improvement in the solution (measured by the L2 norm) is minimal for tolerances beyond `1e-5`, allowing for a trade-off between accuracy and computational cost.

4.  **Physical Interpretation of Heat Flow**: The vector plots generated in `Section 5` and `Section 6` clearly visualize the physics. Heat enters from the hot bottom boundary and the high-flux left boundary, and it exits through the colder top and right boundaries. The heat flux vectors naturally curve upwards, following the path of higher thermal conductivity.
    
