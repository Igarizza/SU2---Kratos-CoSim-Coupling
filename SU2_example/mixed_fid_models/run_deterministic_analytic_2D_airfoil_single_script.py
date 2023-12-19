from dataclasses import dataclass, field
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# In [1], an analytic, nonlinear model of a two dimensional airfoil in airflow as an example for an internally coupled system is presented.
# Here the model is implemented with the numerical values from [1]. 
# The numerical values for u and r were previously back-calculated with the function derive_u_r, since they are not given in the paper.
# 
# [1] SOBIESZCZANSKI SOBIESKI, J. (1990). Sensitivity of complex, internally coupled systems. AIAA Journal, 28(1), 153 160. doi:10.2514/3.10366

@dataclass
class AerodynamicModel():
    # Model of the aerodynamic part of the airfoil

    # Fixed parameters of the model
    B_in_cm: float = 100
    C_in_cm: float = 10
    u: float = 5.111912504285256453e+00
    r: float = 1.905696411667853196e+00
    theta_0_in_rad: float = 0.26

    # Calculate the output of the model (lift L)
    def calc_L_in_N(self, phi_in_rad, q_in_N_per_cm_squared, psi_in_rad):
        S = self.B_in_cm * self.C_in_cm

        q = q_in_N_per_cm_squared
        theta_0 = self.theta_0_in_rad

        u = self.u
        r = self.r

        phi = phi_in_rad
        psi = psi_in_rad

        theta = phi + psi

        C_L = u * theta + r * (1 - np.cos((np.pi / 2) * (theta / theta_0)))

        L_in_N = q * S * C_L

        return L_in_N
    
@dataclass 
class StructuralModel():
    # Model of the structural part of the airfoil

    # Fixed parameters of the model
    a_overline: float = 0.25
    C_in_cm: float = 10
    k_1_in_N_per_cm: float = 4000
    k_2_in_N_per_cm: float = 2000
    z_overline_1: float = 0.2
    z_overline_2: float = 0.7

    @property
    def p(self):
        h_1_overline = self.a_overline - self.z_overline_1
        h_2_overline = self.z_overline_2 - self.a_overline

        return h_1_overline / h_2_overline

    # Calculate the output of the model (elastic support pitch angle phi)
    def calc_phi_in_rad(self, L_in_N):
        k_1 = self.k_1_in_N_per_cm
        k_2 = self.k_2_in_N_per_cm
        z_overline_1 = self.z_overline_1
        z_overline_2 = self.z_overline_2
        p = self.p
        C = self.C_in_cm

        R_1 = L_in_N / (1 + p)
        R_2 = L_in_N * p / (1 + p)

        d_1 = R_1 / k_1
        d_2 = R_2 / k_2

        phi_in_rad = (d_1 - d_2) / (C * (z_overline_2 - z_overline_1))

        return phi_in_rad

@dataclass
class AirfoilModel():
    # Model of the airfoil. It has two parts: An aerodynamic model and a structural model

    aerodynamic_model = AerodynamicModel()
    structural_model = StructuralModel()
    
    def calc_L_and_phi_gauss_seidel(self, proposal, airfoil_model_input, error_threshold):
        # Gauss–Seidel method
        #
        # proposal:             An array containing both proposals for the coupling variables and the deterministic input.
        # error_threshold:      The execution stops when the maximum absolute deviation from the previous iteration step is smaller than "error_threshold".

        print("--- Gauss–Seidel started ---")
        
        solution = [np.copy(proposal)]
        error = np.array([float('inf'), float('inf')])

        while np.any(error > error_threshold):
            solution_old = np.copy(solution[-1])
            solution.append(np.copy(solution[-1]))

            solution[-1][:, 0] = self.aerodynamic_model.calc_L_in_N(solution[-1][:, 1], float(airfoil_model_input["q_in_N_per_cm_squared"]), float(airfoil_model_input["psi_in_rad"]))
            solution[-1][:, 1] = self.structural_model.calc_phi_in_rad(solution[-1][:, 0])

            solution_new = np.copy(solution[-1])
            error = np.squeeze(np.abs(solution_new - solution_old))

            print(f"Absolute deviation from the previous iteration step: {error}")

        solution = np.squeeze(np.array(solution))

        print("--- Gauss–Seidel finished ---")

        self._write_to_df_csv(airfoil_model_input, solution)

        return solution
    
    def _write_to_df_csv(self, airfoil_model_input, solution):
        aerodynamic_model_input = pd.DataFrame(\
            {"q_in_N_per_cm_squared": np.repeat(float(airfoil_model_input["q_in_N_per_cm_squared"]), np.size(solution[:-1, 1])), \
             "psi_in_rad":  np.repeat(float(airfoil_model_input["psi_in_rad"]), np.size(solution[:-1, 1])), \
                "phi_in_rad": solution[:-1, 1]})
        aerodynamic_model_input.index.name='model execution no'
        structural_model_input = pd.DataFrame({"L_in_N": solution[1:, 0]})
        structural_model_input.index.name='model execution no'
        aerodynamic_model_output = pd.DataFrame({"L_in_N": solution[1:, 0]})
        aerodynamic_model_output.index.name='model execution no'
        structural_model_output = pd.DataFrame({"phi_in_rad": solution[1:, 1]})
        structural_model_output.index.name='model execution no'
        airfoil_model_output = pd.DataFrame({"L_in_N": solution[-1, 0]}, index=[0])
        airfoil_model_output.index.name='model execution no'

        airfoil_model_input.to_csv('csv/airfoil_model_input.csv')
        aerodynamic_model_input.to_csv('csv/aerodynamic_model_input.csv')
        structural_model_input.to_csv('csv/structural_model_input.csv')
        aerodynamic_model_output.to_csv('csv/aerodynamic_model_output.csv')
        structural_model_output.to_csv('csv/structural_model_output.csv')
        airfoil_model_output.to_csv('csv/airfoil_model_output.csv')

def main():
    pd.options.display.float_format = '{:,.20f}'.format

    proposal = np.array([[1e3, 1e-2]])
    airfoil_model_input = pd.DataFrame({"q_in_N_per_cm_squared": 1, "psi_in_rad": 0.05}, index=[0])
    airfoil_model_input.index.name='model execution no'
    error_threshold = np.array([1e-8, 1e-12])
    airfoil_model = AirfoilModel()
    sol_gauss_seidel = airfoil_model.calc_L_and_phi_gauss_seidel(proposal, airfoil_model_input, error_threshold)
    
    matplotlib.rcParams.update({'font.size': 18})

    plt.figure()
    plt.plot(sol_gauss_seidel[:, 0], color = 'r', label="$L^{(n)}$ [N]")
    plt.xlabel('Iteration $n$')
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(sol_gauss_seidel[:, 1], color = 'b', label="$\phi^{(n)}$ [rad]")
    plt.xlabel('Iteration $n$')
    plt.legend()
    plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    main()