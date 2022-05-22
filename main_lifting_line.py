import matplotlib.pyplot as plt
import numpy as np
from Lifting_Line import Turbine, GetRelaxationFactor
from read_write import read_from_file, write_to_file


def check_convergence(check_input: int):
    n_rot_wake, n_point_per_rotation, n_blade_elements, convection_speed = 8, 12, 30, 10
    [[cp_BEM, _]] = read_from_file(f'./saved_data/BEM_cp_cT_tsr_10.txt')

    if check_input == 0:
        n_lst = np.arange(2, 11)
    elif check_input == 1:
        n_lst = 3 * np.arange(1, 6)
    elif check_input == 2:
        n_lst = 5 * np.arange(4, 11)
    else:
        raise ValueError("Invalid input number.")

    delta = np.empty(n_lst.shape)

    for i, n in enumerate(n_lst):
        turbineParams = [n_rot_wake, n_point_per_rotation, n_blade_elements, convection_speed]
        turbineParams[check_input] = n
        turb = run_lifting_line(turbineParams=turbineParams)
        delta[i] = abs(cp_BEM - turb.extract_information_N_write(write=False)[2])
        print(delta[i])

    plt.plot(n_lst, delta)
    plt.xlabel('N (-)')
    plt.ylabel('$\\Delta C_p$ (-)')
    plt.show()


def compare_to_BEM():
    linestyles = ('dashed', 'solid', 'dotted')
    for j, tsr in enumerate([10]):#(6, 8, 10)):
        [r_BEM, alpha_BEM, phi_BEM, pn_BEM, pt_BEM, a_BEM, a_prime_BEM] = read_from_file('./saved_data/BEM_r_alpha_phi_pn_pt_a_aprime_tsr_%d.txt'%(tsr))
        [[cp_BEM, cT_BEM]] = read_from_file('./saved_data/BEM_cp_cT_tsr_%d.txt'%(tsr))
        [r_LL, alpha_LL, phi_LL, pn_LL, pt_LL, a_LL, a_prime_LL] = read_from_file('./saved_data/LL_r_alpha_phi_pn_pt_a_aprime_tsr_%d_rotorNumber_0.txt' % (tsr))
        [[cT_LL, cp_LL]] = read_from_file('./saved_data/LL_cp_cT_tsr_%d.txt' % (tsr))

        plt.figure(1, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Angle of attack ($\\alpha$) [$^{\\circ}$]')
        plt.plot(r_BEM, alpha_BEM, linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, np.degrees(alpha_LL), linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

        plt.figure(2, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Inflow angle ($\\phi$) [$^{\\circ}$]')
        plt.plot(r_BEM, np.degrees(phi_BEM), linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, np.degrees(phi_LL), linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

        plt.figure(3, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Normal force ($p_{n}$) [N/m]')
        plt.plot(r_BEM, pn_BEM, linestyle=linestyles[j], color='tab:blue', label='BEM ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_BEM))
        plt.plot(r_LL, pn_LL, linestyle=linestyles[j], color='tab:orange', label='Lifting line ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_LL))

        plt.figure(4, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Tangential force ($p_{t}$) [N/m]')
        plt.plot(r_BEM, pt_BEM, linestyle=linestyles[j], color='tab:blue', label='BEM ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_BEM))
        plt.plot(r_LL, pt_LL, linestyle=linestyles[j], color='tab:orange', label='Lifting line ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_LL))

        plt.figure(5, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Axial induction factor ($a$) [-]')
        plt.plot(r_BEM, a_BEM, linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, a_LL, linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

        plt.figure(6, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Tangential induction factor ($a^\prime$) [-]')
        plt.plot(r_BEM, a_prime_BEM, linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, a_prime_LL, linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

    plt.figure(1, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Angle_of_attack.pdf')

    plt.figure(2, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Inflow_angle.pdf')

    plt.figure(3, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Normal_force.pdf')

    plt.figure(4, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Tangential_force.pdf')

    plt.figure(5, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Axial_induction_factor.pdf')

    plt.figure(6, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Tangential_induction_factor.pdf')

    plt.show()


def compare_w_wo_induction():
    linestyles = ('dashed', 'solid', 'dotted')
    for j, tsr in enumerate([10]):#(6, 8, 10)):
        [r_LL_a, alpha_LL_a, phi_LL_a, pn_LL_a, pt_LL_a, a_LL_a, a_prime_LL_a] = read_from_file('./saved_data/LL_r_alpha_phi_pn_pt_a_aprime_tsr_%d_rotorNumber_0_induction.txt'%(tsr))
        [[cT_LL_a, cp_LL_a]] = read_from_file('./saved_data/LL_cp_cT_tsr_%d_induction.txt'%(tsr))
        [r_LL, alpha_LL, phi_LL, pn_LL, pt_LL, a_LL, a_prime_LL] = read_from_file('./saved_data/LL_r_alpha_phi_pn_pt_a_aprime_tsr_%d_rotorNumber_0.txt' % (tsr))
        [[cT_LL, cp_LL]] = read_from_file('./saved_data/LL_cp_cT_tsr_%d.txt' % (tsr))

        plt.figure(1, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Angle of attack ($\\alpha$) [$^{\\circ}$]')
        plt.plot(r_LL_a, np.degrees(alpha_LL_a), linestyle=linestyles[j+1], color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
        plt.plot(r_LL, np.degrees(alpha_LL), linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

        plt.figure(2, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Inflow angle ($\\phi$) [$^{\\circ}$]')
        plt.plot(r_LL_a, np.degrees(phi_LL_a), linestyle=linestyles[j+1], color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
        plt.plot(r_LL, np.degrees(phi_LL), linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

        plt.figure(3, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Normal force ($p_{n}$) [N/m]')
        plt.plot(r_LL_a, pn_LL_a, linestyle=linestyles[j+1], color='tab:blue', label='Lifting Line ($a_w$ = 0.16) ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_LL_a))
        plt.plot(r_LL, pn_LL, linestyle=linestyles[j], color='tab:orange', label='Lifting line ($a_w$ = 0) ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_LL))

        plt.figure(4, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Tangential force ($p_{t}$) [N/m]')
        plt.plot(r_LL_a, pt_LL_a, linestyle=linestyles[j+1], color='tab:blue', label='Lifting Line ($a_w$ = 0.16) ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_LL_a))
        plt.plot(r_LL, pt_LL, linestyle=linestyles[j], color='tab:orange', label='Lifting line ($a_w$ = 0) ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_LL))

        plt.figure(5, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Axial induction factor ($a$) [-]')
        plt.plot(r_LL_a, a_LL_a, linestyle=linestyles[j+1], color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
        plt.plot(r_LL, a_LL, linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

        plt.figure(6, figsize=(5, 5)); plt.xlabel('$r$ [m]'); plt.ylabel('Tangential induction factor ($a^\prime$) [-]')
        plt.plot(r_LL_a, a_prime_LL_a, linestyle=linestyles[j+1], color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
        plt.plot(r_LL, a_prime_LL, linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

    plt.figure(1, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Angle_of_attack.pdf')

    plt.figure(2, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Inflow_angle.pdf')

    plt.figure(3, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Normal_force.pdf')

    plt.figure(4, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Tangential_force.pdf')

    plt.figure(5, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Axial_induction_factor.pdf')

    plt.figure(6, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Tangential_induction_factor.pdf')

    plt.show()


def run_lifting_line(turbineParams=None, relaxFactor=None, multirotor=False):
    # Ugly but it worksss
    distance_between_turbines = 30

    if turbineParams is None:
        mainTurbine = Turbine(rotation=0, referencePos=(0,0,0))
        if multirotor:
            secondaryTurbine = Turbine(rotation=0, referencePos=(0,distance_between_turbines,0))

    else:
        mainTurbine = Turbine(*turbineParams, rotation=0, referencePos=(0,0,0))
        if multirotor:
            secondaryTurbine = Turbine(*turbineParams, rotation=0, referencePos=(0,distance_between_turbines,0))

    if multirotor:
        turbines = [mainTurbine, secondaryTurbine]
    else:
        turbines = [mainTurbine]
    # This file contains the main iteration over two turbines.
    maxiterations = 100
    iteri = 1
    bConverged = False
    relaxationFactor = .5 if relaxFactor is None else relaxFactor
    highestIndex, highestDeltaGamma, dr = mainTurbine.set_circulations_horseshoes(relaxationFactor)
    if multirotor:
        _, _, _ = secondaryTurbine.set_circulations_horseshoes(relaxationFactor)
    print(f"Initial calculation.\tMax delta gamma = {round(highestDeltaGamma, 3)}.\tRelaxation = {round(relaxationFactor, 3)}")
    while iteri < maxiterations:
        relaxationFactor = GetRelaxationFactor(abs(highestDeltaGamma) / relaxationFactor / dr, dr) if relaxFactor is None else relaxFactor
        mainTurbine.SetInducedVelocityForHorseshoes([turbines])
        if multirotor:
            secondaryTurbine.SetInducedVelocityForHorseshoes([turbines])
        highestIndex, highestDeltaGamma, dr = mainTurbine.set_circulations_horseshoes(relaxationFactor)
        if multirotor:
            highestIndex2, highestDeltaGamma2, dr2 = secondaryTurbine.set_circulations_horseshoes(relaxationFactor)
            # ugly but it works for debugging
            if abs(highestDeltaGamma2) > abs(highestDeltaGamma):
                highestIndex = highestIndex2
                highestDeltaGamma = highestDeltaGamma2
                dr = dr2

        print(f"Finished iteration {iteri}.\tMax delta gamma = {round(highestDeltaGamma, 3)}.\tRelaxation = {round(relaxationFactor, 3)}")

        # Exit criterion based on change in delta gamma
        if abs(highestDeltaGamma) < relaxationFactor * 1e-1:
            bConverged = True
            break

        iteri += 1

    print("Done :D")

    return mainTurbine


if __name__ == '__main__':
    # Sensitivity Analyses
    check_convergence(0)
    check_convergence(1)
    check_convergence(2)

    # Run with default setting for comparison to the BEM results
    turbine = run_lifting_line()
    turbine.extract_information_N_write()
    compare_to_BEM()

    # Run with some induction in the wake
    turbine = run_lifting_line((8, 12, 30, 10 * (1-.1585)))
    turbine.extract_information_N_write(suffix='induction')
    compare_w_wo_induction()
