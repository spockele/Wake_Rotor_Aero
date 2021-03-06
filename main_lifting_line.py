import matplotlib.pyplot as plt
import numpy as np
from Lifting_Line import Turbine, GetRelaxationFactor
from read_write import read_from_file, write_to_file


def read_normalised(path):
    rho = 1.225
    v0 = 10
    r_blade = 50

    [r, alpha, phi, pn, pt, a, a_prime] = read_from_file(path)
    dr = abs(r[1] - r[0])

    pn = 3 * pn / (.5 * rho * v0 ** 2 * np.pi * r_blade)
    pt = 3 * pt / (.5 * rho * v0 ** 2 * np.pi * r_blade)

    r = r / r_blade

    return r, alpha, phi, pn, pt, a, a_prime


def check_convergence(check_input: int):
    n_rot_wake, n_point_per_rotation, n_blade_elements, convection_speed = 8, 12, 30, 10
    [[cp_BEM, _]] = read_from_file(f'./saved_data/BEM_cp_cT_10.txt')

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
        relax = .2 if check_input == 2 or i < 2 else None

        turb = run_lifting_line(turbineParams=turbineParams, relaxFactor=relax)
        delta[i] = abs(cp_BEM - turb.extract_information_N_write(write=False)[2])
        print(delta[i])

    plt.plot(n_lst, delta)
    plt.xlabel('N (-)')
    plt.ylabel('$\\Delta C_p$ (-)')
    plt.savefig(f'figures/Sensitivity_{check_input}.pdf')
    plt.show()


def compare_to_BEM():
    linestyles = ('dashed', 'solid', 'dotted')
    for j, tsr in enumerate((6, 8, 10)):
        [r_BEM, alpha_BEM, phi_BEM, pn_BEM, pt_BEM, a_BEM, a_prime_BEM] = read_normalised(f'./saved_data/BEM_output_{tsr}.txt')
        [[cp_BEM, cT_BEM]] = read_from_file(f'./saved_data/BEM_cp_cT_{tsr}.txt')
        [r_LL, alpha_LL, phi_LL, pn_LL, pt_LL, a_LL, a_prime_LL] = read_normalised(f'./saved_data/LL_output_{tsr}_blade0.txt')
        [[cT_LL, cp_LL]] = read_from_file(f'./saved_data/LL_cT_cp_{tsr}.txt')

        plt.figure(1, figsize=(5, 5))
        plt.plot(r_BEM, alpha_BEM, linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, np.degrees(alpha_LL), linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

        plt.figure(2, figsize=(5, 5))
        plt.plot(r_BEM, np.degrees(phi_BEM), linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, np.degrees(phi_LL), linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

        plt.figure(3, figsize=(5, 5))
        plt.plot(r_BEM, pn_BEM, linestyle=linestyles[j], color='tab:blue', label='BEM ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_BEM))
        plt.plot(r_LL, pn_LL, linestyle=linestyles[j], color='tab:orange', label='Lifting line ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_LL))

        plt.figure(4, figsize=(5, 5))
        plt.plot(r_BEM, pt_BEM, linestyle=linestyles[j], color='tab:blue', label='BEM ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_BEM))
        plt.plot(r_LL, pt_LL, linestyle=linestyles[j], color='tab:orange', label='Lifting line ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_LL))

        plt.figure(5, figsize=(5, 5))
        plt.plot(r_BEM, a_BEM, linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, a_LL, linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

        plt.figure(6, figsize=(5, 5))
        plt.plot(r_BEM, a_prime_BEM, linestyle=linestyles[j], color='tab:blue', label=f'BEM ($\\lambda={tsr}$)')
        plt.plot(r_LL, a_prime_LL, linestyle=linestyles[j], color='tab:orange', label=f'Lifting line ($\\lambda={tsr}$)')

    plt.figure(1, figsize=(5, 5))
    plt.xlabel('$r/R$ [-]')
    plt.ylabel('Angle of attack ($\\alpha$) [$^{\\circ}$]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/Compare_to_BEM/Angle_of_attack.pdf')

    plt.figure(2, figsize=(5, 5))
    plt.xlabel('$r/R$ [-]')
    plt.ylabel('Inflow angle ($\\phi$) [$^{\\circ}$]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/Compare_to_BEM/Inflow_angle.pdf')

    plt.figure(3, figsize=(5, 5))
    plt.xlabel('$r/R$ [-]')
    plt.ylabel('Normal force ($C_n$) [-]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/Compare_to_BEM/Normal_force.pdf')

    plt.figure(4, figsize=(5, 5))
    plt.xlabel('$r/R$ [-]')
    plt.ylabel('Tangential force ($C_t$) [-]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/Compare_to_BEM/Tangential_force.pdf')

    plt.figure(5, figsize=(5, 5))
    plt.xlabel('$r/R$ [-]')
    plt.ylabel('Axial induction factor ($a$) [-]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/Compare_to_BEM/Axial_induction_factor.pdf')

    plt.figure(6, figsize=(5, 5))
    plt.xlabel('$r/R$ [-]')
    plt.ylabel('Tangential induction factor ($a^\prime$) [-]')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./figures/Compare_to_BEM/Tangential_induction_factor.pdf')

    plt.show()


def compare_w_wo_induction():
    tsr = 10

    [r_LL_a, alpha_LL_a, phi_LL_a, pn_LL_a, pt_LL_a, a_LL_a, a_prime_LL_a] = read_normalised(f'./saved_data/LL_output_{tsr}_blade0_induction.txt')
    [[cT_LL_a, cp_LL_a]] = read_from_file(f'./saved_data/LL_cT_cp_{tsr}_induction.txt')
    [r_LL, alpha_LL, phi_LL, pn_LL, pt_LL, a_LL, a_prime_LL] = read_normalised(f'./saved_data/LL_output_{tsr}_blade0.txt')
    [[cT_LL, cp_LL]] = read_from_file(f'./saved_data/LL_cT_cp_{tsr}.txt')

    plt.figure(1, figsize=(5, 5)); plt.xlabel('$r/R$ [-]'); plt.ylabel('Angle of attack ($\\alpha$) [$^{\\circ}$]')
    plt.plot(r_LL_a, np.degrees(alpha_LL_a), color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
    plt.plot(r_LL, np.degrees(alpha_LL), color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

    plt.figure(2, figsize=(5, 5)); plt.xlabel('$r/R$ [-]'); plt.ylabel('Inflow angle ($\\phi$) [$^{\\circ}$]')
    plt.plot(r_LL_a, np.degrees(phi_LL_a), color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
    plt.plot(r_LL, np.degrees(phi_LL), color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

    plt.figure(3, figsize=(5, 5)); plt.xlabel('$r/R$ [-]'); plt.ylabel('Normal force ($C_n$) [-]')
    plt.plot(r_LL_a, pn_LL_a, color='tab:blue', label='Lifting Line ($a_w$ = 0.16) ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_LL_a))
    plt.plot(r_LL, pn_LL, color='tab:orange', label='Lifting line ($a_w$ = 0) ($\\lambda=%d$), $c_T=%.3f$'%(tsr, cT_LL))

    plt.figure(4, figsize=(5, 5)); plt.xlabel('$r/R$ [-]'); plt.ylabel('Tangential force ($C_t$) [-]')
    plt.plot(r_LL_a, pt_LL_a, color='tab:blue', label='Lifting Line ($a_w$ = 0.16) ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_LL_a))
    plt.plot(r_LL, pt_LL, color='tab:orange', label='Lifting line ($a_w$ = 0) ($\\lambda=%d$), $C_P=%.3f$'%(tsr, cp_LL))

    plt.figure(5, figsize=(5, 5)); plt.xlabel('$r/R$ [-]'); plt.ylabel('Axial induction factor ($a$) [-]')
    plt.plot(r_LL_a, a_LL_a, color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
    plt.plot(r_LL, a_LL, color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

    plt.figure(6, figsize=(5, 5)); plt.xlabel('$r/R$ [-]'); plt.ylabel('Tangential induction factor ($a^\prime$) [-]')
    plt.plot(r_LL_a, a_prime_LL_a, color='tab:blue', label=f'Lifting Line ($a_w$ = 0.16) ($\\lambda={tsr}$)')
    plt.plot(r_LL, a_prime_LL, color='tab:orange', label=f'Lifting line ($a_w$ = 0) ($\\lambda={tsr}$)')

    plt.figure(1, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Wake_Speed/Angle_of_attack.pdf')

    plt.figure(2, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Wake_Speed/Inflow_angle.pdf')

    plt.figure(3, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Wake_Speed/Normal_force.pdf')

    plt.figure(4, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Wake_Speed/Tangential_force.pdf')

    plt.figure(5, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Wake_Speed/Axial_induction_factor.pdf')

    plt.figure(6, figsize=(5, 5))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig('./figures/Wake_Speed/Tangential_induction_factor.pdf')

    plt.show()


def multirotor_phased():
    phases = 30 * np.arange(0, 4)
    distance = 1 * 100
    for phase in phases:
        turbs = run_lifting_line(multirotor=True, phase=np.radians(phase), distance=distance)
        out0, ct0, cp0 = turbs[0].extract_information_N_write(suffix=f'_turb0_phase{phase}')
        out1, ct1, cp1 = turbs[1].extract_information_N_write(suffix=f'_turb1_phase{phase}')


def compare_phases():
    phases = 30 * np.arange(0, 4)
    linestyles = ('dashed', 'solid', 'dotted')
    for phase in phases:
        for j in range(3):
            [r_turb0, alpha_turb0, phi_turb0, pn_turb0, pt_turb0, a_turb0, a_prime_turb0] = read_normalised(f'./saved_data/LL_output_10_blade{j}_turb0_phase{phase}.txt')
            [[cT_turb0, cp_turb0]] = read_from_file(f'./saved_data/LL_cT_cp_10_turb0_phase{phase}.txt')
            [r_turb1, alpha_turb1, phi_turb1, pn_turb1, pt_turb1, a_turb1, a_prime_turb1] = read_normalised(f'./saved_data/LL_output_10_blade{j}_turb1_phase{phase}.txt')
            [[cT_turb1, cp_turb1]] = read_from_file(f'./saved_data/LL_cT_cp_10_turb1_phase{phase}.txt')

            plt.figure(1, figsize=(5, 5))
            plt.plot(r_turb0, np.degrees(alpha_turb0), linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, np.degrees(alpha_turb1), linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(2, figsize=(5, 5))
            plt.plot(r_turb0, np.degrees(phi_turb0), linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, np.degrees(phi_turb1), linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(3, figsize=(5, 5))
            plt.plot(r_turb0, pn_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, pn_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(4, figsize=(5, 5))
            plt.plot(r_turb0, pt_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, pt_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(5, figsize=(5, 5))
            plt.plot(r_turb0, a_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, a_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(6, figsize=(5, 5))
            plt.plot(r_turb0, a_prime_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, a_prime_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

        plt.figure(1, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Angle of attack ($\\alpha$) [$^{\\circ}$]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_phases/Angle_of_attack_phase{phase}.pdf')

        plt.figure(2, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Inflow angle ($\\phi$) [$^{\\circ}$]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_phases/Inflow_angle_phase{phase}.pdf')

        plt.figure(3, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Normal force ($C_n$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_phases/Normal_force_phase{phase}.pdf')

        plt.figure(4, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Tangential force ($C_t$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_phases/Tangential_force_phase{phase}.pdf')

        plt.figure(5, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Axial induction factor ($a$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_phases/Axial_induction_factor_phase{phase}.pdf')

        plt.figure(6, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Tangential induction factor ($a^\prime$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_phases/Tangential_induction_factor_phase{phase}.pdf')

        plt.show()


def multirotor_spaced():
    distances = 100 * np.array([1, 2, 5, 1000])
    phase = 60
    for distance in distances:
        turbs = run_lifting_line(multirotor=True, phase=np.radians(phase), distance=distance)
        out0, ct0, cp0 = turbs[0].extract_information_N_write(suffix=f'_turb0_dist{distance}')
        out1, ct1, cp1 = turbs[1].extract_information_N_write(suffix=f'_turb1_dist{distance}')


def compare_distances():
    distances = 100 * np.array([1, 2, 5, 1000])
    linestyles = ('dashed', 'solid', 'dotted')
    for dist in distances:
        for j in range(3):
            [r_turb0, alpha_turb0, phi_turb0, pn_turb0, pt_turb0, a_turb0, a_prime_turb0] = read_normalised(
                f'./saved_data/LL_output_10_blade{j}_turb0_dist{dist}.txt')
            [[cT_turb0, cp_turb0]] = read_from_file(f'./saved_data/LL_cT_cp_10_turb0_dist{dist}.txt')
            [r_turb1, alpha_turb1, phi_turb1, pn_turb1, pt_turb1, a_turb1, a_prime_turb1] = read_normalised(
                f'./saved_data/LL_output_10_blade{j}_turb1_dist{dist}.txt')
            [[cT_turb1, cp_turb1]] = read_from_file(f'./saved_data/LL_cT_cp_10_turb1_dist{dist}.txt')

            plt.figure(1, figsize=(5, 5))
            plt.plot(r_turb0, np.degrees(alpha_turb0), linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, np.degrees(alpha_turb1), linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(2, figsize=(5, 5))
            plt.plot(r_turb0, np.degrees(phi_turb0), linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, np.degrees(phi_turb1), linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(3, figsize=(5, 5))
            plt.plot(r_turb0, pn_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, pn_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(4, figsize=(5, 5))
            plt.plot(r_turb0, pt_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, pt_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(5, figsize=(5, 5))
            plt.plot(r_turb0, a_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, a_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

            plt.figure(6, figsize=(5, 5))
            plt.plot(r_turb0, a_prime_turb0, linestyle=linestyles[j], color='tab:blue',
                     label=f'Turbine 0, blade {j}')
            plt.plot(r_turb1, a_prime_turb1, linestyle=linestyles[j], color='tab:orange',
                     label=f'Turbine 1, blade {j}')

        plt.figure(1, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Angle of attack ($\\alpha$) [$^{\\circ}$]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_distances/Angle_of_attack_dist{dist}.pdf')

        plt.figure(2, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Inflow angle ($\\phi$) [$^{\\circ}$]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_distances/Inflow_angle_dist{dist}.pdf')

        plt.figure(3, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Normal force ($C_n$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_distances/Normal_force_dist{dist}.pdf')

        plt.figure(4, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Tangential force ($C_t$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_distances/Tangential_force_dist{dist}.pdf')

        plt.figure(5, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Axial induction factor ($a$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_distances/Axial_induction_factor_dist{dist}.pdf')

        plt.figure(6, figsize=(5, 5))
        plt.xlabel('$r/R$ [-]')
        plt.ylabel('Tangential induction factor ($a^\prime$) [-]')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./figures/Compare_distances/Tangential_induction_factor_dist{dist}.pdf')

        plt.show()


def run_lifting_line(turbineParams=None, relaxFactor=None, tsr=None, multirotor=False, phase=None, distance=None):
    # Ugly but it worksss
    tsr = 10 if tsr is None else tsr
    distance_between_turbines = 100 if distance is None else distance
    phase = 0 if phase is None else phase
    relaxationFactor = .5 if relaxFactor is None else relaxFactor

    if turbineParams is None:
        mainTurbine = Turbine(tsr=tsr, rotation=0, referencePos=(0, 0, 0))
        turbines = [mainTurbine]
        if multirotor:
            secondaryTurbine = Turbine(tsr=tsr, rotation=phase, referencePos=(distance_between_turbines, 0, 0))
            turbines.append(secondaryTurbine)
            # relaxationFactor = .1 if relaxFactor is None else relaxFactor

    else:
        mainTurbine = Turbine(*turbineParams, tsr=tsr, rotation=0, referencePos=(0, 0, 0))
        turbines = [mainTurbine]
        if multirotor:
            secondaryTurbine = Turbine(*turbineParams, tsr=tsr, rotation=phase, referencePos=(distance_between_turbines, 0, 0))
            turbines.append(secondaryTurbine)
            # relaxationFactor = .1 if relaxFactor is None else relaxFactor

    maxiterations = 100
    iteri = 1
    bConverged = False

    highestIndex, highestDeltaGamma, dr = turbines[0].set_circulations_horseshoes(relaxationFactor)
    if multirotor:
        _, _, _ = turbines[1].set_circulations_horseshoes(relaxationFactor)
    print(f"Initial calculation.\tMax delta gamma = {round(highestDeltaGamma, 3)}.\tRelaxation = {round(relaxationFactor, 3)}\tIndex = {highestIndex}")

    while iteri < maxiterations:
        relaxationFactor = GetRelaxationFactor(abs(highestDeltaGamma) / relaxationFactor / dr, dr) if relaxFactor is None else relaxFactor

        turbines[0].SetInducedVelocityForHorseshoes(turbines)
        if multirotor:
            print(f'Turbine 1 induced')
            turbines[1].SetInducedVelocityForHorseshoes(turbines)
            print(f'Turbine 2 induced')

        highestIndex, highestDeltaGamma, dr = turbines[0].set_circulations_horseshoes(relaxationFactor)
        ind = 0

        if multirotor:
            highestIndex2, highestDeltaGamma2, dr2 = turbines[1].set_circulations_horseshoes(relaxationFactor)
            # ugly but it works for debugging
            if abs(highestDeltaGamma2) > abs(highestDeltaGamma):
                highestIndex = highestIndex2
                highestDeltaGamma = highestDeltaGamma2
                dr = dr2
                ind = 1

        print(f"Finished iteration {iteri}.\tMax delta gamma = {round(highestDeltaGamma, 3)}.\tRelaxation = {round(relaxationFactor, 3)}\tIndex = {ind}, {highestIndex}")

        # Exit criterion based on change in delta gamma
        if abs(highestDeltaGamma) < relaxationFactor * 1e-1:
            bConverged = True
            break

        iteri += 1

    print("Done :D")

    return turbines if multirotor else mainTurbine


if __name__ == '__main__':
    # Sensitivity Analyses
    print("--- Running Sensitivity Analysis ---")
    check_convergence(0)
    check_convergence(1)
    check_convergence(2)

    # Run with default setting for comparison to the BEM results
    print("--- Running with default settings ---")
    for lamda in (6, 8, 10):
        turbine = run_lifting_line(tsr=lamda)
        turbine.extract_information_N_write()

    # Run with some induction in the wake
    print("--- Running with induction in the wake ---")
    turbine = run_lifting_line((8, 12, 30, 10 * (1-.1585)))
    turbine.extract_information_N_write(suffix='_induction')
    compare_w_wo_induction()

    # Multirotor Stuff
    # Different phases between the 2 rotors
    print("--- Running multiple rotors at different phase angles ---")
    multirotor_phased()
    compare_phases()
    # Different distances between the 2 rotors
    print("--- Running multiple rotors at different rotor spacings ---")
    multirotor_spaced()
    compare_distances()
