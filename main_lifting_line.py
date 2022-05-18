import matplotlib.pyplot as plt
import numpy as np
from Lifting_Line import Turbine, GetRelaxationFactor, compare_to_BEM
from read_write import read_from_file, write_to_file


def check_convergence():
    n_rot_wake, n_point_per_rotation, n_blade_elements, convection_speed = 8, 12, 30, 10
    [[cp_BEM, _]] = read_from_file(f'./saved_data/BEM_cp_cT_tsr_10.txt')

    n_lst = 5 * np.arange(4, 11)
    delta = np.empty(n_lst.shape)

    for i, n in enumerate(n_lst):
        turbineParams = (n_rot_wake, n_point_per_rotation, n, convection_speed)
        turbine = run_lifting_line(turbineParams=turbineParams)
        delta[i] = abs(cp_BEM - turbine.extract_information_N_write(write=False)[2])
        print(delta[i])

    plt.plot(n_lst, delta)
    plt.xlabel('N (-)')
    plt.ylabel('$\\Delta C_p$ (-)')
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
    #
    # plt.figure(1)
    # out, _, _ = mainTurbine.extract_information_N_write()
    # plt.plot(out[0, 0], np.degrees(out[0, 1]), label='alpha 0')
    # plt.plot(out[0, 0], np.degrees(out[0, 2]), label='phi 0')
    # plt.plot(out[1, 0], np.degrees(out[1, 1]), label='alpha 1')
    # plt.plot(out[1, 0], np.degrees(out[1, 2]), label='phi 1')
    # plt.plot(out[2, 0], np.degrees(out[2, 1]), label='alpha 2')
    # plt.plot(out[2, 0], np.degrees(out[2, 2]), label='phi 2')
    # plt.legend()
    #
    # plt.figure(2)
    # plt.plot(out[0, 0], out[0, 3], label='pn 0')
    # plt.plot(out[0, 0], out[0, 4], label='pt 0')
    # plt.plot(out[1, 0], out[1, 3], label='pn 1')
    # plt.plot(out[1, 0], out[1, 4], label='pt 1')
    # plt.plot(out[2, 0], out[2, 3], label='pn 2')
    # plt.plot(out[2, 0], out[2, 4], label='pt 2')
    # plt.legend()
    #
    # plt.figure(3)
    # plt.plot(out[0, 0], out[0, 5], label='a 0')
    # plt.plot(out[0, 0], out[0, 6], label='ap 0')
    # plt.plot(out[1, 0], out[1, 5], label='a 1')
    # plt.plot(out[1, 0], out[1, 6], label='ap 1')
    # plt.plot(out[2, 0], out[2, 5], label='a 2')
    # plt.plot(out[2, 0], out[2, 6], label='ap 2')
    # plt.legend()
    #
    # plt.show()


if __name__ == '__main__':
    run_lifting_line()
    #compare_to_BEM()
