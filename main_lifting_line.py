import matplotlib.pyplot as plt
import numpy as np
from Lifting_Line import Turbine, GetRelaxationFactor


def check_convergence():
    pass


def run_lifting_line():
    # Ugly but it worksss

    mainTurbine = Turbine(0, (0,0,0))

    # This file contains the main iteration over two turbines.
    maxiterations = 100
    iteri = 1
    bConverged = False
    relaxationFactor = .5
    highestIndex, highestDeltaGamma = mainTurbine.set_circulations_horseshoes(relaxationFactor)
    print(f"Initial calculation. max delta gamma = {highestDeltaGamma} at idx = {highestIndex}")
    while iteri < maxiterations:
        relaxationFactor = GetRelaxationFactor(highestDeltaGamma/relaxationFactor)
        mainTurbine.SetInducedVelocityForHorseshoes()
        highestIndex, highestDeltaGamma = mainTurbine.set_circulations_horseshoes(relaxationFactor)



        # print ("Finished iteration {}. max delta gamma = {:.3f}".format(iteri, highestDeltaGamma))
        print(f"Finished iteration {iteri}. max delta gamma = {highestDeltaGamma} at idx = {highestIndex}")
        # Exit criterion based on change in delta gamma
        if abs(highestDeltaGamma) < 1e-1:
            bConverged = True
            break;

        iteri += 1

    print("Done :D")

    mainTurbine.extract_information_N_write()
    
    # plt.figure(1)
    # out = mainTurbine.extract_information_N_write()
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
