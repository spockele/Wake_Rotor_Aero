from Lifting_Line import *
# Ugly but it worksss

mainTurbine = Turbine(0, (10,10,10))

# This file contains the main iteration over two turbines.
maxiterations = 100
iteri = 1
bConverged = False
mainTurbine.set_circulations_horseshoes()
while iteri < maxiterations:
    mainTurbine.SetInducedVelocityForHorseshoes()
    highestDeltaGamma = mainTurbine.set_circulations_horseshoes()



    print ("Finished iteration {}. max delta gamma = {:.3f}".format(iteri, highestDeltaGamma))
    # Exit criterion based on change in delta gamma
    if highestDeltaGamma < 0.1:
        bConverged = True
        break;

    iteri += 1