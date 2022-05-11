from Lifting_Line import *
# Ugly but it worksss

mainTurbine = Turbine(0, (10,10,10))

# This file contains the main iteration over two turbines.
maxiterations = 10
iteri = 1
mainTurbine.set_circulations_horseshoes()
while iteri < maxiterations:
    mainTurbine.SetInducedVelocityForHorseshoes()
    mainTurbine.set_circulations_horseshoes()

    iteri += 1