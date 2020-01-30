
from crystalpy.util.Vector import Vector
import numpy


who = "Sean"
# Sean
if who == "Sean":
    S  = Vector.initializeFromComponents([0,0,0])
    M1 = Vector.initializeFromComponents([0,1.3,0])
    M2 = Vector.initializeFromComponents([0,1.3,0.28])
    M3 = Vector.initializeFromComponents([-0.293348	,1.3,0.78])
    M4 = Vector.initializeFromComponents([0.896652	,1.3,0.78])
    M5 = Vector.initializeFromComponents([0.896652	,0.29,0.78])
    DW = Vector.initializeFromComponents([2.45789,-0.302286	,0.78])
else:
    # Manuel
    S  = Vector.initializeFromComponents([0,0,0])
    M1 = Vector.initializeFromComponents([0,1.3,0])
    M2 = Vector.initializeFromComponents([0,1.3,0.28])
    M3 = Vector.initializeFromComponents([-0.29215859372,1.3,0.77797226692])
    M4 = Vector.initializeFromComponents([0.83781386898,1.3,0.78586109105])
    M5 = Vector.initializeFromComponents([0.83781386898,0.275,0.78586109105])
    DW = Vector.initializeFromComponents([ 2.45486318225, -0.30081939205, 0.79715040834])



S_M1  = M1.duplicate()
M1_M2 = M2.duplicate()
M2_M3 = M3.duplicate()
M3_M4 = M4.duplicate()
M4_M5 = M5.duplicate()
M5_DW = DW.duplicate()



S_M1  = M1.subtractVector(S)
M1_M2 = M2.subtractVector(M1)
M2_M3 = M3.subtractVector(M2)
M3_M4 = M4.subtractVector(M3)
M4_M5 = M5.subtractVector(M4)
M5_DW = DW.subtractVector(M5)



print("Angle M1: ",0.5*180/numpy.pi*S_M1.angle(M1_M2))
print("Angle M2: ",0.5*180/numpy.pi*M1_M2.angle(M2_M3))
print("Angle M3: ",0.5*180/numpy.pi*M2_M3.angle(M3_M4))
print("Angle M4: ",0.5*180/numpy.pi*M3_M4.angle(M4_M5))
print("Angle M5: ",0.5*180/numpy.pi*M4_M5.angle(M5_DW))

