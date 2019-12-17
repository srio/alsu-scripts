
import numpy

def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm

def get_q(p_foc,R,incidence):
    mm = (2.0 / ( numpy.cos(incidence * numpy.pi / 180) * R) - (1.0 / p_foc))
    return 1/mm

if __name__ == "__main__":


    print("ALS M2 Radius: ",get_R(1.580, 5.87 - 1.58, 90 - 22.3))
    print("ALS M5 Radius: ", get_R(2.76, 5.87 - 4.25, 45))

    print("ALSU M2 Radius: ",get_R(1.580, 5.87 - 1.58, 74.933333))
    print("ALSU M4 Radius: ", get_R(2.76, 5.87 - 4.25, 45))

    print("q0",get_q(1.580,8.884408370766353,74.933333),5.87 - 1.58 )

    drift = 75.0 * 1e-3
    lengthBM = 500.0 * 1e-3
    lengthAB = 305 * 1e-3
    p0 = 1.58 - ( lengthBM/2 + drift*2 + lengthAB + lengthBM/2)
    print("p0 for Mag 8: ",p0)


    print("q0 for Mag8",get_q(p0,8.88440,74.933333))