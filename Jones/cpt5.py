import numpy as np
import scipy.constants as cte
import math
import matplotlib.pyplot as plt

def stress_strain_for_CL(e, n, T):

    stress = n*cte.k*T*((1+e) - 1/math.pow((1+e),2))

    return stress

def question5_2():
    strain = [
        0.000, 0.162, 0.270, 0.433, 0.678, 0.950,
        1.358, 1.657, 2.338, 2.964, 3.480, 4.350,
        4.973, 5.461, 6.190, 6.403, 6.699, 6.914,
        7.019, 7.151, 7.256, 7.361, 7.489
    ]

    stress = [  # Force / unstrained area (N mm^-2)
        0.000, 0.152, 0.246, 0.327, 0.420, 0.489,
        0.605, 0.697, 0.882, 1.067, 1.253, 1.613,
        1.986, 2.313, 3.050, 3.448, 3.811, 4.151,
        4.503, 4.878, 5.242, 5.605, 6.321
    ]

    for idx in range(len(stress)):
        stress[idx] = 1e6*stress[idx] #stress to N per m2

    T=293
    n_density = stress[1]/(3*cte.k*T*strain[1])
    print(n_density)


    strain_curve = np.linspace(0, 8, 100)
    pred_stress = []

    #print(stress_strain_for_CL(1.627, n_density, T))
    for value in strain_curve:
        pred_stress.append(stress_strain_for_CL(value, n_density, T))

    # print(strain_curve)
    # print(pred_stress)

    plt.plot(strain_curve, pred_stress, marker='', linestyle='-')
    plt.plot(strain, stress, marker='o', linestyle='')
    plt.show()


question5_2()
