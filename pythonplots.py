from numpy import *
import matplotlib.pyplot as plt


N10 = loadtxt("numbers-10.txt", skiprows = 3)
N100 = loadtxt("numbers2-100.txt", skiprows = 3)
N1000 = loadtxt("numbers3-1000.txt", skiprows = 3)

N1 = 12
N2 = 102
N3 = 1002

x1 = linspace(0, 1, N1)
x2 = linspace (0, 1, N2)
x3 = linspace(0, 1, N3)

u1 = N10[:,  1]
u2 = N100[:,  1]
u3 = N1000[:,  1]

def anal_u(x):
    return  1 - (1 - exp(-10))*x - exp(-10*x)

plt.figure(1)
plt.subplot(3, 1, 1)
plt.plot(x1, u1)
plt.plot(x1, anal_u(x1), 'k--')
plt.legend(["$N = 10$", "$Analytic$ $u(x)$"], fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$u(x)$', fontsize = 20)
plt.title('$Analytic$ $and$ $numerical$ $solution$ $of$ $u(x)$ $for$ $N = 10$', fontsize = 20)

plt.subplot(3, 1, 2)
plt.plot(x2, u2)
plt.plot(x2, anal_u(x2), 'k.-')
plt.legend(["$N = 100$","$Analytic$ $u(x)$"], fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$u(x)$', fontsize = 20)
plt.title('$Analytic$ $and$ $numerical$ $solution$ $of$ $u(x)$ $for$ $N = 100$', fontsize = 20)


plt.subplot(3, 1, 3)
plt.plot(x3, u3)
plt.plot(x3, anal_u(x3))
plt.legend(["$N = 1000$", "$Analytic$ $u(x)$"], fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$u(x)$', fontsize = 20)
plt.title('$Analytic$ $and$ $numerical$ $solution$ $of$ $u(x)$ $for$ $N = 1000$', fontsize = 20)

plt.figure(1).set_tight_layout(True)
plt.show()