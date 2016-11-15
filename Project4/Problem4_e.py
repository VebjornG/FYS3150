from numpy import *
from matplotlib.pyplot import *

L40 = loadtxt("Data_MPI_40.txt", skiprows = 1)
L60 = loadtxt("Data_MPI_60.txt", skiprows = 1)
L100 = loadtxt("Data_MPI_100.txt", skiprows = 1)
L140 = loadtxt("Data_MPI_140.txt", skiprows = 1)

Temp = L40[:, 0]

Heat_capacity_40 = L40[:, 2]
Heat_capacity_60 = L60[:, 2]
Heat_capacity_100 = L100[:, 2]
Heat_capacity_140 = L140[:, 2]

Susceptibility_40 = L40[:, 4]
Susceptibility_60 = L60[:, 4]
Susceptibility_100 = L100[:, 4]
Susceptibility_140 = L140[:, 4]

Energy_40 = L40[:, 1]
Energy_60 = L60[:, 1]
Energy_100 = L100[:, 1]
Energy_140 = L140[:, 1]

Magnetization_40 = L40[:, 5]
Magnetization_60 = L60[:, 5]
Magnetization_100 = L100[:, 5]
Magnetization_140 = L140[:, 5]


#print amax(Heat_capacity_140)
max_x = Temp[Heat_capacity_140.argmax()]
print "%.6f" %max_x

plot(Temp, Heat_capacity_140, '.-r')
plot(Temp, Heat_capacity_100, '.-b')
plot(Temp, Heat_capacity_60, '.-g')
plot(Temp, Heat_capacity_40, '.-y')
axvline(x=2.269, ymin=0.0, ymax = 1, linewidth=0.5, color='k')
xlabel(r'Temperature, $k_BT/J$', fontsize = 18)
ylabel(r'Heat capacity, $\langle C_v \rangle$', fontsize = 18)
legend([r"$L = 140$", r"$L = 100$", r"$L = 60$", r"$L = 40$", r"$T_C \approx 2.269$"], loc='upper left', fancybox='True', fontsize = 20)
savefig("Heat_capacity.png")
show()   
hold('on')


plot(Temp, Susceptibility_140, '.-r')
plot(Temp, Susceptibility_100, '.-b')
plot(Temp, Susceptibility_60, '.-g')
plot(Temp, Susceptibility_40, '.-y')
axvline(x=2.269, ymin=0.0, ymax = 1, linewidth=0.5, color='k')
xlabel(r'Temperature, $k_BT/J$', fontsize = 18)
ylabel(r'Susceptibility, $\langle \chi \rangle$', fontsize = 18)
legend([r"$L = 140$", r"$L = 100$", r"$L = 60$", r"$L = 40$", r"$T_C \approx 2.269$"], loc='upper left', fancybox='True', fontsize = 20)
savefig("Susceptibility.png")
show()
hold('on')


plot(Temp, Energy_40, '.-r')
plot(Temp, Energy_60, '.-b')
plot(Temp, Energy_100,'.-g')
plot(Temp, Energy_140, '.-y')
axvline(x=2.269, ymin=0.0, ymax = 1, linewidth=0.5, color='k')
xlabel(r'Temperature, $k_BT/J$', fontsize = 18)
ylabel(r'Energy, $\langle E \rangle$', fontsize = 18)
legend([r"$L = 40$", r"$L = 60$", r"$L = 100$", r"$L = 140$", r"$T_C \approx 2.269$"], loc='upper left', fancybox='True', fontsize = 20)
savefig("Energy.png")
show()
hold('on')

plot(Temp, Magnetization_140, '.-r')
plot(Temp, Magnetization_100, '.-b')
plot(Temp, Magnetization_60, '.-g')
plot(Temp, Magnetization_40, '.-y')
axvline(x=2.269, ymin=0.0, ymax = 1, linewidth=0.5, color='k')
xlabel(r'Temperature, $k_BT/J$', fontsize = 18)
ylabel(r'Magnetization, $\langle |M| \rangle$', fontsize = 18)
legend([r"$L = 140$", r"$L = 100$", r"$L = 60$", r"$L = 40$", r"$T_C \approx 2.269$"], loc='lower left', fancybox='True', fontsize = 20)
savefig("Magnetization.png")
show()