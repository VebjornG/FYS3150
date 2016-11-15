from numpy import *
from matplotlib.pyplot import *

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
File with T = 1.0
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
File_1 = loadtxt("Probability_E_T_1.txt", skiprows = 2)
T_1 = 1.0
E_variance_1 = 5.8552633E-05
E_1 = File_1[:]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
File with T = 2.4
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
File_2 = loadtxt("Probability_E_T_24.txt", skiprows = 2)
T_2 = 2.4
E_variance_2 = 0.020417803
E_2 = File_2[:]
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

#fig, ax = figure(1)
#rcParams.update({'font.size': 18})

subplot(121)
hist(E_2, 150, normed=1, facecolor='blue', alpha=0.75, label='$\sigma_E^2 = %f$, $T = 2.4$' % E_variance_2)
xlabel('Energy per spin, $E/n_{\mathrm{spins}}$', fontsize = 18)
ylabel('Probability, $P(E)$', fontsize = 18)
legend(loc='upper right', fancybox='True', fontsize = 20)
grid(True)

subplot(122)
hist(E_1, 5, normed=1, facecolor='blue', alpha=0.75, label='$\sigma_E^2 = %f$, $T = 1$' % E_variance_1)
xlabel('Energy per spin, $E/n_{\mathrm{spins}}$', fontsize = 18)
ylabel('Probability, $P(E)$', fontsize = 18)
legend(loc='upper right', fancybox='True', fontsize = 20)
grid(True)
#title('Energy histogram after reached equilibrium, T = %.1f kT/J' % T, fontsize = 15)

#savefig('Project4_d_probability_lowT.eps', format='eps', dpi=1000)
#savefig('Project4_d_probability_highT.eps', format='eps', dpi=1000)
savefig('Project4_d_probability.eps', format='eps', dpi=1000)
show()