from numpy import *

L = [40, 60, 100, 140]

T_C_CV = [2.295, 2.285, 2.280, 2.275]
T_C_X = [2.320, 2.300, 2.280, 2.280]

Theoretical = 2.269

#x1 = zeros(len(T_C_CV) - 1)
#x2 = zeros(len(T_C_X) - 1)

#for i in range(len(T_C_CV)-1):
#    x1[i] = (T_C_CV[i]*L[i] - T_C_CV[i+1]*L[i+1])/(L[i] - L[i+1])
#    x2[i] = (T_C_X[i]*L[i] - T_C_X[i+1]*L[i+1])/(L[i] - L[i+1])

#print "Numerical T_C = ", (mean(x1) + mean(x2))/2.

#print "Theoretical T_C = ", Theoretical

#print "Error in estimate = ", abs((mean(x1) + mean(x2))/2. - Theoretical)

T_limit_Cv = zeros(len(T_C_CV) - 1)
T_limit_X = zeros(len(T_C_X) - 1)

for i in range(len(T_C_CV) - 1):
    T_limit_Cv[i] = 1./L[i] - T_C_CV[i]
    T_limit_X[i] = 1./L[i] - T_C_X[i]

T_lim_Cv = abs(mean(T_limit_Cv))
T_lim_X = abs(mean(T_limit_X))
#print T_lim_Cv, T_lim_X
#print  abs(T_lim_Cv + T_lim_X)/2 - Theoretical

print "Numerical T_C = ", (T_lim_Cv+ T_lim_X)/2.

print "Theoretical T_C = ", Theoretical

print "Error in estimate = ", abs((T_lim_Cv + T_lim_X)/2. - Theoretical)