from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sherif'})


data = loadtxt("results.txt")
datax = loadtxt("results_x.txt")
data1 = loadtxt("results1.txt")
data2 = loadtxt("results2.txt")
data3 = loadtxt("results3.txt")
data4 = loadtxt("results4.txt")
data5 = loadtxt("results5.txt")

"""
######### load data ################
x = datax
phi_1 = data[:,0]
phi_2 = data[:,1]
phi_3 = data[:,2]
phi_4 = data[:,3]
phi_5 = data[:,4]
phi_6 = data[:,5]


phi_3[-2] = phi_3[-3]
phi_2[-2] = phi_2[-3]

phi_3[-1] = phi_3[-2]
phi_2[-1] = phi_2[-2]
phi_1[-1] = phi_1[-2]

psi_1 = data1[:,0]
psi_2 = data1[:,1]
psi_3 = data1[:,2]
psi_4 = data1[:,3]
psi_5 = data1[:,4]
psi_6 = data1[:,5]

for i in range(590, 600):
	psi_1[i] = psi_1[i-1]
	psi_2[i] = psi_2[i-1]
	psi_3[i] = psi_3[i-1]
	psi_4[i] = psi_4[i-1]
	psi_5[i] = psi_5[i-1]
	psi_6[i] = psi_6[i-1]

psi_3[-2] = psi_3[-3]
psi_2[-2] = psi_2[-3]

psi_3[-1] = psi_3[-2]
psi_2[-1] = psi_2[-2]
psi_1[-1] = psi_1[-2]

v_1 = data2[:,0]
v_2 = data2[:,1]
v_3 = data2[:,2]
v_4 = data2[:,3]
v_5 = data2[:,4]
v_6 = data2[:,5]

for i in range(590, 600):
	v_1[i] = v_1[i-1]
	v_2[i] = v_2[i-1]
	v_3[i] = v_3[i-1]
	v_4[i] = v_4[i-1]
	v_5[i] = v_5[i-1]
	v_6[i] = v_6[i-1]
v_3[-1] = v_3[-2]
v_2[-1] = v_2[-2]

delta_1 = data3[:,0]
delta_2 = data3[:,1]
delta_3 = data3[:,2]
delta_4 = data3[:,3]
delta_5 = data3[:,4]
delta_6 = data3[:,5]


deltab_1 = data4[:,0]
deltab_2 = data4[:,1]
deltab_3 = data4[:,2]
deltab_4 = data4[:,3]
deltab_5 = data4[:,4]
deltab_6 = data4[:,5]


theta_1 = data5[:,0]
theta_2 = data5[:,1]
theta_3 = data5[:,2]
theta_4 = data5[:,3]
theta_5 = data5[:,4]
theta_6 = data5[:,5]

for i in range(590, 600):
	theta_1[i] = theta_1[i-1]
	theta_2[i] = theta_2[i-1]
	theta_3[i] = theta_3[i-1]
	theta_4[i] = theta_4[i-1]
	theta_5[i] = theta_5[i-1]
	theta_6[i] = theta_6[i-1]


theta_3[-1] = theta_3[-2]
theta_2[-1] = theta_2[-2]


########## Phi #####################
plt.plot(x, phi_1)
plt.plot(x, phi_2)
plt.plot(x, phi_3)
plt.plot(x, phi_4)
plt.plot(x, phi_5)
plt.plot(x, phi_6)
plt.ylabel('$d\phi/dx$')
plt.xlabel('x = ln a ')
plt.xlim([-18, -0])
plt.title('The derivative of the gravitational potential $\phi$')
#plt.yscale('log')
#plt.xscale('log')
"""
######################## milestone 4############
S_data = loadtxt("results_S.txt")
theta_data = loadtxt("results_theta.txt")
theta_data2 = loadtxt("results_theta2.txt")
theta_data3 = loadtxt("results_theta3.txt")

phi_data = loadtxt("results_phi.txt")
cl_data = loadtxt("results_cl.txt")

plank_lo = loadtxt("plankhi.txt", unpack=True, skiprows=3)
plank_hi = loadtxt("planklo.txt", unpack=True, skiprows=3)

j_data = loadtxt("results_j.txt")

lhi = plank_hi[0]
llo = plank_lo[0]

plank_l = hstack((llo,lhi))
cl_plank = hstack([plank_lo[1], plank_hi[1]])

#print shape(plank_l), shape(cl_plank)
H_0 = 0.7 * 100. * 1000. / 3.08568025e22
c = 2.99792458e8

x = phi_data[:,0]
S = phi_data[:,1]
j_l = phi_data[:,2]

x_low = S_data[:,0]
S_lores = S_data[:,1]

k = theta_data[:,0]
theta = theta_data[:,1]
#j_l = theta_data[:,2]

theta_2 = theta_data2[:,0]
theta_50 = theta_data2[:,1]
theta_200 = theta_data2[:,2]

theta_500 = theta_data3[:,0]
theta_800 = theta_data3[:,1]
theta_1200 = theta_data3[:,2]

l = cl_data[:,0]
l_int = [ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200]
cl = cl_data[:,1]
#cl2 = cl_data[:,2]

j1 = j_data[:,0]
j2 = j_data[:,1]
j3 = j_data[:,2]

#print k[4000]
	
#print theta
#print phi2[0], phi2[1], phi2[2], phi2[-1]
#phi2[0] = 1

#print H_0

#for i in range(550,600 ):
	#S_lores[i] = S_lores[i-1]
	#S[i] = S[i-1]
	#theta[i] = theta[i-1]
	#phi2[i] = phi2[i-1]
	#theta[i] = theta[i-1]


#plt.yscale('log')
#plt.plot(l, cl/max(cl)*5775.)
#plt.xlabel('l')
#plt.ylabel('$l(l+1)C_l / 2\pi$')

#plt.plot(plank_l, cl_plank)

#plt.plot(c*k/H_0,theta**2/(c*k)/(10**(-6)*H_0**(-1)))
#plt.xlim([8, 0])
#plt.plot(c*k/H_0, theta**2/(c*k)*H_0)
#plt.plot(x,S*j_l/10.**3)
#plt.xlim([-10,0])
#plt.ylabel('$S(k,x)j_l(\eta_0 - \eta(x)) / 10^3$')

#plt.plot(x, j1/1e3)
plt.plot(x, j2/1e3)
plt.plot(x, j3/1e3)
plt.xlabel('x')
plt.ylabel('$j_l/ 10^3$')
plt.title('Spherical Bessel functions')

#plt.plot(c*k/H_0, theta_2/(c*k)*H_0, label='$l = 2$') #2.*(2.+1.)*
#plt.plot(c*k/H_0, theta_50/(c*k)*H_0, label='$l = 50$') #50.*(50.+1.)*
#plt.plot(c*k/H_0, theta_200/(c*k)*H_0, label='$l = 200$') #200.*(200.+1.)*
#plt.plot(c*k/H_0, theta_500/(c*k)*H_0, label='$l = 500$')
#plt.plot(c*k/H_0, theta_800/(c*k)*H_0, label='$l = 800$')
#plt.plot(c*k/H_0, theta_1200/(c*k)*H_0, label='$l = 1200$')
#plt.legend(loc= 'best')
#plt.xlabel('$ck/H_0$')
#plt.ylabel('$\Theta_l H_0/ck$')


plt.show()
