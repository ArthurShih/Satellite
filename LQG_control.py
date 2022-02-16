# TODO: Controlling a satellite in circular orbit
# LQG
import numpy as np 
import matplotlib.pyplot as plt
import control
from scipy.linalg import eig
from scipy.integrate import odeint

# Model Parameters
m = 100 # Satellite's mass
Rad = 6.37*1e3 # Radius of Earth
r = (Rad + 300)*1e3 # distance between earth and satellite
G = 6.673*1e-11 # universal gravitational constant
M = 5.98*1e24 # Earth's mass
k = G*M
omega = np.sqrt(k/r**3)
delta1 = 0.1 # invariance of noise w1
delta2 = 0.1 # invariance of noise w2
delta3 = 0.1 # invariance of noise v

# SYSTEM: m(rdd - r*thetad**2) = u1 - (GMm/r**2) + w1
#		  m(2*rd*thetad + r*thetadd) = u2 +w2
# noise measurement of theta(x2)
# Linearized system
A = np.array([[0,0,1,0],
	 		  [0,0,0,1],
	 		  [3*omega**2,0,0,2*r*omega],
	 		  [0,0,-2*omega/r,0]])
Bu = np.array([[0,0],
	  		   [0,0],
	  		   [1/m,0],
	  		   [0,1/(m*r)]])
Bw = np.array([[0,0],
	  		   [0,0],
	  		   [1/m,0],
	  		   [0,1/(m*r)]])
Cy = np.array([[0,1,0,0]])
Dyw = np.array([[0,0,1]])
W = np.array([[delta1,0,0],
			  [0,delta2,0],
			  [0,0,delta3]])
Q = np.eye(4)
R = np.eye(2)
###############################################################
# Augment matrices
Bwa = np.array([[0,0,0],
				[0,0,0],
				[1/m,0,0],
				[0,1/(m*r),0]])
Cz = np.array([[1,0,0,0],
			   [0,1,0,0],
			   [0,0,1,0],
			   [0,0,0,1],
			   [0,0,0,0],
			   [0,0,0,0]])
Dzu = np.array([[0,0],
				[0,0],
				[0,0],
				[0,0],
				[1,0],
				[0,1]]) # NOTICE: Cz'*Dzu = 0

# Find K
Qk = np.transpose(Cz).dot(Cz)
Rk = np.transpose(Dzu).dot(Dzu)
[K,X,E] = control.lqr(A,Bu,Q,R)
K = -K
print("K =\n",K)
eigvalue, eigvector = eig(A+Bu.dot(K))
print("\nEigenvalue of A+Bu*K =\n",eigvalue)


# Find F
Af = np.transpose(A)
Cf = np.transpose(Cy)
Qf = Bwa.dot(W).dot(np.transpose(Bwa))
Rf = Dyw.dot(W).dot(np.transpose(Dyw))
[F,Y,E] = control.lqr(Af,Cf,Qf,Rf)
F = -np.transpose(F)
print("\nF =\n",F)
eigvalue, eigvector = eig(A+F.dot(Cy))
print("\nEigenvalue of A+F*Cy =\n",eigvalue)
J1 = np.trace(Cz.dot(Y).dot(np.transpose(Cz)))
print("\nJ1 =\n",J1)
Bbar = X.dot(np.transpose(Cy)).dot((Dyw.dot(W).dot(np.transpose(Dyw)))**-0.5)
J2 = np.trace(np.transpose(Bbar).dot(X).dot(Bbar))
print("\nJ2 =\n",J2)
print("\nCost =\n",J1+J2)

# Zeros and poles
Acl = A+Bu.dot(K)+F.dot(Cy)
Bcl = -F
Ccl1 = K[0,:]
Ccl2 = K[1,:]
Dcl = 0
eigvalue, eigvector = eig(Acl)
print("\nEigenvalue of A+Bu*K+F*Cy =\n",eigvalue)

system1 = control.ss(Acl,Bcl,Ccl1,Dcl)
system2 = control.ss(Acl,Bcl,Ccl2,Dcl)
print("\n#1\n")
plt.figure(1)
pole, zero = control.pzmap(system1,plot=True)
print("\nZeros =\n",zero)
print("\nPoles =\n",pole)
print("\n#2\n")
plt.figure(2)
pole, zero = control.pzmap(system2,plot=True)
print("\nZeros =\n",zero)
print("\nPoles =\n",pole)

# Simulation parameters
t_start = 0 
t_end = 100000
t_step = 0.1
t = np.arange(t_start, t_end, t_step)
x0 = [1,0,0,1]
xhat0 = [0,0,0,0]
f0 = np.array(x0+xhat0)

# SYSTEM
def sys(f,t):
	x = f[0:4]
	xhat = f[4:8]
	xdot = A.dot(x) + Bu.dot(K).dot(xhat)
	xhatdot = -F.dot(Cy).dot(x) + (A+Bu.dot(K)+F.dot(Cy)).dot(xhat)
	return np.concatenate((xdot,xhatdot))

f = odeint(sys, f0, t)

x = f[:,0:4]
xhat = f[:,4:8]
e = x - xhat
e1 = e[:,0]
e2 = e[:,1]
e3 = e[:,2]
e4 = e[:,3]
plt.figure(3)
plt.subplot(221)
plt.plot(t,e1,"r-")
plt.title("r error")

plt.subplot(222)
plt.plot(t,e2,"r-")
plt.title("theta error")

plt.subplot(223)
plt.plot(t,e3,"r-")
plt.title("rdot error")
plt.xlabel("t")

plt.subplot(224)
plt.plot(t,e4,"r-")
plt.title("theta dot error")
plt.xlabel("t")
plt.show()
