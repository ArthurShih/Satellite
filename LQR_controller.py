# TODO: Controlling a satellite in circular orbit
# LQR
import numpy as np 
import matplotlib.pyplot as plt
import control
from scipy.linalg import eig

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

# SYSTEM: m(rdd - r*thetad**2) = u1 - (GMm/r**2) + w1
#		  m(2*rd*thetad + r*thetadd) = u2 +w2
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
W = 0.1*np.eye(2)

#####################################################################################################
# Control(by u2 only)
print("CONTROL BY U2 ONLY")
Q = np.eye(4) # Q = I
rho = 1e6
R = rho # R = rho*I
Bu2 = np.array(Bu[:,1]).reshape((4,1)) # u2 only
[K2,X2,E2] = control.lqr(A,Bu2,Q,R) # returns K that means A-BK
print("K by using u2 only =\n",K2)
Ju2 = np.trace(X2.dot(Bw).dot(W).dot(np.transpose(Bw)))
print("\nLQR controller by using u2 only, cost J =\n"+str(Ju2))

# closed loop eigenvalue
Acl_u2 = A-Bu2.dot(K2)
eigenvalue_u2, eigenvector_u2 = eig(Acl_u2)
print("\nEigenvalue =\n",np.reshape(eigenvalue_u2, (4,1)))

# Simulation parameter
t = np.arange(0,10000,0.1)

# TODO: Solve system
from scipy.integrate import odeint

# closed loop B,C,D
Bcl = np.array([[0],
				[0],
				[0],
				[0]])
Ccl = np.eye(4)
Dcl = np.array([[0],
				[0],
				[0],
				[0]])

def system(x,t,Acl_u2):
	xdot = Acl_u2.dot(x)
	return xdot

x0 = np.array([1,0,0,1])

y_u2 = odeint(system, x0, t, args=(Acl_u2,))

satellite_radius = y_u2[:,0]
theta = y_u2[:,1]
rdot = y_u2[:,2]
thetadot = y_u2[:,3]

plt.figure(1)
plt.suptitle("State: by u2 only")
plt.subplot(221)
plt.title("rhat")
plt.plot(t, satellite_radius)

plt.subplot(222)
plt.title("theta")
plt.plot(t, theta)

plt.subplot(223)
plt.title("rhat dot")
plt.plot(t, rdot)
plt.xlabel("t")

plt.subplot(224)
plt.title("theta dot")
plt.plot(t, thetadot)
plt.xlabel("t")
print("------------------------------------------------------------------")

#####################################################################################################
# Control by u1 only
print("CONTROL BY U1 ONLY")
Bu1 = np.array(Bu[:,0]).reshape((4,1))
[K1,X1,E1] = control.lqr(A,Bu1,Q,R)
print("K by using u1 only =\n",K1)

Acl_u1 = A-Bu1.dot(K1)
Ju1 = np.trace(X1.dot(Bw).dot(W).dot(np.transpose(Bw)))
print("\nLQR controller by useing u1 only, cost J =\n"+str(Ju1))

eigenvalue_u1, eigenvector_u1 = eig(Acl_u1)
print("\nEigenvalue =\n",np.reshape(eigenvalue_u1, (4,1)))

def system(x,t,Acl_u1):
	xdot = Acl_u1.dot(x)
	return xdot

x0 = np.array([1,1,1,1])

y_u1 = odeint(system, x0, t, args=(Acl_u1,))

satellite_radius = y_u1[:,0]
theta = y_u1[:,1]
rdot = y_u1[:,2]
thetadot = y_u1[:,3]

plt.figure(2)
plt.suptitle("State: by u1 only")
plt.subplot(221)
plt.title("rhat")
plt.plot(t, satellite_radius)

plt.subplot(222)
plt.title("theta")
plt.plot(t, theta)

plt.subplot(223)
plt.title("rhat dot")
plt.plot(t, rdot)
plt.xlabel("t")

plt.subplot(224)
plt.title("theta dot")
plt.plot(t, thetadot)
plt.xlabel("t")
print("------------------------------------------------------------------")

# #####################################################################################################
# Control by u1 and u2
print("CONTROL BY U1 AND U2")
R = rho*np.eye(2)
[K,X,E] = control.lqr(A,Bu,Q,R)
print("K by using u1 and u2 =\n",K)

Acl = A-Bu.dot(K)
J = np.trace(X.dot(Bw).dot(W).dot(np.transpose(Bw)))
print("\nLQR controller by useing u1 and u2, cost J =\n"+str(J))

eigenvalue, eigenvector = eig(Acl)
print("\nEigenvalue =\n",np.reshape(eigenvalue, (4,1)))

def system(x,t,Acl):
	xdot = Acl.dot(x)
	return xdot

x0 = np.array([1,0,0,1])

y = odeint(system, x0, t, args=(Acl,))

satellite_radius = y[:,0]
theta = y[:,1]
rdot = y[:,2]
thetadot = y[:,3]

plt.figure(3)
plt.suptitle("State: by u1 and u2")
plt.subplot(221)
plt.title("rhat")
plt.plot(t, satellite_radius)

plt.subplot(222)
plt.title("theta")
plt.plot(t, theta)

plt.subplot(223)
plt.title("rhat dot")
plt.plot(t, rdot)
plt.xlabel("t")

plt.subplot(224)
plt.title("theta dot")
plt.plot(t, thetadot)
plt.xlabel("t")
plt.show()