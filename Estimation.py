# TODO: Controlling a satellite in circular orbit
# LQR
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

# SYSTEM: m(rdd - r*thetad**2) = u1 - (GMm/r**2) + w1
#		  m(2*rd*thetad + r*thetadd) = u2 +w2
# Linearized system
A = np.array([[0,0,1,0],
	 		  [0,0,0,1],
	 		  [3*omega**2,0,0,2*r*omega],
	 		  [0,0,-2*omega/r,0]])
Bw = np.array([[0,0,0],
	  		   [0,0,0],
	  		   [1/m,0,0],
	  		   [0,1/(m*r),0]])
Cy = np.array([[0,1,0,0]])
Dyw = np.array([0,0,1])
W = 0.1*np.eye(3)
#################################################################
# Find F
At = np.transpose(A)
C = np.transpose(Cy)
Q = Bw.dot(W).dot(np.transpose(Bw))
R = Dyw.dot(W).dot(np.transpose(Dyw))
[F,X,E] = control.lqr(At,C,Q,R)
F = -np.transpose(F)
print(F)
eigenvalue, eigenvector = eig(A+F.dot(Cy))
print(eigenvalue)

# Simulation parameters
t_start = 0
t_end = 1500000
t_step = 0.1
t = np.arange(t_start, t_end, t_step)
x0 = np.array([0,0,0,0])
xhat0 = np.array([0,0,-6,0])
f0 = np.concatenate((x0,xhat0))
w = np.random.randn(1)*np.eye(3)

# System
def sys(f,t):
	x = f[0:4]
	xhat = f[4:8]
	xdot = A.dot(x)
	y = Cy.dot(x)
	xhatdot = A.dot(xhat) + np.reshape(F.dot(Cy.dot(xhat)-y),(4,))
	return np.concatenate((xdot,xhatdot))

f = odeint(sys, f0, t)

x = f[:,0:4]
xhat = f[:,4:8]
plt.figure()
plt.subplot(421)
plt.plot(t,x[:,0])
plt.plot(t,xhat[:,0])
plt.subplot(422)
plt.plot(t,x[:,1])
plt.plot(t,xhat[:,1])
plt.subplot(423)
plt.plot(t,x[:,2])
plt.plot(t,xhat[:,2])
plt.subplot(424)
plt.plot(t,x[:,3])
plt.plot(t,xhat[:,3])
e = xhat - x
plt.subplot(425)
plt.plot(t,e[:,0])
plt.subplot(426)
plt.plot(t,e[:,1])
plt.subplot(427)
plt.plot(t,e[:,2])
plt.subplot(428)
plt.plot(t,e[:,3])
plt.show()