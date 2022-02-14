# TODO: Controlling a satellite in circular orbit
# LQR
import numpy as np 
import matplotlib.pyplot as plt
import control
import scipy.signal as sig 
from scipy.linalg import eig
from scipy.integrate import odeint

# Model parameter
m = 100 # 100kg
W1 = 100**2 # covariance of w1
W2 = 10**2 # covariance of w2
W = np.array([[W1,0],
			  [0,W2]]) 

# SYSTEM:
# m*xdoubledot = w1
# y1 = x + w2
A = np.array([[0,1],
			  [0,0]])
Bw = np.array([[0,0],
			   [1/m,0]])
Cy = np.array([1,0])
Dyw = np.array([0,1])

# Find F
# by LQR function
At = np.transpose(A)
C = np.transpose(Cy).reshape(2,1)
Q = Bw.dot(W).dot(np.transpose(Bw))
R = Dyw.dot(W).dot(np.transpose(Dyw))
[F,X,E] = control.lqr(At,C,Q,R)
F = -np.transpose(F)
print("F =",F)
print("size of F is",np.shape(F))

# Simulation parameter
t_start = 0 
t_end = 100 
t_step = 0.1
t = np.arange(t_start, t_end, t_step)
N = int(t_end/t_step)
w = np.array([np.random.randn(1)*np.sqrt(W1), np.random.randn(1)*np.sqrt(W2)])

# System
f0 = np.array([10,0,0,0])
def model(f,t):
	x = f[0:2]
	xhat = f[2:4]
	[x1dot,x2dot] = A.dot(x)
	y = Cy.dot(x)
	F_new = np.reshape(F.dot(Cy.dot(xhat)-y), (2,)) #!!!!!!reshape to (2,)!!!!!
	[x1hatdot,x2hatdot] = A.dot(xhat) + F_new	
	fdot = np.array([x1dot,x2dot,x1hatdot,x2hatdot])
	return fdot

# Simulation
f = odeint(model, f0, t)

x = f[:,0:2]
xhat = f[:,2:4]
e = xhat - x
plt.figure()
plt.subplot(211)
plt.plot(t,x)
plt.plot(t,xhat)
plt.legend(["x1","x2","x1hat","x2hat"])
plt.subplot(212)
plt.plot(t,e)
plt.legend(["e1","e2"])
plt.show()