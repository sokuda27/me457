import numpy as np
import matplotlib.pyplot as plt
import integrators as intg


# Nonlinear state space form:
#  xdot = f(t, x, u), 
#   t: time
#   x: state vector
#   u: input
def f(t, x, u):
    return np.array([x[1], -x[0]-(0.25)*x[1]])
        

t = 0; x = np.array([0, 1]); u = 0
dt = 0.1; n = 100

# integrator = intg.Euler(dt, f)
integrator = intg.Heun(dt, f)
integrator1 = intg.RungeKutta4(dt, f)

t_history = [0]
x_history = [x]
x_history1 = [x]
for i in range(n):
    
    x = integrator.step(t, x, u)
    x1 = integrator1.step(t, x, u)
    t = (i+1) * dt

    t_history.append(t)
    x_history.append(x)
    x_history1.append(x1)

intg.__doc__
plt.figure()
plt.plot(t_history, x_history)
plt.legend(["heun x", "heun x_dot"])

plt.figure()
plt.plot(t_history, x_history1)
plt.legend(["rk4 x", "rk4 x_dot"])
plt.show()