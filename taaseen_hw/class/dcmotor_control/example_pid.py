import numpy as np
import parameters as P
from integrators import get_integrator
from pid import PIDControl
from matplotlib import pyplot as plt


class Controller:
    def __init__(self):
        """Initialize the controller with parameters."""
        # Initialize PID parameters
        self.kp = P.kp
        self.ki = P.ki
        self.kd = P.kd
        self.umax = P.umax
        self.udead = P.udead
        self.sigma = P.sigma
        self.Ts = P.Ts
        self.integrator = get_integrator(P.Ts, None, integrator="Heun")
        self.pid = PIDControl(self.kp, self.ki, self.kd, self.umax, self.sigma, self.Ts)


    def update(self, r, y):
        u = self.pid.PID(r, y)
        return u
    
class System:
    def __init__(self):
        """Initialize the system with parameters."""
        # Initialize system parameters
        self.K = P.K
        self.tau = P.tau
        self.umax = P.umax
        self.integrator = get_integrator(P.Ts, None, integrator="Heun")
        self.state = np.zeros(1)
        pass        
   
    def update(self, u):
        state_dot = -(1/self.tau)*self.state + (self.K/self.tau)*u
        self.state += state_dot*P.Ts
        return self.state[0]

# Init system and feedback controller
system = System()
controller = Controller()

# Simulate step response
t_history = [0]
y_history = [0]
u_history = [0]

r = 1
y = 0
t = 0
for i in range(P.nsteps):
    u = controller.update(r, y) 
    y = system.update(u) 
    t += P.Ts

    t_history.append(t)
    y_history.append(y)
    u_history.append(u)

# Plot response y due to step change in r
plt.figure(figsize=(10, 5))
plt.subplot(2, 1, 1)
plt.plot(t_history, y_history, label="System Output (y)")
plt.title("Step Response of the System")
plt.xlabel("Time [s]")
plt.ylabel("Output (y)")
plt.grid(True)

# Plot actuation signal (u)
plt.subplot(2, 1, 2)
plt.plot(t_history, u_history, label="Control Signal (u)", color="orange")
plt.title("Control Signal (u) over Time")
plt.xlabel("Time [s]")
plt.ylabel("Control Signal (u) [V]")
plt.grid(True)

plt.tight_layout()
plt.show()