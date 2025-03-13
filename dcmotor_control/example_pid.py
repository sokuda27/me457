import numpy as np
import parameters as P
from integrators import get_integrator
from pid import PIDControl


class Controller:
    def __init__(self):
        """Initialize the controller with parameters."""
        # Initialize PID parameters
        self.kp = P.kp
        self.ki = P.ki
        self.kd = P.kd
        self.Ts = P.Ts
        # Initialize integrator
        self.integrator = get_integrator(P.Ts, None, integrator="Heun")
        # Initialize state
        self.state = np.zeros(2)  # [integral, previous error]
        pass


    def update(self, r, y):
        """Update the controller with new reference and measurement.
        Args:
            r (float): Reference value.
            y (float): Measurement value.
        Returns:
            u (float): Control signal.
        """
        # Compute error
        e = r - y
        # Compute control signal
        u = PIDControl(e, P.kp, P.ki, P.kd, P.Ts)
        return u
        pass
    
class System:
    def __init__(self):
        """Initialize the system with parameters."""
        # Initialize system parameters
        self.K = P.K
        self.tau = P.tau
        self.umax = P.umax
        # Initialize integrator
        self.integrator = get_integrator(P.Ts, None, integrator="Heun")
        # Initialize state
        self.state = np.zeros(1)
        pass        
    
   
    def update(self, u):
        pass

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


# Plot actuation signal