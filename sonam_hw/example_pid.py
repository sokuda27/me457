import numpy as np
import parameters as P
from integrators_sonam import get_integrator
from pid import PIDControl


class Controller:
    def __init__(self):
        pass
        
    def update(self, r, y):
        pass
    
class System:
    def __init__(self):
        self._tau = 0.01
        self._k = 10
        self._y = 0        
    
    def state(y, u):
        y_dot = super().k/super().tau*u - y/super().tau
        return y_dot
   
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