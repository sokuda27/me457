class Integrator:
    """Integrator for a system of first-order ordinary differential equations
    of the form \dot x = f(t, x, u).
    """
    def __init__(self, dt, f):
        self.dt = dt
        self.f = f

    def step(self, t, x, u):
        raise NotImplementedError

class Euler(Integrator):
    def step(self, t, x, u):
        return x + self.dt * self.f(t, x, u)

class Heun(Integrator):
    def step(self, t, x, u):
        intg = Euler(self.dt, self.f)
        xe = intg.step(t, x, u) # Euler predictor step (k1)
        return x + 0.5*self.dt * (self.f(t, x, u) + self.f(t+self.dt, xe, u))
    
class RungeKutta4(Integrator):
    def step(self, t, x, u):
        k1 = self.f(t, x, u) # first predictor step
        k2 = self.f(t + self.dt/2, x + self.dt*k1/2, u)
        k3 = self.f(t + self.dt/2, x + self.dt*k2/2, u)
        k4 = self.f(t + self.dt, x + self.dt*k3, u)
        return x + self.dt*(k1/6 + k2/3 + k3/3 + k4/6)




        

