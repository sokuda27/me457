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
        xe = intg.step(t, x, u) # Euler predictor step
        return x + 0.5*self.dt * (self.f(t, x, u) + self.f(t+self.dt, xe, u))
    
class RungeKutta3(Integrator):
    def step(self, t, x, u):
        intg = Heun(self.dt, self.f)
        x2 = intg.step(t, x, u) # Heun step (rk2)
        return x + 0.5*self.dt * (self.f(t, x, u) + self.f(t+self.dt, x2, u))
    
class RungeKutta4(Integrator):
    def step(self, t, x, u):
        intg = RungeKutta3(self.dt, self.f)
        x3 = intg.step(t, x, u) # RK3 step
        return x + self.dt * (self.f(t, x, u) + self.f(t+self.dt, x3, u))