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

class RK3(Integrator):
    def step(self, t, x, u):
        intE = Euler(self.dt, self.f)
        xE = intE.step(t, x, u)
        intH = Heun(self.dt, self.f)
        xH = intH.step(t, x, u)
        return x + (self.dt/6) * (self.f(t, x, u) + 4*self.f(t+self.dt/2, xE, u) + self.f(t+self.dt, xH, u))

class RK4(Integrator):
    def step(self, t, x, u):
        intE = Euler(self.dt, self.f)
        xE = intE.step(t, x, u)
        intH = Heun(self.dt, self.f)
        xH = intH.step(t, x, u)
        int3 = RK3(self.dt, self.f)
        x3 = int3.step(t, x, u)
        return x + (self.dt/6) * (self.f(t, x, u) + 2*self.f(t+self.dt/2, xE, u) + 2*self.f(t+self.dt/2, xH, u) + self.f(t+self.dt, x3, u))

# class RK4(Integrator):
#     def step(self, t, x, u):
#         k1 = self.f(t, x, u)
#         k2 = self.f(t + self.dt/2, x + 0.5 * self.dt * k1, u)
#         k3 = self.f(t + self.dt/2, x + 0.5 * self.dt * k2, u)
#         k4 = self.f(t + self.dt, x + self.dt * k3, u)
#         return x + (self.dt/6) * (k1 + 2*k2 + 2*k3 + k4)