import numpy as np

class SelfTuningPIDControl:
    def __init__(self, kp=0.5, ki=0.05, kd=0.1, Ts=0.01, limit=np.inf, learning_rate=0.01):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.Ts = Ts
        self.limit = limit
        self.learning_rate = learning_rate

        self.integrator = 0.0
        self.prev_error = 0.0
        self.differentiator = 0.0
        self.prev_measurement = 0.0

    def update(self, setpoint, measurement):
        error = setpoint - measurement

        # Integrator
        self.integrator += 0.5 * self.Ts * (error + self.prev_error)

        # Derivative (band-limited)
        self.differentiator = (2.0 * (measurement - self.prev_measurement)
                               + (2.0 - self.Ts) * self.differentiator) / (2.0 + self.Ts)

        # PID output
        u = self.kp * error + self.ki * self.integrator - self.kd * self.differentiator
        u = self._saturate(u)

        # Online gradient descent (naive)
        self.kp += self.learning_rate * error * error
        self.ki += self.learning_rate * error * self.integrator
        self.kd -= self.learning_rate * error * self.differentiator

        # Optional: clamp gains
        self.kp = np.clip(self.kp, 0.0, 10.0)
        self.ki = np.clip(self.ki, 0.0, 5.0)
        self.kd = np.clip(self.kd, 0.0, 5.0)

        self.prev_error = error
        self.prev_measurement = measurement
        return u

    def _saturate(self, u):
        return max(min(u, self.limit), -self.limit)
