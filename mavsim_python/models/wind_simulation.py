"""
Class to determine wind velocity at any given moment,
calculates a steady wind speed and uses a stochastic
process to represent wind gusts. (Follows section 4.4 in uav book)
"""
from tools.transfer_function import TransferFunction
import numpy as np


class WindSimulation:
    def __init__(self, Ts, gust_flag = True, steady_state = np.array([[0., 0., 0.]]).T):
        # steady state wind defined in the inertial frame
        self._steady_state = steady_state

        #   Dryden gust model parameters (pg 56 UAV book)
        L_u = 200
        L_v = L_u
        L_w = 50
        V_a = 25
        gust_flag = True
        if gust_flag == True:
            sigma_u = 1.06
            sigma_v = sigma_u
            sigma_w = 0.7
        
        # Dryden transfer functions (section 4.4 UAV book) - Fill in proper num and den
        self.u_w = TransferFunction(num=np.array([[sigma_u*np.sqrt(2*V_a/(np.pi*L_u))]]), den=np.array([[1,V_a/L_u]]),Ts=Ts)
        self.v_w = TransferFunction(num=np.array([[sigma_v*np.sqrt(3*V_a/(np.pi*L_v)), sigma_v*np.sqrt(3*V_a/(np.pi*L_v))*V_a/(np.sqrt(3)*L_v)]]), den=np.array([[1, 2*V_a/L_v, (V_a/L_v)**2]]),Ts=Ts)
        self.w_w = TransferFunction(num=np.array([[sigma_w*np.sqrt(3*V_a/(np.pi*L_w)), sigma_w*np.sqrt(3*V_a/(np.pi*L_w))*V_a/(np.sqrt(3)*L_w)]]), den=np.array([[1, 2*V_a/L_w, (V_a/L_w)**2]]),Ts=Ts)
        self._Ts = Ts

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        gust = np.array([[self.u_w.update(np.random.randn())],
                         [self.v_w.update(np.random.randn())],
                         [self.w_w.update(np.random.randn())]])
<<<<<<< HEAD
        return np.concatenate(( self._steady_state, gust ))
=======
        gust = np.array([[0.],[0.],[0.]])
        return np.concatenate(( self._steady_state, gust ))

>>>>>>> a4cbd1edfb47c55225c7e222198175d833e75f75
