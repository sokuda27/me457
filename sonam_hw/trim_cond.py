"""
compute_trim 
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        12/29/2018 - RWB
"""
import numpy as np
from scipy.optimize import minimize
from rotations import Euler2Quaternion
from msg_delta import MsgDelta
import q4_parameters as param

# def objective(x):
#     return(x[0]**2 + x[1]**2)

# x0 = np.array([[5], [2]])

# cons = ({'type': 'eq',
#          'fun': lambda x: np.array([
#              x[0] + x[1] - 2
#             ])
#          })

# res = minimize(objective, x0, method='SLSQP', constraints=cons)
# print("xstar =", res.x)

# inital conditions x
p_n0 = 1; p_e0 = 0; p_d0 = 2
u_0 = 1; v_0 = 1; w_0 = 0
phi_0 = 0; theta_0 = 0; psi_0 = 0
p_0 = 0; q_0 = 0; r_0 = 0

x0 = np.array([[p_n0], [p_e0], [p_d0], [u_0], [v_0], [w_0], [phi_0], [theta_0], [psi_0], [p_0], [q_0], [r_0]])

def compute_trim(x, Va, gamma):
    # define initial state and input
    e0 = Euler2Quaternion(0., gamma, 0.)
    # state0 = x
    # delta0 = MsgDelta(elevator=x[12])
    # x0 = np.concatenate((state0, delta0.to_array()), axis=0)
    # define equality constraints
    cons = ({'type': 'eq',
             'fun': lambda x: np.array([
                                x[3]**2 + x[4]**2 + x[5]**2 - Va**2,  # magnitude of velocity vector is Va
                                x[4],  # v=0, force side velocity to be zero
                                x[6]**2 + x[7]**2 + x[8]**2 - 1.,  # force quaternion to be unit length
                                x[6],  # e1=0  - forcing e1=e3=0 ensures zero roll and zero yaw in trim
                                x[8],  # e3=0
                                x[9],  # p=0  - angular rates should all be zero
                                x[10],  # q=0
                                x[11],  # r=0
                                ]),
             'jac': lambda x: np.array([
                                [0., 0., 0., 2*x[3], 2*x[4], 2*x[5], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 2*x[6], 2*x[7], 2*x[8], 2*x[9], 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
                                [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
                                ])
             }) 
    # solve the minimization problem to find the trim states and inputs

    res = minimize(trim_objective_fun, x, method='SLSQP', args=(Va, gamma),
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = MsgDelta(elevator=res.x.item(13),
                          aileron=res.x.item(14),
                          rudder=res.x.item(15),
                          throttle=res.x.item(16))
    trim_input.print()
    print('trim_state=', trim_state.T)
    return trim_state, trim_input


def trim_objective_fun(x, Va, gamma):
    # objective function to be minimized
    desired_trim_state_dot = np.array([[0., 0., -Va*np.sin(gamma), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]]).T
    curr_dot = param.f(0, x, 0)
    tmp = desired_trim_state_dot - curr_dot
    J = np.linalg.norm(tmp[2:13])**2.0
    print(J)
    return J

compute_trim(x0, 25, 0)
