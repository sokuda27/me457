#%%
import numpy as np

# Part a
psi = np.radians(2)
theta = np.radians(-10)
phi = np.radians(20)

A = np.array([[np.cos(psi), -np.sin(psi), 0], 
              [np.sin(psi), np.cos(psi), 0], 
              [0, 0, 1]])

B = np.array([[np.cos(theta), 0, np.sin(theta)], 
              [0, 1, 0], 
              [-np.sin(theta), 0, np.cos(theta)]])

C = np.array([[1, 0, 0], 
              [0, np.cos(phi), -np.sin(phi)], 
              [0, np.sin(phi), np.cos(phi)]])

R = np.dot(np.dot(A, B), C)
r = np.dot(R, [0.2, 0, 0]) + [0,0,-10]
print('Part a answer:')
print(r)

# Part b
V = np.array([15, 1, 0.5])
Vg = np.dot(R, V)
print('Part b answer:')
print(Vg)

# Part c
gamma = np.arctan(Vg[2]/(np.sqrt(Vg[0]**2 + Vg[1]**2)))
gamma = np.degrees(gamma)
print('Part c answer:')
print(gamma)

# part d
alpha = np.arctan(Vg[2]/Vg[0])
alpha = np.degrees(alpha)
print('Part d answer:')
print(alpha)

# part e
X = np.arctan(Vg[1]/Vg[0])
X = np.degrees(X)
print('Part e answer:')
print(X)