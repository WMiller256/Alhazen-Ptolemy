import numpy as np
from cmath import sqrt, sin, cos
from math import acos, atan

tolerance = 1e-9

def f_E7(obs, c, b):
# Convenience function which calculates the C7 coefficient
# for a given observer angle and radius and source radius
	s_obs = sin(obs)
	c_obs = cos(obs)
	c_obs2 = c_obs**2
	s_obs2 = s_obs**2
	c2 = c**2
	b2 = b**2

	E0 = c2 - 4 + b2 + 2 * b * c * s_obs
	E1 = 24 * b * c_obs2 * (b - c * s_obs) - 48 * (-1 + c2) * c_obs2 + E0**2 
	E2 = -432 * b2 * (-1 + c2) * c_obs2**2 + 432 * c_obs2 * (b - c * s_obs)**2 + 72 * b * c_obs2 * (b - c * s_obs) * E0 + 288 * (-1 + c2) * c_obs2 * E0 + 2 * E0**3 
	E3 = (E2 + sqrt(-4 * E1**3 + E2**2))**(1.0 / 3.0) 
	E4 = E1 / (6 * 2**(2.0 / 3.0) * E3) + E3 / (12 * 2**(1.0 / 3.0)) 
	E5 = sqrt((b2 * c_obs2) / 4.0 - E0 / 6.0 + E4)	

	return ((b2 * c_obs2) / 2 - E0 / 3.0 - E4 + (b**3 * c_obs**3 - 4 * c_obs * (b - c * s_obs) - b * c_obs * E0) / (4 * E5)).real

def onefinite(obs, c, branch=0):
# Calculates a specific solution for the one-finite case without 
# branch deduction. This is more useful for plotting all the 
# solutions over a given interval than for finding the physical
# specular point

	s_obs = sin(2 * obs)
	c_obs = cos(2 * obs)
	c2 = c**2 + 0j
	c2_4 = c2 - 4
	c2_43 = c2_4**3
	c4 = c**4

	D0 = 1152 * c2_4 * (-1 + c2 - c_obs + c2 * c_obs) + 864 * c2 * s_obs**2
	D1 = sqrt(-4 * (160 - 128 * c2 + 4 * c4 + 96 * c_obs - 96 * c2 * c_obs)**3 + (-16 * c2_43 - D0)**2)
	D2 = (16 * c2_43 + D0 - D1)**(1.0 / 3.0)
	D3 = (40 - 32 * c2 + c4 + 24 * c_obs - 24 * c2 * c_obs) / (3 * 2**(2.0 / 3.0) * D2) + D2 / (24 * 2**(1.0 / 3.0))
	D4 = sqrt(-0.25 * c2_4 + (1.0 / 12.0) * c2_4 + D3)

	p = -0.5 * D4 - 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 - (c * s_obs) / (2 * D4))
	if branch == 0 or branch == 7: p = -0.5 * D4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 - (c * s_obs) / (2 * D4))
	if branch == 1 or branch == 6: p = -0.5 * D4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 - (c * s_obs) / (2 * D4))
	if branch == 2 or branch == 5: p =  0.5 * D4 - 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 + (c * s_obs) / (2 * D4))
	if branch == 3 or branch == 4: p =  0.5 * D4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 + (c * s_obs) / (2 * D4))

	if branch > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def twofinite(obs, c, b, branch=0):
# Calculates a specific solution for the both-finite case without 
# branch deduction. This is more useful for plotting all the 
# solutions over a given interval than for finding the physical
# specular point

	s_obs = sin(obs)
	c_obs = cos(obs)
	c_obs2 = c_obs**2
	s_obs2 = s_obs**2
	c2 = c**2
	b2 = b**2

	E0 = c2 - 4 + b2 + 2 * b * c * s_obs
	E1 = 24 * b * c_obs2 * (b - c * s_obs) - 48 * (-1 + c2) * c_obs2 + E0**2 
	E2 = -432 * b2 * (-1 + c2) * c_obs2**2 + 432 * c_obs2 * (b - c * s_obs)**2 + 72 * b * c_obs2 * (b - c * s_obs) * E0 + 288 * (-1 + c2) * c_obs2 * E0 + 2 * E0**3 
	E3 = (E2 + sqrt(-4 * E1**3 + E2**2))**(1.0 / 3.0) 
	E4 = E1 / (6 * 2**(2.0 / 3.0) * E3) + E3 / (12 * 2**(1.0 / 3.0)) 
	E5 = sqrt((b2 * c_obs2) / 4.0 - E0 / 6.0 + E4)
	E6 = (b**3 * c_obs**3 - 4 * c_obs * (b - c * s_obs) - b * c_obs * E0) / (4 * E5)

	p = (b * c_obs) * 0.25 + E4 * 0.5 - sqrt(2 * c_obs2 - E3 + E6)*0.5
	if branch == 0 or branch == 7: p = (b * c_obs) * 0.25 - E4 * 0.5 - sqrt(2 * c_obs2 - E3 - E6)*0.5
	if branch == 1 or branch == 6: p = (b * c_obs) * 0.25 - E4 * 0.5 + sqrt(2 * c_obs2 - E3 - E6)*0.5
	if branch == 2 or branch == 5: p = (b * c_obs) * 0.25 + E4 * 0.5 - sqrt(2 * c_obs2 - E3 + E6)*0.5
	if branch == 3 or branch == 4: p = (b * c_obs) * 0.25 + E4 * 0.5 + sqrt(2 * c_obs2 - E3 + E6)*0.5

	if branch > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def branchdeducing_onefinite(obs, c):
# Calculates the specular point in the one-finite case,
# performing the necessary branch deductions automatically

	s_obs = sin(2 * obs)
	c_obs = cos(2 * obs)
	c2 = c**2 + 0j
	c2_4 = c2 - 4
	c2_43 = c2_4**3
	c4 = c**4

	D0 = 1152 * c2_4 * (-1 + c2 - c_obs + c2 * c_obs) + 864 * c2 * s_obs**2
	D1 = sqrt(-4 * (160 - 128 * c2 + 4 * c4 + 96 * c_obs - 96 * c2 * c_obs)**3 + (-16 * c2_43 - D0)**2)
	D2 = (16 * c2_43 + D0 - D1)**(1.0 / 3.0)
	D3 = (40 - 32 * c2 + c4 + 24 * c_obs - 24 * c2 * c_obs) / (3 * 2**(2.0 / 3.0) * D2) + D2 / (24 * 2**(1.0 / 3.0))
	D4 = sqrt(-0.25 * c2_4 + (1.0 / 12.0) * c2_4 + D3)

	p1 = acos((0.5 * D4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 + (c * s_obs) / (2 * D4))).real)
	p2 = acos((0.5 * D4 - 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - D3 + (c * s_obs) / (2 * D4))).real)
	# Use the slope of the line between (-pi/2, pi/2) to ({l}, 0) 
	# to deduce the correct branch
	l = (np.pi * (np.arctan(1 / c) - obs)) / (np.pi + 2 * np.arctan(1 / c))
	if abs(p1 - l) >= abs(p2 - l):
		return p1
	else:
		return p2
	

def branchdeducing_twofinite(obs, c, b):
# Calculates the specular point in the both-finite case,
# performing the necessary branch deductions automatically

	s_obs = sin(obs)
	c_obs = cos(obs)
	c_obs2 = c_obs**2
	s_obs2 = s_obs**2
	c2 = c**2
	b2 = b**2

	E0 = c2 - 4 + b2 + 2 * b * c * s_obs
	E1 = 24 * b * c_obs2 * (b - c * s_obs) - 48 * (-1 + c2) * c_obs2 + E0**2 
	E2 = -432 * b2 * (-1 + c2) * c_obs2**2 + 432 * c_obs2 * (b - c * s_obs)**2 + 72 * b * c_obs2 * (b - c * s_obs) * E0 + 288 * (-1 + c2) * c_obs2 * E0 + 2 * E0**3 
	E3 = (E2 + sqrt(-4 * E1**3 + E2**2))**(1.0 / 3.0) 
	E4 = E1 / (6 * 2**(2.0 / 3.0) * E3) + E3 / (12 * 2**(1.0 / 3.0)) 
	E5 = sqrt((b2 * c_obs2) / 4.0 - E0 / 6.0 + E4)
	E6 = (b**3 * c_obs**3 - 4 * c_obs * (b - c * s_obs) - b * c_obs * E0) / (4 * E5)

	l1 = -np.arctan2((c + b * sqrt(1 - b2 + c2)).real, (b2 - c2).real)
	l2 = -np.arctan2((c - b * sqrt(1 - b2 + c2)).real, (b2 - c2).real)
	# This approximation for the derivative of E7 is very robust, 
	# assuming the machine precision is larger than tolerance
	deriv = f_E7(obs, c, b) - f_E7(obs - tolerance, c, b)

	if - np.pi / 2 <= obs and obs < l1:
		return -acos(((b * c_obs) * 0.25 + E5 * 0.5 + sqrt(b2 * c_obs2 / 2.0 - E0 / 3.0 - E4 + E6)*0.5).real)
	elif l1 <= obs and obs <= l2 and deriv <= 0:
		return  acos(((b * c_obs) * 0.25 + E5 * 0.5 + sqrt(b2 * c_obs2 / 2.0 - E0 / 3.0 - E4 + E6)*0.5).real)
	elif (l2 < obs and obs < np.pi/2) or deriv > 0:
		return  acos(((b * c_obs) * 0.25 + E5 * 0.5 - sqrt(b2 * c_obs2 / 2.0 - E0 / 3.0 - E4 + E6)*0.5).real)

def numerical(obs, c, b, src=np.pi*0.5, rt=2575.0):
# Numerically determines the specular point for an arbitrary 
# configuration of source and observer for a sphere of 
# arbitrary size. c and b are mutated for notational consistency
# with the other routines.

	c = rt / c
	b = rt / b
	tolerance = 1e-6
	spec = src        # Start at the source
	t1 = np.arctan2(b * np.sin(obs) - rt * np.sin(spec), b * np.cos(obs) - rt * np.cos(spec))
	t2 = np.arctan2(c * np.sin(spec) - rt * np.sin(spec), c * np.cos(src) - rt * np.cos(spec))
	e = src - (spec - t1)
	i = src - (t2 - spec)
	ie_diff = e - i
	n = 0
	# Until the difference between incidence and emission is below the tolerance,
	# perform downhill minimization
	while abs(ie_diff) > tolerance:
		n += 1
		t1 = np.arctan2(b * np.sin(obs) - rt * np.sin(spec), b * np.cos(obs) - rt * np.cos(spec))
		t2 = np.arctan2(c * np.sin(src) - rt * np.sin(spec), c * np.cos(src) - rt * np.cos(spec))
		spec += ie_diff * 0.05
		e = src - (spec - t1)
		i = src - (t2 - spec)
		ie_diff = e - i
		if n > 100: break
	return spec

obs = np.pi*0.25
c = 0.85
b = 0.95
print(numerical(obs, c, b),
onefinite(obs, c),
twofinite(obs, c, b),
branchdeducing_onefinite(obs, c),
branchdeducing_twofinite(obs, c, b))
