import numpy as np
from cmath import sqrt, sin, cos
from math import acos, atan

tolerance = 1e-9

def f_C7(d, c, b):
	sd = sin(d)
	cd = cos(d)
	cd2 = cd**2
	sd2 = sd**2
	c2 = c**2
	b2 = b**2

	C0 = c2 - 4 + b2 + 2 * b * c * sd
	C1 = 24 * b * cd2 * (b - c * sd) - 48 * (-1 + c2) * cd2 + C0**2 
	C2 = -432 * b2 * (-1 + c2) * cd2**2 + 432 * cd2 * (b - c * sd)**2 + 72 * b * cd2 * (b - c * sd) * C0 + 288 * (-1 + c2) * cd2 * C0 + 2 * C0**3 
	C3 = (C2 + sqrt(-4 * C1**3 + C2**2))**(1.0 / 3.0) 
	C4 = C1 / (6 * 2**(2.0 / 3.0) * C3) + C3 / (12 * 2**(1.0 / 3.0)) 
	C5 = sqrt((b2 * cd2) / 4.0 - C0 / 6.0 + C4)	

	return  ((b2 * cd2) / 2 - C0 / 3.0 - C4 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5)).real
	print(cd)
	print(sd)
	print(cd2)
	print(sd2)
	print(C0)
	print(C1)
	print(C2)
	print(C3)
	print(C4)
	print(C5)

def onefinite(d, c, opt=0):
	sd = sin(2 * d)
	cd = cos(2 * d)
	c2 = c**2 + 0j
	c2_4 = c2 - 4
	c2_43 = c2_4**3
	c4 = c**4

	C0 = 1152 * c2_4 * (-1 + c2 - cd + c2 * cd) + 864 * c2 * sd**2
	C1 = sqrt(-4 * (160 - 128 * c2 + 4 * c4 + 96 * cd - 96 * c2 * cd)**3 + (-16 * c2_43 - C0)**2)
	C2 = (16 * c2_43 + C0 - C1)**(1.0 / 3.0)
	C3 = (40 - 32 * c2 + c4 + 24 * cd - 24 * c2 * cd) / (3 * 2**(2.0 / 3.0) * C2) + C2 / (24 * 2**(1.0 / 3.0))
	C4 = sqrt(-0.25 * c2_4 + (1.0 / 12.0) * c2_4 + C3)

	p = -0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 - (c * sd) / (2 * C4))
	if opt == 0 or opt == 7: p = -0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 - (c * sd) / (2 * C4))
	if opt == 1 or opt == 6: p = -0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 - (c * sd) / (2 * C4))
	if opt == 2 or opt == 5: p = 0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 + (c * sd) / (2 * C4))
	if opt == 3 or opt == 4: p = 0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 + (c * sd) / (2 * C4))

	if opt > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def twofinite(d, c, b, opt=0):
	sd = sin(d)
	cd = cos(d)
	cd2 = cd**2
	sd2 = sd**2
	c2 = c**2
	b2 = b**2

	C0 = c2 - 4 + b2 + 2 * b * c * sd
	C5 = 24 * b * cd2 * (b - c * sd) - 48 * (-1 + c2) * cd2 + C0**2 
	C6 = -432 * b2 * (-1 + c2) * cd2**2 + 432 * cd2 * (b - c * sd)**2 + 72 * b * cd2 * (b - c * sd) * C0 + 288 * (-1 + c2) * cd2 * C0 + 2 * C0**3 
	C1 = (C6 + sqrt(-4 * C5**3 + C6**2))**(1.0 / 3.0) 
	C2 = (b2 * cd2) / 4.0 - C0 / 6.0 
	C3 = C5 / (6 * 2**(2.0 / 3.0) * C1) + C1 / (12 * 2**(1.0 / 3.0)) 
	C4 = sqrt(C2 + C3)

	p = (b * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C4))*0.5
	if opt == 0 or opt == 7: p = (b * cd) * 0.25 - C4 * 0.5 - sqrt(2 * C2 - C3 - (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C4))*0.5
	if opt == 1 or opt == 6: p = (b * cd) * 0.25 - C4 * 0.5 + sqrt(2 * C2 - C3 - (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C4))*0.5
	if opt == 2 or opt == 5: p = (b * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C4))*0.5
	if opt == 3 or opt == 4: p = (b * cd) * 0.25 + C4 * 0.5 + sqrt(2 * C2 - C3 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C4))*0.5

	if opt > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def branchdeducing_onefinite(d, c):
	sd = sin(2 * d)
	cd = cos(2 * d)
	c2 = c**2 + 0j
	c2_4 = c2 - 4
	c2_43 = c2_4**3
	c4 = c**4

	C0 = 1152 * c2_4 * (-1 + c2 - cd + c2 * cd) + 864 * c2 * sd**2
	C1 = sqrt(-4 * (160 - 128 * c2 + 4 * c4 + 96 * cd - 96 * c2 * cd)**3 + (-16 * c2_43 - C0)**2)
	C2 = (16 * c2_43 + C0 - C1)**(1.0 / 3.0)
	C3 = (40 - 32 * c2 + c4 + 24 * cd - 24 * c2 * cd) / (3 * 2**(2.0 / 3.0) * C2) + C2 / (24 * 2**(1.0 / 3.0))
	C4 = sqrt(-0.25 * c2_4 + (1.0 / 12.0) * c2_4 + C3)

	p1 = acos((0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 + (c * sd) / (2 * C4))).real)
	p2 = acos((0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * c2_4 - C3 + (c * sd) / (2 * C4))).real)
	l = (np.pi * (np.arctan(1 / c) - d)) / (np.pi + 2 * np.arctan(1 / c))
	if abs(p1 - l) >= abs(p2 - l):
		return p1
	else:
		return p2
	

def branchdeducing_twofinite(d, c, b):
	sd = sin(d)
	cd = cos(d)
	cd2 = cd**2
	sd2 = sd**2
	c2 = c**2
	b2 = b**2

	C0 = c2 - 4 + b2 + 2 * b * c * sd
	C1 = 24 * b * cd2 * (b - c * sd) - 48 * (-1 + c2) * cd2 + C0**2 
	C2 = -432 * b2 * (-1 + c2) * cd2**2 + 432 * cd2 * (b - c * sd)**2 + 72 * b * cd2 * (b - c * sd) * C0 + 288 * (-1 + c2) * cd2 * C0 + 2 * C0**3 
	C3 = (C2 + sqrt(-4 * C1**3 + C2**2))**(1.0 / 3.0) 
	C4 = C1 / (6 * 2**(2.0 / 3.0) * C3) + C3 / (12 * 2**(1.0 / 3.0)) 
	C5 = sqrt((b2 * cd2) / 4.0 - C0 / 6.0 + C4)
	C6 = (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5)

	l1 = -np.arctan2((c + b * sqrt(1 - b2 + c2)).real, (b2 - c2).real)
	l2 = -np.arctan2((c - b * sqrt(1 - b2 + c2)).real, (b2 - c2).real)
	deriv = f_C7(d, c, b) - f_C7(d - tolerance, c, b)

	if - np.pi / 2 <= d and d < l1:
#		print('b: one  d: {: 6.4f} deriv: {: 6.4e}, l1: {: 6.4f}, l2: {: 6.4f}'.format(d, deriv, l1, l2))
		return -acos(((b * cd) * 0.25 + C5 * 0.5 + sqrt(b2 * cd2 / 2.0 - C0 / 3.0 - C4 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5))*0.5).real)
	elif l1 <= d and d <= l2 and deriv <= 0:
#		print('b: two  d: {: 6.4f} deriv: {: 6.4e}, l1: {: 6.4f}, l2: {: 6.4f}'.format(d, deriv, l1, l2))
		return  acos(((b * cd) * 0.25 + C5 * 0.5 + sqrt(b2 * cd2 / 2.0 - C0 / 3.0 - C4 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5))*0.5).real)
	elif (l2 < d and d < np.pi/2) or deriv > 0:
#		print('b: thr  d: {: 6.4f} deriv: {: 6.4e}, l1: {: 6.4f}, l2: {: 6.4f}'.format(d, deriv, l1, l2))
		return  acos(((b * cd) * 0.25 + C5 * 0.5 - sqrt(b2 * cd2 / 2.0 - C0 / 3.0 - C4 + (b**3 * cd**3 - 4 * cd * (b - c * sd) - b * cd * C0) / (4 * C5))*0.5).real)

def numerical(td, c, b, ts=np.pi*0.5, rt=2575.0):
	tolerance = 1e-6
	tp = ts        # Start at the source
	t1 = np.arctan2(b * np.sin(td) - rt * np.sin(tp), b * np.cos(td) - rt * np.cos(tp))
	t2 = np.arctan2(c * np.sin(tp) - rt * np.sin(tp), c * np.cos(ts) - rt * np.cos(tp))
	e = ts - (tp - t1)
	i = ts - (t2 - tp)
	ie_diff = e - i
	n = 0
	while abs(ie_diff) > tolerance:
		n += 1
		t1 = np.arctan2(b * np.sin(td) - rt * np.sin(tp), b * np.cos(td) - rt * np.cos(tp))
		t2 = np.arctan2(c * np.sin(ts) - rt * np.sin(tp), c * np.cos(ts) - rt * np.cos(tp))
		tp += ie_diff * 0.05
		e = ts - (tp - t1)
		i = ts - (t2 - tp)
		ie_diff = e - i
		if n > 100: break

import matplotlib.pyplot as plt
x = np.linspace(-np.pi*0.5, np.pi*0.5, 200)
y = [branchdeducing_twofinite(d, 0.85, 0.95) for d in x]
plt.plot(x, y)
plt.show()
