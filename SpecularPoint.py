import numpy as np
from cmath import sqrt, sin, cos
from math import acos, atan

def onefinite(d, rs, opt=0):
	sd = sin(2 * d)
	cd = cos(2 * d)
	rs2 = rs**2 + 0j
	rs2_4 = rs2 - 4
	rs2_43 = rs2_4**3
	rs4 = rs**4

	C0 = 1152 * rs2_4 * (-1 + rs2 - cd + rs2 * cd) + 864 * rs2 * sd**2
	C1 = sqrt(-4 * (160 - 128 * rs2 + 4 * rs4 + 96 * cd - 96 * rs2 * cd)**3 + (-16 * rs2_43 - C0)**2)
	C2 = (16 * rs2_43 + C0 - C1)**(1.0 / 3.0)
	C3 = (40 - 32 * rs2 + rs4 + 24 * cd - 24 * rs2 * cd) / (3 * 2**(2.0 / 3.0) * C2) + C2 / (24 * 2**(1.0 / 3.0))
	C4 = sqrt(-0.25 * rs2_4 + (1.0 / 12.0) * rs2_4 + C3)

	p = -0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))
	if opt == 0 or opt == 7: p = -0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))
	if opt == 1 or opt == 6: p = -0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))
	if opt == 2 or opt == 5: p = 0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))
	if opt == 3 or opt == 4: p = 0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))

	if opt > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def twofinite(d, rs, rd, opt=0):
	sd = sin(d)
	cd = cos(d)
	cd2 = cd**2
	sd2 = sd**2
	rs2 = rs**2
	rd2 = rd**2

	C0 = rd2 - 4 * cd2 + rs2 * cd2 + 2 * rd * rs * sd - 4 * sd2 + rs2 * sd2 
	C5 = 24 * rd * cd2 * (rd - rs * sd) - 48 * (-1 + rs2) * cd2 + C0**2 
	C6 = -432 * rd2 * (-1 + rs2) * cd2**2 + 432 * cd2 * (rd - rs * sd)**2 + 72 * rd * cd2 * (rd - rs * sd) * C0 + 288 * (-1 + rs2) * cd2 * C0 + 2 * C0**3 
	C1 = (C6 + sqrt(-4 * C5**3 + C6**2))**(1.0 / 3.0) 
	C2 = (rd2 * cd2) / 4.0 - C0 / 6.0 
	C3 = C5 / (6 * 2**(2.0 / 3.0) * C1) + C1 / (12 * 2**(1.0 / 3.0)) 
	C4 = sqrt(C2 + C3)

	p = (rd * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5
	if opt == 0 or opt == 7: p = (rd * cd) * 0.25 - C4 * 0.5 - sqrt(2 * C2 - C3 - (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5
	if opt == 1 or opt == 6: p = (rd * cd) * 0.25 - C4 * 0.5 + sqrt(2 * C2 - C3 - (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5
	if opt == 2 or opt == 5: p = (rd * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5
	if opt == 3 or opt == 4: p = (rd * cd) * 0.25 + C4 * 0.5 + sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5

	if opt > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def branchdeducing_onefinite(d, rs):
	sd = sin(2 * d)
	cd = cos(2 * d)
	rs2 = rs**2 + 0j
	rs2_4 = rs2 - 4
	rs2_43 = rs2_4**3
	rs4 = rs**4

	C0 = 1152 * rs2_4 * (-1 + rs2 - cd + rs2 * cd) + 864 * rs2 * sd**2
	C1 = sqrt(-4 * (160 - 128 * rs2 + 4 * rs4 + 96 * cd - 96 * rs2 * cd)**3 + (-16 * rs2_43 - C0)**2)
	C2 = (16 * rs2_43 + C0 - C1)**(1.0 / 3.0)
	C3 = (40 - 32 * rs2 + rs4 + 24 * cd - 24 * rs2 * cd) / (3 * 2**(2.0 / 3.0) * C2) + C2 / (24 * 2**(1.0 / 3.0))
	C4 = sqrt(-0.25 * rs2_4 + (1.0 / 12.0) * rs2_4 + C3)

	p1 = acos((0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))).real)
	p2 = acos((0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))).real)
	l = (np.pi * (np.arctan(1 / rs) - d)) / (np.pi + 2 * np.arctan(1 / rs))
	if abs(p1 - l) >= abs(p2 - l):
		return p1
	else:
		return p2
	

flag = True
def branchdeducing_twofinite(d, rs, rd):
	global flag
	sd = sin(d)
	cd = cos(d)
	cd2 = cd**2
	sd2 = sd**2
	rs2 = rs**2
	rd2 = rd**2

	C0 = rd2 - 4 * cd2 + rs2 * cd2 + 2 * rd * rs * sd - 4 * sd2 + rs2 * sd2 
	C5 = 24 * rd * cd2 * (rd - rs * sd) - 48 * (-1 + rs2) * cd2 + C0**2 
	C6 = -432 * rd2 * (-1 + rs2) * cd2**2 + 432 * cd2 * (rd - rs * sd)**2 + 72 * rd * cd2 * (rd - rs * sd) * C0 + 288 * (-1 + rs2) * cd2 * C0 + 2 * C0**3 
	C1 = (C6 + sqrt(-4 * C5**3 + C6**2))**(1.0 / 3.0) 
	C2 = (rd2 * cd2) / 4.0 - C0 / 6.0 
	C3 = C5 / (6 * 2**(2.0 / 3.0) * C1) + C1 / (12 * 2**(1.0 / 3.0)) 
	C4 = sqrt(C2 + C3)

	l1 = atan(((rs + rd * sqrt(1 - rd2 + rs2)) / (rs2 - rd2)).real)
	l2 = np.pi * (7.000e-4 * sin(np.pi * (3.600 * rs + 2.0 / 3.0)) * rd**4 + 
	           3.500e-4 * sin(np.pi * (3.600 * rs + 5.0 / 3.0)) * rd**3 +
	           4.725e-4 * sin(np.pi * (3.550 * rs + 2.0 / 3.0)) * rd2 + 
	           3.750e-4 * sin(np.pi * (3.525 * rs - 1.0 / 3.0)) * rd +
	           1.550e-4 * sin(np.pi * (211.0 / 60.0 * rs + 2.0 / 3.0)))
	if flag:
		print(l1, l2.real)
		flag = False
	if - np.pi / 2 <= d and d < l1:
		return -acos(((rd * cd) * 0.25 + C4 * 0.5 + sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real)
	elif l1 <= d and d <= l2:
		return acos(((rd * cd) * 0.25 + C4 * 0.5 + sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real)
	elif d <= np.pi / 2:
		return acos(((rd * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5).real)


def numerical(td, rs, rd, ts=np.pi*0.5, rt=2575.0):
	tolerance = 1e-6
	tp = ts        # Start at the source
	t1 = np.arctan2(rd * np.sin(td) - rt * np.sin(tp), rd * np.cos(td) - rt * np.cos(tp))
	t2 = np.arctan2(rs * np.sin(tp) - rt * np.sin(tp), rs * np.cos(ts) - rt * np.cos(tp))
	e = ts - (tp - t1)
	i = ts - (t2 - tp)
	ie_diff = e - i
	n = 0
	while abs(ie_diff) > tolerance:
		n += 1
		t1 = np.arctan2(rd * np.sin(td) - rt * np.sin(tp), rd * np.cos(td) - rt * np.cos(tp))
		t2 = np.arctan2(rs * np.sin(ts) - rt * np.sin(tp), rs * np.cos(ts) - rt * np.cos(tp))
		tp += ie_diff * 0.05
		e = ts - (tp - t1)
		i = ts - (t2 - tp)
		ie_diff = e - i
		if n > 100: break

	print("Iters    : ", n)
	print("Detector : ", np.degrees(td))
	print("Scatter  : ", np.degrees(np.pi*0.5))
	print("specular : ", np.degrees(tp))
	print("abs(i-e) : ", abs(ie_diff))

import matplotlib.pyplot as plt
x = np.linspace(-np.pi*0.5, np.pi*0.5, 20)
y = [branchdeducing_onefinite(d, 0.85) for d in x]
plt.plot(x, y)
plt.show()
