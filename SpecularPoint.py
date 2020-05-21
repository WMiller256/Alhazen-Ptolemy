from cmath import sqrt, sin, cos
from math import acos

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
	if opt == 0 || opt == 7: p = -0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 - (rs * sd) / (2 * C4))
	if opt == 1 || opt == 6: p = 0.5 * C4 - 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))
	if opt == 2 || opt == 5: p = 0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))
	if opt == 3 || opt == 4: p = 0.5 * C4 + 0.5 * sqrt(-(1.0 / 3.0) * rs2_4 - C3 + (rs * sd) / (2 * C4))

	if opt > 3:
		return -acos(p.real)
	else:
		return acos(p.real)

def twofinite(d, rs, rd, opt=0):
	sd = sin(d);
	cd = cos(d);
	cd2 = cd**2;
	sd2 = sd**2;
	rs2 = rs**2;
	rd2 = rd**2;

	C0 = rd2 - 4 * cd2 + rs2 * cd2 + 2 * rd * rs * sd - 4 * sd2 + rs2 * sd2; 
	C5 = 24 * rd * cd2 * (rd - rs * sd) - 48 * (-1 + rs2) * cd2 + C0**2; 
	C6 = -432 * rd2 * (-1 + rs2) * cd2**2 + 432 * cd2 * (rd - rs * sd)**2 + 72 * rd * cd2 * (rd - rs * sd) * C0 + 288 * (-1 + rs2) * cd2 * C0 + 2 * C0**3; 
	C1 = (C6 + sqrt(-4 * C5**3 + C6**2))**(1.0 / 3.0); 
	C2 = (rd2 * cd2) / 4.0 - C0 / 6.0; 
	C3 = C5 / (6 * 2**(2.0 / 3.0) * C1) + C1 / (12 * 2**(1.0 / 3.0)); 
	C4 = sqrt(C2 + C3);

	p = (rd * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5;
	if opt == 0 or opt == 7: p = (rd * cd) * 0.25 - C4 * 0.5 - sqrt(2 * C2 - C3 - (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5;
	if opt == 1 or opt == 6: p = (rd * cd) * 0.25 - C4 * 0.5 + sqrt(2 * C2 - C3 - (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5;
	if opt == 2 or opt == 5: p = (rd * cd) * 0.25 + C4 * 0.5 - sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5;
	if opt == 3 or opt == 4: p = (rd * cd) * 0.25 + C4 * 0.5 + sqrt(2 * C2 - C3 + (rd**3 * cd**3 - 4 * cd * (rd - rs * sd) - rd * cd * C0) / (4 * C4))*0.5;

	if opt > 3:
		return -acos(p.real)
	else:
		return acos(p.real)
