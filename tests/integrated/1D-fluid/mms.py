from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t
from math import pi

# Length of the y domain
Ly = 10.

# Atomic mass number
AA = 2.0

# metric tensor
metric = Metric()  # Identity

# Define solution in terms of input x,y,z

n = 1 + 0.1*sin(2*y - t)
p = 1 + 0.1*cos(3*y + t)
mnv = AA * 0.1*sin(y + 2*t)

# Turn solution into real x and z coordinates
replace = [ (y, metric.y*2*pi/Ly) ]

n = n.subs(replace)
p = p.subs(replace)
mnv = mnv.subs(replace)

##############################
# Calculate time derivatives

nv = mnv / AA
v = nv / n
gamma = 5./3

# Density equation
dndt = - Div_par(nv)

# Pressure equation
dpdt = - Div_par(p*v) - (gamma-1.0)*p*Div_par(v)

# Momentum equation
dmnvdt = - Div_par(mnv*v) - Grad_par(p)

#############################
# Calculate sources

Sn = diff(n, t) - dndt
Sp = diff(p, t) - dpdt
Smnv = diff(mnv, t) - dmnvdt

# Substitute back to get input y coordinates
replace = [ (metric.y, y*Ly/(2*pi) ) ]

n = n.subs(replace)
p = p.subs(replace)
mnv = mnv.subs(replace)

Sn = Sn.subs(replace)
Sp = Sp.subs(replace)
Smnv = Smnv.subs(replace)

print("[n]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))

print("\n[p]")
print("solution = " + exprToStr(p))
print("\nsource = " + exprToStr(Sp))

print("\n[nv]")
print("solution = " + exprToStr(mnv))
print("\nsource = " + exprToStr(Smnv))
