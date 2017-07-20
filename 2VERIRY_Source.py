import math as m
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys

R_UNIVERSAL = 8314.472; # J/kmol/K
Wc = 29. # kg/kmol
Wv = 18. # kg/kmol
Wr = Wc/Wv

def MD(d,rhol):
  return rhol*m.pi/6.*d**3

def calcHvap(T): 
  #latent heat of vaporization in J/Kg...
  AL = 2.257E6;
  BL = 2.595E3;
  CL = 373.15;
  DL = 1.0;
  nL = 1;
  Hvap = AL+BL*pow((CL-T/DL),nL);
  #print "Hvap: ", Hvap
  return( Hvap ); #  J/kg

def calcPv(T,pg):
  Tboil = calcTboil(pg)
  pv = pg*np.exp(calcHvap(T)*Wv/R_UNIVERSAL*(1.0/Tboil-1.0/T))
  #print "pv: ", pv, " , Tboil:", Tboil, " , T:", T
  #print "ext arg: ", calcHvap(T)*Wv/R_UNIVERSAL*(1.0/Tboil-1.0/T)
  return(pv);

def FindD(Mw,N,rhol):
  d = (Mw/(N*rhol*m.pi/6.))**(1./3.)
  return d;

def FindND(Mw,D,rhol):
  nd = Mw/(rhol*m.pi/6.*D**3)
  return nd

def FindD_mdot(N , dmw , rhol , d_old):
  assert(dmw >= 0),"dmw should be positive"
  print "d_old**3: " , d_old**3, " , dmw/(N*rhol/6.*m.pi: " , dmw/(N*rhol/6.*m.pi)
  d = (d_old**3 - dmw/(N*rhol/6.*m.pi ))**(1./3.)
  return d;

def calcTboil(P):
  # using Clausius-Clapeyron eqn: dP/dT = ds/dv = L/Tdv
  # assuming Hvap is constant, we use its value at T0 = 373K
  Hvap = calcHvap(373) # must have the same unit as Rv*T (J/kg)
  Rv = R_UNIVERSAL/Wv # J/kg/K
  T0 = 373.
  P0 = 1e5
  Tb = 1.0/(1.0/T0 - Rv*m.log(P/P0)/Hvap)
  return Tb

fm = open("Source_xme.dat",'r')
line = fm.readline();
line = line.strip();
col = line.split();
ilm = int(col[0]);
print "ilm: " , ilm 

x = np.zeros(ilm)
sm = np.zeros(ilm)
se = np.zeros(ilm)

ix = 0
for line in fm:
  line = line.strip()
  col = line.split()
  x[ix] = float(col[0])
  sm[ix] = float(col[1])
  se[ix] = float(col[2])
  ix += 1

print "x: " , x
print "sm: " , sm
print "se: " , se

Nx = ilm
dx = x[2] - x[1]
print "dx: ", dx, " , x[ilm-1]", x[ilm-1]

def getSM(x0):
  s0 = -1.0
  for ix in range(len(x)-1):
    if (x0 >= x[ix] and x0 < x[ix+1]):
      dx = x[ix+1] - x[ix]
      d = x0 - x[ix]
      s0 = sm[ix]*(1.0-d/dx) + sm[ix+1]*d/dx
      break
  assert(s0>0), "couldn't find x point in getFv!"
  return s0

def getSE(x0):
  s0 = -1.0
  for ix in range(len(x)-1):
    if (x0 >= x[ix] and x0 < x[ix+1]):
      dx = x[ix+1] - x[ix]
      d = x0 - x[ix]
      s0 = sm[ix]*(1.0-d/dx) + sm[ix+1]*d/dx
      break
  assert(s0>0), "couldn't find x point in getQ!"
  return s0

time.sleep(5)
#================
# inputs
#================
pg0 = 9.0873e5
Tg0 = 1395.43
mdot_g = 25.5867 # kg/s
gamma = 1.29843

Rg = R_UNIVERSAL/Wc
Rg = 292.066
Rv = R_UNIVERSAL/Wv
rhog0 = pg0/Rg/Tg0
cpg = gamma/(gamma-1.0)*Rg
cpv = 2000.0 # water vapor cp

#Area = np.zeros(Nx)
#A0 = 0.249861
#A1 = 0.519
#A2 = 0.519
#A3 = 0.1297
#A4 = 0.1297
#x0 = 0.0
#x1 = 0.515 # m
#x2 = 1.010 # m
#x3 = 1.365 # m
#x4 = 2.75
#slope0 = (A1-A0)/(x1-x0)
#slope1 = (A2-A1)/(x2-x1)
#slope2 = (A3-A2)/(x3-x2)
#slope3 = (A4-A3)/(x4-x3)
#for ix in range(Nx):
#  if (x[ix]>=x0 and x[ix]<x1):
#    Area[ix] = A0 + (x[ix]-x0)*slope0
#  elif (x[ix]>=x1 and x[ix]<x2):
#    Area[ix] = A1 + (x[ix]-x1)*slope1
#  #elif (x[ix]>=x2):
#  #  Area[ix] = A2 + (x[ix]-x2)*slope2
#  elif (x[ix]>=x2 and x[ix]<x3):
#    Area[ix] = A2 + (x[ix]-x2)*slope2
#  elif (x[ix]>=x3 and x[ix]<x4):
#    Area[ix] = A3 + (x[ix]-x3)*slope3
#  else:
#    print "not enough options in the slopes!"
#    sys.exit()
#  assert(Area[ix]>0.0),"Area should be possitive"

Area = np.zeros(Nx)
Lx = 1.0
Amin = 0.1936
Amax = Amin*1.0
Lmax = 0.05
IXR = Nx*Lmax/Lx
for ix in range(Nx):
  Area[ix] = Amin + ix*(Amax-Amin)/IXR
  Area[ix] = min(Area[ix],Amax)

Vx0 = mdot_g/rhog0/Area[0]
print "Area: ", Area

# liquid and vapor
Mv = np.zeros(Nx)

# compressible gas
ug = np.zeros(Nx)
rhog = np.zeros(Nx)
pg = np.zeros(Nx)
Tg = np.zeros(Nx)
q = np.zeros(Nx)
C = np.zeros(Nx)

#===============
# initializing
#===============
ug[0] = Vx0
pg[0] = pg0
Tg[0] = Tg0
rhog[0] = rhog0
print "pg0: " ,pg[0], " , Tg0: " , Tg[0], " , mdot_g", mdot_g, " , Vx0: " , Vx0, " , rhog0" , rhog[0], " , cpg: " , cpg, " , cpv: " , cpv

Mv[0] = 0.0

# Source terms
Sm = np.zeros(Nx)
Su = np.zeros(Nx)
Se = np.zeros(Nx)
IX = ilm - 1
for ix in range(ilm-1):

  #Fv = getSM(x[ix])*(Area[ix]*dx)
  #q[ix] = getSE(x[ix])*(Area[ix]*dx)
  Fv = (sm[ix])*(Area[ix]*dx)
  q[ix] = (-se[ix])*(Area[ix]*dx)/10.0
  Mv[ix+1] = Mv[ix] + Fv*dx/ug[ix]

  # cp_mix calc
  m_mix = rhog[ix]*Area[ix]*dx
  mass_gas = m_mix - Mv[ix] 
  #print "m_mix: " , m_mix, " , Mv[ix]: ", Mv[ix], " , mass_gas: " , mass_gas
  assert(mass_gas>0.0), "mg should be positive!"
  Y_g = mass_gas/m_mix
  Y_v = Mv[ix]/m_mix
  cp_mix = cpg*Y_g + cpv*Y_v
  W_mix_inv = Y_g/Wc + Y_v/Wv
  R_mix = R_UNIVERSAL*W_mix_inv
  #print "R_mix: " , R_mix, " , cp_mix: " , cp_mix
  a1 = rhog[ix]*ug[ix]*Area[ix] + Fv
  a2 = pg[ix]*Area[ix] + rhog[ix]*ug[ix]**2*Area[ix] + (Area[ix+1]-Area[ix])*pg[ix]
  a3 = rhog[ix]*ug[ix]*Area[ix]*(cp_mix*Tg[ix] + (ug[ix]**2)/2.0) - q[ix]
  #a3 = rhog[ix]*ug[ix]*Area[ix]*cp_mix*Tg[ix] - q[ix]

  delta = a2**2 - 4.0*(a1*(1.0 - R_mix/(2.0*cp_mix))-Fv)*(R_mix/cp_mix)*a3*a1/(a1-Fv)
  #delta = a2**2 - 4.0*(a1-Fv)*(R_mix/cp_mix)*a3*a1/(a1-Fv)
  ug[ix+1] = (a2 - m.sqrt(delta))/(2.0*(a1*(1.0-R_mix/(2.0*cp_mix)) - Fv))
  #ug[ix+1] = (a2 - m.sqrt(delta))/(2.0*(a1 - Fv))
  rhog[ix+1] = a1/ug[ix+1]/Area[ix+1]
  Tg[ix+1] = (a3/(a1-Fv) - ug[ix+1]**2/2.0)/cp_mix
  #Tg[ix+1] = (a3/(a1-Fv))/cp_mix
  pg[ix+1] = rhog[ix+1]*R_mix*Tg[ix+1]

C = np.sqrt(gamma*pg/rhog)
M = ug/C

print "mdot_g: " , mdot_g, " , cpg: " , cpg, " , Tg0: " , Tg0, " , cpv: " , cpv
print "Tg[end]", Tg[IX]

#=============
# I/O section
#=============

plt.figure()
plt.plot(x[0:IX],Mv[0:IX],label="Vapor mass")
plt.legend()
plt.xlabel("x(m)")
plt.ylabel("vapor mass(kg)")
plt.show(block=False)
plt.savefig("mass.png")

plt.figure()
plt.plot(x[0:IX],ug[0:IX])
plt.xlabel("x(m)")
plt.ylabel("ug(m/s)")
plt.show(block=False)
plt.savefig("ug.png")

plt.figure()
plt.plot(x[0:IX],rhog[0:IX])
plt.xlabel("x(m)")
plt.ylabel("rhog(kg/m^3)")
plt.show(block=False)
plt.savefig("rhog.png")

plt.figure()
plt.plot(x[0:IX],Tg[0:IX])
plt.xlabel("x(m)")
plt.ylabel("Tg(K)")
plt.show(block=False)
plt.savefig("Tg.png")

plt.figure()
plt.plot(x[0:IX],pg[0:IX])
plt.xlabel("x(m)")
plt.ylabel("pg(Pa)")
plt.show(block=False)
plt.savefig("pg.png")

plt.figure()
plt.plot(x[1:IX],q[1:IX])
plt.xlabel("x(m)")
plt.ylabel("q(J/s)")
plt.show(block=False)
plt.savefig("q.png")

plt.figure()
plt.plot(x[0:IX],C[0:IX])
plt.xlabel("x(m)")
plt.ylabel("C(m/s)")
plt.show(block=False)
plt.savefig("C.png")

plt.figure()
plt.plot(x[0:IX],Area[0:IX])
plt.xlabel("x(m)")
plt.ylabel("Area(m^2)")
plt.show(block=False)
plt.savefig("Area.png")

plt.figure()
plt.plot(x[0:IX],M[0:IX])
plt.xlabel("x(m)")
plt.ylabel("Mach #")
plt.show(block=False)
plt.savefig("Mach.png")

plt.figure()
plt.plot(x[0:IX],Sm[0:IX])
plt.xlabel("x(m)")
plt.ylabel("mass source (kg/s/m^3)")
plt.show(block=False)

plt.figure()
plt.plot(x[1:IX],Se[1:IX])
plt.xlabel("x(m)")
plt.ylabel("energy source (kg/s/m^3)")
plt.show()
