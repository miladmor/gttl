import math as m
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys
import read_probe as rp

R_UNIVERSAL = 8314.472; # J/kmol/K
Wc = 29. # kg/kmol
Wv = 18. # kg/kmol
#Wv = Wc # kg/kmol
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

#print "P: ", 1e5, " , Tb: " , calcTboil(1e5), " P: ", 10e5, " , Tb: " , calcTboil(10e5), " P: ", 20e5, " , Tb: " , calcTboil(20e5)
#time.sleep(10)
#print "calcHvap(460): ", calcHvap(460)

Nx = 20000
#Lx = 1.365
#Lx = 2.75
#Lx = 1.
Lx = 1.85
#Lx = 1.461957 # almost chokes with this Lx and no D increase
x = np.array(range(Nx))
x = x*Lx/Nx
dx = x[2] - x[1]
print "dx: ", dx

#================
# inputs
#================
pg0 = 9.0873e5
Tg0 = 1395.43
mdot_g = 25.5867 # kg/s
gamma = 1.29843

mdot_w0 = 14.35 # kg/s
Tl0 = 294.
rhol = 1000
cl = 4184.0 # J/(kg*K)
We_cr = 10
sigma = 0.072

Rg = R_UNIVERSAL/Wc
Rg = 292.066
Rv = R_UNIVERSAL/Wv
rhog0 = pg0/Rg/Tg0
cpg = gamma/(gamma-1.0)*Rg
cpv = 2000.0 # water vapor cp
#cpv = cpg # water vapor cp

Area = np.zeros(Nx)
A0 = 0.249861
A1 = 0.519
A2 = 0.519
A3 = 0.1297
A4 = 0.1297
x0 = 0.0
x1 = 0.515 # m
x2 = 1.010 # m
x3 = 1.365 # m
x4 = 2.75
slope0 = (A1-A0)/(x1-x0)
slope1 = (A2-A1)/(x2-x1)
slope2 = (A3-A2)/(x3-x2)
slope3 = (A4-A3)/(x4-x3)
for ix in range(Nx):
  if (x[ix]>=x0 and x[ix]<x1):
    Area[ix] = A0 + (x[ix]-x0)*slope0
  elif (x[ix]>=x1 and x[ix]<x2):
    Area[ix] = A1 + (x[ix]-x1)*slope1
  #elif (x[ix]>=x2):
  #  Area[ix] = A2 + (x[ix]-x2)*slope2
  elif (x[ix]>=x2 and x[ix]<x3):
    Area[ix] = A2 + (x[ix]-x2)*slope2
  elif (x[ix]>=x3 and x[ix]<x4):
    Area[ix] = A3 + (x[ix]-x3)*slope3
  else:
    print "not enough options in the slopes!"
    sys.exit()
  assert(Area[ix]>0.0),"Area should be possitive"

##########################################
# area, length of domain (Lx), length of probe, file of probe should be checked when comparing charles to this 
#Area = np.zeros(Nx)
#Amin = 0.1936
#Amax = Amin*1.0
#Lmax = 0.05
#IXR = Nx*Lmax/Lx
#for ix in range(Nx):
#  Area[ix] = Amin + ix*(Amax-Amin)/IXR
#  Area[ix] = min(Area[ix],Amax)
##########################################

Vx0 = mdot_g/rhog0/Area[0]
print "Area: ", Area

# liquid and vapor
D = np.zeros(Nx)
Mw = np.zeros(Nx)
Mv = np.zeros(Nx)
mdot_w = np.zeros(Nx)
mdot_v = np.zeros(Nx)
Tl = np.zeros(Nx)

# compressible gas
ug = np.zeros(Nx)
rhog = np.zeros(Nx)
pg = np.zeros(Nx)
Tg = np.zeros(Nx)
q = np.zeros(Nx)
q_mod = np.zeros(Nx)
Q = np.zeros(Nx)
C = np.zeros(Nx)
cpmix = np.zeros(Nx)

#===============
# initializing
#===============
ug[0] = Vx0
pg[0] = pg0
Tg[0] = Tg0
rhog[0] = rhog0
cpmix[0] = cpg
print "pg0: " ,pg[0], " , Tg0: " , Tg[0], " , mdot_g", mdot_g, " , Vx0: " , Vx0, " , rhog0" , rhog[0], " , cpg: " , cpg, " , cpv: " , cpv

D[0] = We_cr/(rhog0*Vx0**2/sigma)
mdot_w[0] = mdot_w0
mdot_v[0] = 0.0
Mw[0] = mdot_w0*dx/Vx0
Mv[0] = 0.0
Tl[0] = Tl0
ND = m.floor(Mw[0]/MD(D[0],rhol))
print "mdot_w0: " , mdot_w0, " , We_cr: " , rhog0*Vx0**2*D[0]/sigma, " , D[0]: ", D[0] , " , ND: " , ND, "Mw[0]", Mw[0]

#Tr = Tg
Tr = calcTboil(pg0)
mug = 6.109e-06 +  4.604e-08*Tr -  1.051e-11*Tr*Tr; # kg/m/s
#kg = 3.227e-03 + 8.3894e-05*Tr - 1.9858e-08*Tr*Tr; # W/m/K
#if (Tr > 600):
#  cpg = (0.647 + 5.5e-05*Tr)*kg/mug;
#else:
#  cpg = (0.815 - 4.958e-04*Tr + 4.514e-07*Tr*Tr)*kg/mug; # J/kg/K

# Source terms
Sm = np.zeros(Nx)
Su = np.zeros(Nx)
Se = np.zeros(Nx)

# checks and balances
m_evap = 0.0
mg = np.zeros(Nx)
Q_gas = 0.0
Q_wat = 0.0
m_gas = 0.0
mdot_water = np.zeros(Nx)
mdot_vapor = np.zeros(Nx)
Q_LH = 0.0
Q_VH = 0.0

IX = Nx-1
for ix in range(Nx-1):
  #print "ix: " , ix, " , Tg: " , Tg[ix]
  #print "md: " , MD(D[ix],rhol)
  md = MD(D[ix],rhol)
  tau = (rhol*D[ix]**2)/(18.*mug)
  Scg = 1.
  Red = rhog[ix]*ug[ix]*D[ix]/mug
  Sh = 2. + 0.552*Red**(0.5)*Scg**(1./3.)
  pv = calcPv(Tl[ix],pg[ix])
  Xs = pv/pg[ix]
  Ys = Xs/(Xs + (1.-Xs)*Wr)
  Bm = Ys/(1.-Ys)
  #print "Bm: ", Bm, " Xs: " , Xs , " Ys: " , Ys, " pv: " , pv
  if Xs > 1.0:
    IX = ix
    break
  Hm = m.log(1+Bm)
  k = Sh/(3.*Scg)*Hm
  mdot = ND*k/tau*md
  mdot_w[ix+1] = mdot_w[ix] - mdot
  mdot_v[ix+1] = mdot_v[ix] + mdot
  # used ug[ix] instead of ug[ix+1], this is a first order approximation
  Mw[ix+1] = mdot_w[ix+1]*dx/ug[ix]
  Mv[ix+1] = mdot_v[ix+1]*dx/ug[ix]
  #Mw[ix+1] = max(0.0,Mw[ix+1])
  mdot_w[ix+1] = max(0.0,mdot_w[ix+1])
  if mdot_w[ix+1] == 0.0:
    IX = ix
    break
  D[ix+1] = FindD(Mw[ix+1] , ND , rhol)
  D[ix+1] = min(D[ix+1],D[ix])
  ND = FindND(Mw[ix+1], D[ix+1], rhol)
  #dmw = Mw[ix] - Mw[ix+1]
  #print "ix: " , ix, " , D[ix]: " , D[ix] , " , md: " , md
  #print "tau: " , tau, " , Red: " , Red, " , Sh: " , Sh
  #print "mdot_w[ix]: ", mdot_w[ix]," , mdot_v[ix]: ", mdot_v[ix]
  #print "mdot: ", mdot , " , mdot_w[ix+1]: ", mdot_w[ix+1], " , Mw[ix+1]: " , Mw[ix+1], " , mdot_v[ix+1]: ", mdot_v[ix+1], " , Mv[ix+1]: ", Mv[ix+1], " , dmw: ", dmw
  #print "D[ix+1]: ", D[ix+1]
  #print "D[ix+1]: " , D[ix+1], " Mw[ix+1]: ", Mw[ix+1], " dMw: ", Mw[ix+1]-Mw[ix], " tau: ",  tau
  #print "md:", md, " tau: ",  tau, " Scg: ", Scg, " Red: ", Red, " Sh: " , Sh, " pv: " , pv, " Xs: ", Xs, " Ys: ", Ys, " Bm: ", Bm, " Hm: " , Hm, " k: ", k
  
  Prg = 0.7
  Nu = 2. + 0.552*Red**(0.5)*Prg**(1./3.)
  f2 = 1.
  m_mix = rhog[ix]*Area[ix]*dx
  mass_gas = m_mix - Mv[ix]
  assert(mass_gas>0.0), "mg should be positive!"
  cp_mix = (cpg*mass_gas + cpv*Mv[ix])/m_mix
  rhsTl_d = f2*Nu/3./Prg*(cp_mix/cl)/tau*(-1.0)/ug[ix]
  rhsTl = f2*Nu/3./Prg*(cp_mix/cl)/tau*(Tg[ix])/ug[ix] - calcHvap(Tl[ix])/cl*(k/tau)/ug[ix]
  #rhsTl = f2*Nu/3./Prg*(cp_mix/cl)/tau*(Tg[ix])/ug[ix] - calcHvap(350)/cl*(k/tau)/ug[ix]
  Tl[ix+1] = (Tl[ix] + dx*rhsTl)/(1.0 - dx*rhsTl_d)

  #q[ix+1] = calcHvap(350)*mdot # J/s
  #q[ix+1] = calcHvap(Tl[ix])*mdot # J/s
  q[ix+1] = f2*Nu/3./Prg*(cp_mix/cl)/tau*(Tg[ix]-Tl[ix+1])*Mw[ix+1]*cl
  q[ix+1] += mdot*cpv*(Tg[ix]-Tl[ix+1])
  #q[ix+1] += mdot*cpv*(Tg[ix]-Tl[0])

  Q[ix+1] = Q[ix] + q[ix+1]
  Q_LH += f2*Nu/3./Prg*(cp_mix/cl)/tau*(Tg[ix]-Tl[ix+1])*Mw[ix+1]*cl
  Q_VH += mdot*cpv*(Tg[ix]-Tl[ix+1])

  if Mw[ix+1] == 0.0:
    IX = ix
    break
  Fv = mdot

  # cpg constant calc
  #a1 = rhog[ix]*ug[ix]*Area[ix] + Fv
  #a2 = pg[ix]*Area[ix] + rhog[ix]*ug[ix]**2*Area[ix] + (Area[ix+1]-Area[ix])*pg[ix]
  #a3 = rhog[ix]*ug[ix]*Area[ix]*(cpg*Tg[ix] + (ug[ix]**2)/2.0) - q[ix]

  ##delta = a2**2 - 4.0*(a1*(gamma+1.0)/(2.0*gamma)-Fv)*(gamma-1.0)/gamma*a3
  ##ug[ix+1] = (a2 - m.sqrt(delta))/((a1*(gamma+1.0)/gamma) - 2.0*Fv)
  #delta = a2**2 - 4.0*(a1*(1.0 - Rg/(2.0*cpg))-Fv)*(Rg/cpg)*a3
  #ug[ix+1] = (a2 - m.sqrt(delta))/(2.0*(a1*(1.0-Rg/(2.0*cpg)) - Fv))
  #rhog[ix+1] = a1/ug[ix+1]/Area[ix+1]
  #Tg[ix+1] = (a3/a1 - ug[ix+1]**2/2.0)/cpg
  #pg[ix+1] = rhog[ix+1]*Rg*Tg[ix+1]

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

  # Tg implicit calc
  # q should only include the explicit parts
  #a1 = rhog[ix]*ug[ix]*Area[ix] + Fv
  #a2 = pg[ix]*Area[ix] + rhog[ix]*ug[ix]**2*Area[ix] + (Area[ix+1]-Area[ix])*pg[ix]
  #a3 = rhog[ix]*ug[ix]*Area[ix]*(cpg*Tg[ix] + (ug[ix]**2)/2.0) - q[ix] + mdot*cpv*(Tl[ix])

  #cp_avg = cpg + mdot/a1*cpv
  ##cp_avg = cpg
  #delta = a2**2 - 4.0*(a1*(1.0 - Rg/(2.0*cp_avg))-Fv)*(Rg/cp_avg)*a3
  #ug[ix+1] = (a2 - m.sqrt(delta))/(2.0*(a1*(1.0-Rg/(2.0*cp_avg)) - Fv))
  #rhog[ix+1] = a1/ug[ix+1]/Area[ix+1]
  #Tg[ix+1] = (a3/a1 - ug[ix+1]**2/2.0)/cp_avg
  #pg[ix+1] = rhog[ix+1]*Rg*Tg[ix+1]
  ##cpmix[ix+1] = (rhog[ix+1]*Area[ix+1]*dx)/(rhog[ix+1]*Area[ix+1]*dx + Mv[ix+1])*cpg + Mv[ix+1]/(rhog[ix+1]*Area[ix+1]*dx + Mv[ix+1])*cpv

  #print "rhog[ix+1]: " , rhog[ix+1], " , ug[ix+1]: " , ug[ix+1], " , Tg[ix+1]: " , Tg[ix+1], " , pg[ix+1]: " , pg[ix+1]

  # Source terms
  heat_flux_2 = rhog[ix+1]*ug[ix+1]*Area[ix+1]*(cpg*Tg[ix+1]+ug[ix+1]**2/2.0)
  heat_flux_1 = rhog[ix]*ug[ix]*Area[ix]*(cpg*Tg[ix]+ug[ix]**2/2.0)
  q_mod[ix] = (heat_flux_2 - heat_flux_1)
  Sm[ix] = mdot/(Area[ix]*dx) # kg/s/vol
  Se[ix] = q_mod[ix]/(Area[ix]*dx) # J/s/vol
  #Se[ix] = -q[ix]/(Area[ix]*dx) # J/s/vol

  # checks and balances
  m_evap += mdot*dx/ug[ix]
  mg[ix] = rhog[ix]*Area[ix]*ug[ix]
  Q_gas += rhog[ix+1]*Area[ix+1]*dx*cpg*(Tg[ix+1]-Tg[0]) + rhog[ix+1]*Area[ix+1]*dx*(ug[0]**2 - ug[ix+1]**2)*0.5
  Q_wat += q[ix+1]*rhog[ix+1]*Area[ix+1]*dx
  m_gas += rhog[ix]*Area[ix]*dx
  mdot_water[ix] = Mw[ix]/(dx/ug[ix])
  mdot_vapor[ix] = Mv[ix]/(dx/ug[ix])
  
  #sys.stdin.read(1)

# * (rhog[ix]*Area[ix]*dx)
#calcHvap(Tl[ix])*mdot*dx/ug[ix]

print "m_evap: " , m_evap, " , Mw[0]: " , Mw[0]
print "Q_gas: ", Q_gas 
print "Q_wat: " , Q_wat
#print "latent heat: " , m_evap*calcHvap(Tl[0])
print "m_gas: ", m_gas
print  "Q_gas/m_gas: ",Q_gas/m_gas

C = np.sqrt(gamma*pg/rhog)
M = ug/C
heat_flow = rhog*ug*Area*(cpg*Tg+ug**2/2.0)
print "heat_flow[0]: " , heat_flow[0], " , heat_flow[IX]: ", heat_flow[IX], " , DH: " , heat_flow[IX] - heat_flow[0]
#for ix in range(len(q_mod)-1):
#  q_mod[ix] = (heat_flow[ix+1] - heat_flow[ix])/(Area[ix]*dx)

Tl_mid = 420
Q_total = mdot_w[0]*calcHvap(350)
intake_E = mdot_g*(cpg*Tg[0]) + mdot_w[0]*(cpv*Tl[0])
rhs = intake_E - Q_total
Tf = (rhs)/(mdot_g*cpg+mdot_w[0]*cpv)
print "mdot_g: " , mdot_g, " , cpg: " , cpg, " , Tg0: " , Tg0, " , mdot_w[0]: " ,mdot_w[0], " , cpv: " , cpv, " , Tl[0]: ", Tl[0], " , calcHvap(Tl[0]): " , calcHvap(350)
print "Tf_theory: " , Tf, " , Tg[end]", Tg[IX]

#=============
# I/O section
#=============
T_E = 0.0
T_m = 0.0
T_E_mod = 0.0
for ix in range(IX):
  T_E += Se[ix]*(Area[ix]*dx)
  T_m += Sm[ix]*(Area[ix]*dx)
  T_E_mod += q_mod[ix]*(Area[ix]*dx)
print "T_E: " ,T_E, " , T_m: " , T_m , " , T_E_mod: " , T_E_mod

fm = open("Source_xme.dat",'w')
#fe = open("Source_engy.dat",'w')
fm.write("%d" % IX)
fm.write("\n")

#=================
# read probe file
#=================
nx_probe = 200
x_probe = np.array(range(nx_probe))
Lx_probe = 1.85
x_probe = Lx_probe*x_probe/nx_probe

#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe/centerline.T"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_BOX_1D_CHEMTABLE3D/centerline.T"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_BOX_1D_caloricallyPerfect/centerline.T"
filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_GTTL_ChemTable3D/centerline.T"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_GTTL_caloricallyPerfect/centerline.T"
pr = rp.lineProbe(filename)
#T_probe  = pr.readLastLine()
#T_probe  = pr.averageLines()
T_probe  = pr.averageFromLines(10)
print "x_probe: " , x_probe, " , T_probe: " ,T_probe

#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe/centerline.rho"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_BOX_1D_CHEMTABLE3D/centerline.rho"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_BOX_1D_caloricallyPerfect/centerline.rho"
filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_GTTL_ChemTable3D/centerline.rho"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_GTTL_caloricallyPerfect/centerline.rho"
pr = rp.lineProbe(filename)
#rho_probe  = pr.readLastLine()
rho_probe  = pr.averageFromLines(10)

#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe/centerline.U-x"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_BOX_1D_CHEMTABLE3D/centerline.U-x"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_BOX_1D_caloricallyPerfect/centerline.U-x"
filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_GTTL_ChemTable3D/centerline.U-y"
#filename = "/Users/miladmortazavi/CASCADE_CODES/NEXTGEN/nextgen/src/charles/probe_GTTL_caloricallyPerfect/centerline.U-y"
pr = rp.lineProbe(filename)
#U_probe  = pr.readLastLine()
U_probe  = pr.averageFromLines(10)

rhoU_probe = rho_probe*U_probe

for i in range(IX):
  fm.write("%.16e" % x[i])
  fm.write("  ")
  fm.write("%.16e" % Sm[i])
  fm.write("  ")
  fm.write("%.16e" % Se[i])
  fm.write('\n')
fm.close()

plt.figure()
#plt.subplot(1,2,1)
plt.plot(x[0:IX],D[0:IX]**2)
plt.xlabel("x(m)")
plt.ylabel("d^2(m^2)")
#plt.show(block=False)
plt.savefig("d2.png")

plt.figure()
#plt.subplot(1,2,2)
n = 20
plt.plot(x[0:IX-n],Tl[0:IX-n])
plt.xlabel("x(m)")
plt.ylabel("T(K)")
#plt.show(block=False)
plt.savefig("T.png")

plt.figure()
plt.plot(x[0:IX],Mw[0:IX],label="Water mass")
plt.plot(x[0:IX],Mv[0:IX],label="Vapor mass")
plt.plot(x[0:IX],Mw[0:IX]+Mv[0:IX],label="Water+Vapor mass")
plt.legend()
plt.xlabel("x(m)")
plt.ylabel("mass(kg)")
#plt.show(block=False)
plt.savefig("mass.png")

plt.figure()
plt.plot(x[0:IX],ug[0:IX],label='1D model')
plt.plot(x_probe,U_probe,label='Charles with source terms')
plt.xlabel("x(m)")
plt.ylabel("ug(m/s)")
plt.legend()
plt.show(block=False)
plt.savefig("ug.png")

plt.figure()
plt.plot(x[0:IX],rhog[0:IX],label='1D model')
plt.plot(x_probe,rho_probe,label='Charles with source terms')
plt.xlabel("x(m)")
plt.ylabel("rhog(kg/m^3)")
plt.legend()
plt.show(block=False)
plt.savefig("rhog.png")

plt.figure()
plt.plot(x[0:IX],Tg[0:IX],label='1D model')
plt.plot(x_probe,T_probe,label='Charles with source terms')
plt.xlabel("x(m)")
plt.ylabel("Tg(K)")
plt.legend()
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

print "mass flow rate: " , mg[IX-1], " , mdot_g+mdot_w[0]: ", mdot_g+mdot_w[0] , " , x final: ", x[IX], " , Mach final:" , M[IX]
plt.figure()
plt.plot(x[0:IX],mg[0:IX])
plt.xlabel("x(m)")
plt.ylabel("mass flow(kg/s)")
plt.show(block=False)

plt.figure()
plt.plot(x[0:IX],mdot_w[0:IX],label="mdot water")
plt.plot(x[0:IX],mdot_v[0:IX],label="mdot vapor")
mdot_total = mdot_w+mdot_v
plt.plot(x[0:IX],mdot_total[0:IX],label="mdot total")
plt.xlabel("x(m)")
plt.ylabel("mdot(kg/s)")
plt.legend()
plt.show(block=False)
plt.savefig("mdot.png")

plt.figure()
plt.plot(x[0:IX],heat_flow[0:IX],label="heat flow gas")
plt.plot(x[0:IX],Q[0:IX],label="heat transfer")
plt.xlabel("x(m)")
plt.ylabel("heat flow(J/s)")
plt.legend()
plt.show(block=False)

plt.figure()
plt.plot(x[0:IX],Sm[0:IX])
plt.xlabel("x(m)")
plt.ylabel("mass source (kg/s/m^3)")
plt.show(block=False)
plt.savefig("Source_m.png")

plt.figure()
plt.plot(x[1:IX],Se[1:IX])
plt.xlabel("x(m)")
plt.ylabel("energy source (J/s/m^3)")
#plt.show()
plt.show(block=False)
plt.savefig("Source_e.png")

rhoug = rhog*ug
plt.figure()
plt.plot(x[0:IX-1],rhoug[0:IX-1],label='1D model')
plt.plot(x_probe,rhoU_probe,label='Charles with source terms')
plt.xlabel("x(m)")
plt.ylabel("rho*u")
plt.legend()
plt.savefig("RhoTimesU.png")
plt.show()

