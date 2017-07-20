import math as m
import numpy as np
import matplotlib.pyplot as plt

class lineProbe(object):

  def __init__(self,filename):
    self.filename = filename

  def readLastLine(self):
    fh = open(self.filename,'r')
    lines = fh.readlines()
    line = lines[-1].strip()
    col = line.split()
    nstep = int(col[0])
    time = float(col[1])
    nx = int(col[2])
    print "nstep: " , nstep, " , time: " , time , " , nx: " , nx
    probe = np.array([float(i) for i in col[3:nx+3]])
    fh.close()
    return probe

  def averageLines(self):
    fh = open(self.filename,'r')
    lines = fh.readlines()
    line = lines[0].strip()
    col = line.split()
    nstep0 = int(col[0])
    time0 = float(col[1])

    line = lines[-1].strip()
    col = line.split()
    nstep1 = int(col[0])
    time1 = float(col[1])
    nx = int(col[2])

    probe = np.zeros(nx)

    for line in lines:
      line = line.strip()
      col = line.split()
      probe += np.array([float(i) for i in col[3:nx+3]])

    probe /= float(len(lines))
    fh.close()
    print "averaging nsteps " , nstep0 , ":" , nstep1, " , time " , time0 , ":" , time1
    return probe

  def averageFromLines(self,nl0,nl1=None):
    fh = open(self.filename,'r')
    lines = fh.readlines()
    line = lines[nl0].strip()
    col = line.split()
    nstep0 = int(col[0])
    time0 = float(col[1])

    if (nl1 == None):
      nl1 = len(lines) - 1
    line = lines[nl1].strip()
    col = line.split()
    nstep1 = int(col[0])
    time1 = float(col[1])
    nx = int(col[2])

    probe = np.zeros(nx)

    for line in lines[nl0:nl1+1]:
      line = line.strip()
      col = line.split()
      probe += np.array([float(i) for i in col[3:nx+3]])

    probe /= float(len(lines[nl0:nl1+1]))
    fh.close()
    print "averaging lines ", nl0 , ":" , nl1, " , length of averaging: ", len(lines[nl0:nl1+1]) , " , nsteps " , nstep0 , ":" , nstep1, " , time " , time0 , ":" , time1
    return probe

