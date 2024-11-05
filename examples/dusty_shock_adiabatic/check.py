ed""" @package ./examples/shocktube_1d/check
Code that checks results of 1d shocktube problem

created by Rainer Weinberger, last modified: 19.02.2019
"""

""" load libraries """
import sys    # needed for exit codes
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os      # file specific calls
import matplotlib.pyplot as plt    ## needs to be active for plotting!
from dusty_shock_adiabatic import shock
plt.rcParams['text.usetex'] = True

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("dusty_shock_adiabatic: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
    
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = np.sqrt(NumberOfCells) ## 2d sim


""" loop over all output files """
status = 0
i_file = 5
errors = []
GAMMA = 7./5.

M = 10     #mach number
dust_gas_ratio = 0.5 
shock_length = 1.0
P_0 = 1.0
FB = 1.0
K = 3.0


shock_steps = 1000.
offset=50-0.05

while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break

    
    """ get simulation data """
    
    xc = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    vel = np.array(data["PartType0"]["Velocities"], dtype = FloatType)[:,0]
    vel_dust = np.array(data["PartType0"]["DustVelocities"], dtype=FloatType)[:,0]
    rho = np.array(data["PartType0"]["Density"], dtype = FloatType)
    d_rho = np.array(data["PartType0"]["DustDensity"], dtype=FloatType)
    mass = np.array(data["PartType0"]["Masses"], dtype=FloatType)
    d_mass = np.array(data["PartType0"]["DustMasses"], dtype=FloatType)
    u = np.array(data['PartType0']['InternalEnergy'], dtype = FloatType)
    
    Pressure = u * rho * (GAMMA - 1)
    time = FloatType(data["Header"].attrs["Time"])
    
    total_d_mass = np.sum(d_mass)

    true = shock(M, dust_gas_ratio, {'drag_type':'power_law', 'drag_const':1.0}, shock_length, shock_steps, GAMMA, P_0, t=time, FB=FB, Kin=K, offset=offset)
    

    """ True solution """
    #sol = DustyWaveSolver(K=100.0, delta=1e-6, feedback=1.0)(time)

    x = xc[:,0]# - 1.0*time
    #rho_gas = sol.rho_gas(x)
    
    #rho_dust = sol.rho_dust(x)
    #velocity_dust = sol.v_dust(x)
    #velocity_gas = sol.v_gas(x)

    #fig = plt.figure(figsize=np.array([7.0,7.0]), dpi=300)
    #ax = plt.axes([0.1,0.1,0.8,0.8])
    f, subs = plt.subplots(3, 1, sharex=True)
    subs[0].plot(x, vel, c="r", label="Gas")
    subs[0].plot(x, vel_dust, c="k", label="Dust")
    subs[0].plot(true["xi"], true["wd"], c="gray", ls="--", label="True Dust")
    subs[0].plot(true["xi"], true["wg"], c="pink", ls="--", label="True Gas" )
    subs[0].set_xlim([0, 60.0])

    subs[0].set_ylabel("Velocity")
    subs[0].legend(loc="best")
    f.suptitle("Adiabatic DustyShock, Mach 10    High Resolution (N=500)    Low Resolution (N=20)")
    
    #plt.plot(x, velocity_dust, label="Analytical")
    subs[1].plot(x, rho, c="r", label="Gas")
    subs[1].plot(x, d_rho, c="k", label="Dust")
       
    subs[1].plot(true["xi"], true["rhog"], c="pink", ls="--", label="True Gas")
    subs[1].plot(true["xi"], true["rhod"], c="gray", ls="--", label="True Dust")
    subs[1].set_ylabel("Density")
    subs[1].set_ylim([0, 6.0])
    subs[1].set_xlim([0, 60.0])
        
    subs[2].plot(x, Pressure, c="red", label="Pressure")
    subs[2].plot(true['xi'], true['P'], c="pink", ls="--", label="True Pressure")
    subs[2].set_ylabel("Pressure")
    subs[2].set_xlabel("position")
    subs[2].set_xlim([0, 60.0])
    subs[2].set_ylim([0, 50])
    f.savefig(simulation_directory+"/figure_all_%03d.png" % (i_file), dpi=300)
    plt.show()
    plt.close()
    i_file += 1



