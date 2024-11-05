""" @package examples/dusty_shock_adiabatic/create.py
Code that creates 1d adiabatic dusty shock test.

created by Maggie Celeste, 2024
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/dusty_shock_adiabatic/create.py: creating ICs in directory " + simulation_directory)

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(50.0)
NumberOfCells = IntType(20)

FilePath = simulation_directory + '/IC.hdf5'

""" initial condition parameters """
M = 10     #mach number
dust_gas_ratio = 0.5 
shock_length = 1.0
rho_g = 1.0
P_0 = 1.0
FB = 1.0

rho_d = rho_g * dust_gas_ratio
GAMMA = FloatType(7./5.0)
gamma_minus_one = FloatType(GAMMA - 1.0)
c_s = np.sqrt(GAMMA * P_0 / rho_g)
v_s = M * c_s

A = (1+FB*dust_gas_ratio)* (GAMMA+1) / (2*GAMMA)
B = - ( v_s*(1+FB*dust_gas_ratio) + P_0 / (rho_g*v_s))
C = (GAMMA-1)/(2*GAMMA) *(v_s**2 * (1+FB*dust_gas_ratio) ) + P_0/rho_g

v_post = (-B - np.sqrt(B**2 - 4*A*C)) / (2*A)

dv = v_s - v_post


""" set up grid """
## spacing
dx = Boxsize / FloatType(NumberOfCells)
## position of first and last cell
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx
## set up grid
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = np.linspace(pos_first, pos_last, NumberOfCells, dtype=FloatType)
Volume = np.full(NumberOfCells, dx, dtype=FloatType)

## left state
Mass = np.full(NumberOfCells, rho_g*dx, dtype=FloatType)
DustMass = np.full(NumberOfCells, rho_d*dx, dtype=FloatType)
Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
Velocity[:,0] = dv
Uthermal = np.full(NumberOfCells, (P_0/rho_g/gamma_minus_one), dtype=FloatType)





""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells

## write header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)
part0.create_dataset("DustVelocities", data=Velocity)
part0.create_dataset("DustMasses", data = DustMass)
## close file
IC.close()

""" normal exit """
sys.exit(0) 
