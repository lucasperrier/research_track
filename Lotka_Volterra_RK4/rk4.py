import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from plotting import *


alpha = 1.1
beta =0.4
gamma = 0.1
delta = 0.4

def lotka_volterra_system(t,y, parameters_list):
	alpha, beta, gamma, delta = parameters_list
	x,y = y
	'''
	computes dx/dt and dy/dt in lotka volterra model
	IN -
	y: x and y value, np.array, shape = (1,2)
	parameters_list: list of values taken up by model [alpha, beta, gamma, delta]
	OUT - 
	value of dy/dt, dx/dt, np.array, shape = (1,2)
	'''
	dxdt = alpha * x - beta * x * y
	dydt = delta * x * y - gamma * y
	return np.array([dxdt, dydt])




def rk4_step(model, dt, t0, y0, parameters_list):
	'''
	performs one runge kutta 4 step
	IN -
	model: model in vectorized format that you want to integrate, python function
	dt: timestep, float
	t0: value of t at current step, float, shape = (1)
	y0: value of y at current step ( x(t),y(t) ), float, shape = (1,2)
	parameters_list: list of values taken up by model [alpha, beta, gamma, delta]
	OUT - 
	yout: value of y at next step ( x(t+1), y(t+1) ), float, shape = (1,2)
	'''
	f1 = model(t0,y0, parameters_list)
	f2 = model(t0+(dt/2), y0+(dt/2)*f1, parameters_list)
	f3 = model(t0+(dt/2), y0+(dt/2)*f2, parameters_list)
	f4 = model(t0+dt, y0+dt*f3, parameters_list)
	yout = y0 + (dt/6) * (f1+ 2*f2 + 2*f3 +f4)
	return yout

def simulation(method, system, dt, boundary, initial_conditions, plotting_function, parameters_list):
	'''
	performs one simulation of lotka volterra model using runge kutta 4 method
	IN:
	method: numerical method used (example: rk4_step)
	system: model you want to integrate (example: lotka_volterra_system)
	dt: size step size h to use in simulation, float
	boundary: time boundary simulation runs until
	initial_conditions: initial values of x and y
	plotting_function: function used to plot simulation
	parameters_list: parameters used in simulation
	OUT:
	Y: solution matrix containting time evolution of x and y,  shape: (2, number of time points))
	t: time points vector, shape: (1, number of time points)
	'''
	# parameter values
	alpha = parameters_list[0]
	beta = parameters_list[1]
	gamma = parameters_list[2]
	delta = parameters_list[3]

	T = boundary
	# initial conditions 
	y0 = initial_conditions
	# discrete time vector
	number_of_timepoints = int(T/dt)
	# Generate the time array with N+1 points to include the initial time
	t = np.linspace(0, T, number_of_timepoints+1)
	# initialize solution matrix (time evolution of x and y)
	Y = np.zeros((2,number_of_timepoints+1))
	Y[:,0] = y0
	# use initial conditions as starting point in loop
	yin = y0

	# simulation
	for i in tqdm(range(number_of_timepoints), desc="Simulation Progress"):
		yout = method(system,dt,t[i],yin, parameters_list)
		Y[:, i+1] = yout
		yin = yout
	# plot solution found
	plotting_function(Y,t,parameters_list, initial_conditions)
	return Y, t





