from rk4 import *
from convergence_test import *
from generate_data import *
from plotting import *

# parameter values
alpha = 1.1
beta = 0.4
gamma = 0.4
delta = 0.1
parameters_list = [alpha,beta,gamma,delta]
#initial conditions
x_0 = 10
y_0 = 10
y0 = [x_0, y_0]
# boundary definition
T = 100

# Define suitable ranges for parameters
alpha_range = np.linspace(1.1, 2.1, 1)  # Adjust the number of points as needed
beta_range = np.linspace(0.4, 1.4, 1)
gamma_range = np.linspace(0.1, 1.1, 1)
delta_range = np.linspace(0.4, 1.4, 4)

# Define suitable ranges for initial conditions
x0_range = np.linspace(10, 20, 1)
y0_range = np.linspace(1.0, 30, 8)

# Generate parameter combinations
parameter_combinations_list = generate_parameter_combinations(alpha_range, beta_range, gamma_range, delta_range)

# Generate initial conditions
initial_conditions_list = generate_initial_conditions(x0_range, y0_range)

#equilibirum initial conditions
equilibrium_x_0 = parameters_list[2]/parameters_list[3]
equilibrium_y_0 = parameters_list[0]/parameters_list[1]
equilibrium_y0 = [equilibrium_x_0, equilibrium_y_0]

# step size values for convergence test
h_values = generate_h_values(8)
def main():
	#generate_time_data(parameter_combinations_list, initial_conditions_list, time_boundary=50, dt=0.0001, output_folder='simulation_data')
	'''
	simulation(method=rk4_step, 
					system=lotka_volterra_system, 
					dt=0.001, 
					boundary=T, 
					initial_conditions=equilibrium_y0,
					plotting_function = plot_equilibirum,
					parameters_list=parameters_list)
	
	for initial_condition in initial_conditions_list:
		Y_y0, t_h = simulation(method=rk4_step, 
					system=lotka_volterra_system, 
					dt=0.001, 
					boundary=T, 
					initial_conditions=initial_condition,
					plotting_function = plot_x_vs_y,
					parameters_list=parameters_list)
	plt.show()
	for parameter_combination in parameter_combinations_list:
		Y_h, t_h = simulation(method=rk4_step, 
					system=lotka_volterra_system, 
					dt=0.001, 
					boundary=T, 
					initial_conditions=y0,
					plotting_function = plot_side_by_side,
					parameters_list=parameter_combination)
	'''
	error_vs_h, prey_error_vs_h, predator_error_vs_h = convergence_test(h_values=h_values, 
				  method=rk4_step, 
				  system=lotka_volterra_system, 
				  boundary=T, 
				  initial_conditions=y0, 
				  plotting_function=plot_time, 
				  parameters_list=parameters_list)
	
	plot_convergence_test(error_vs_h)
	plot_convergence_test(prey_error_vs_h)
	plot_convergence_test(predator_error_vs_h)
	
	return None
	
if __name__ == '__main__':
	main()