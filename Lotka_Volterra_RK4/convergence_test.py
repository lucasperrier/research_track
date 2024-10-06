from rk4 import * 
from scipy.interpolate import interp1d


def convergence_test(h_values, method, system, boundary, initial_conditions, plotting_function, parameters_list):
	'''
	performs convergence test for a runge kutta 4 method. 
	Plots log(error(h)) vs h, and a line of best fit and its slope.
	IN:
	h_values: list of step size values ussed in convergence test
	method: numerical method tested
	system: model numerical method solves
	boundary: time boundary simulations run until
	initial_conditions: initial conditions used for all simulations
	plotting_function: plotting function for log log plot
	parameters_list: parameters used in all simulations.
	OUT:
	error_vs_h: array of total error vs h data, shape (2, number of simulations)
	prey_error_vs_h: array of error for prey curve vs h data, shape (2, number of simulations)
	predator_error_vs_h: array of error for predator curve vs h data, shape (2, number of simulations)
	'''
	error_vs_h = np.zeros((2,len(h_values)))
	prey_error_vs_h = np.zeros((2,len(h_values)))
	predator_error_vs_h = np.zeros((2,len(h_values)))

	Y_ref, t_ref = simulation(method, 
			system, 
			dt=0.0001, 
			boundary=boundary, 
			initial_conditions=initial_conditions,
			plotting_function = plotting_function,
			parameters_list=parameters_list)

	prey_ref = Y_ref[0]
	predator_ref = Y_ref[1]

	for i,h in enumerate(h_values):
		Y_h, t_h = simulation(method, 
			system, 
			dt=h, 
			boundary=boundary, 
			initial_conditions=initial_conditions,
			plotting_function = plotting_function,
			parameters_list=parameters_list)
		t_ref[-1] = t_h[-1]

		prey_h_data = Y_h[0]
		predator_h_data = Y_h[1]

		interpolator_prey = interp1d(t_h, prey_h_data, kind='cubic', bounds_error=False, fill_value='extrapolate')
		prey_h_data_interpolated = interpolator_prey(t_ref)
		error_prey = np.linalg.norm(prey_ref - prey_h_data_interpolated, ord=2) 
		interpolator_predator = interp1d(t_h, predator_h_data, kind='cubic', bounds_error=False, fill_value='extrapolate')
		predator_h_data_interpolated = interpolator_predator(t_ref)

		error_predator = np.linalg.norm(predator_ref - predator_h_data_interpolated, ord=2) 

		total_error_h = error_prey + error_predator

		prey_error_vs_h[0,i] = h
		predator_error_vs_h[0,i] = h
		error_vs_h[0,i] = h

		prey_error_vs_h[1,i] = error_prey
		predator_error_vs_h[1,i] = error_predator
		error_vs_h[1,i] = total_error_h
	return error_vs_h, prey_error_vs_h, predator_error_vs_h

def generate_h_values(n):
	h_values = [.1]
	h_0 = h_values[0]
	for i in range(n):
		h = h_0 / 2
		h_values.append(h)
		h_0 = h
	return h_values


def plot_convergence_test(error_vs_h):
	h_values = error_vs_h[0, :]
	error_values = error_vs_h[1, :]
	# Filter out zero or negative errors
	valid_indices = error_values > 0
	h_values = h_values[valid_indices]
	error_values = error_values[valid_indices]
	log_h = np.log(h_values)
	log_error = np.log(error_values)
	slope, intercept = np.polyfit(log_h, log_error, 1)  # 1 means a linear fit (degree 1)
	y_fit = slope * log_h + intercept
	plt.plot(log_h, log_error, marker='o', label='results')
	plt.plot(log_h, y_fit, label=f'Best Fit Line (slope={slope:.2f})', color='red')
	plt.xlabel('Step size h (log)')
	plt.ylabel('Error E(h) (log)')
	plt.title('Error vs Step Size')
	plt.legend()
	plt.text(0.05, 0.95, f'Order of convergence = {slope:.2f}', transform=plt.gca().transAxes,
         fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
	plt.grid(True)
	plt.show()