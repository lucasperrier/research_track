import matplotlib.pyplot as plt
import numpy as np


def plot_time(Y,t, parameters_list, initial_conditions):
	plt.plot(t, Y[0], label=f'x(t)')
	plt.plot(t, Y[1], label=f'y(t)')
	plt.xlabel('Time t [seconds]')
	plt.ylabel('Population')
	plt.suptitle(f'alpha={parameters_list[0]}, beta={parameters_list[1]}, gamma={parameters_list[2]}, delta={parameters_list[3]}', y=0.99, fontsize=10)
	plt.title('Time Evolution of x and y')
	plt.legend()
	plt.grid(True)
	plt.show()

def plot_x_vs_y(Y,t, parameters_list, initial_conditions):
	show_plot=False
	plt.plot(Y[0], Y[1], label=f'y0={initial_conditions[1]:.2f}')
	plt.xlabel('x(t)')
	plt.ylabel('y(t)')
	plt.suptitle(f'alpha={parameters_list[0]}, beta={parameters_list[1]}, gamma={parameters_list[2]}, delta={parameters_list[3]}, x0={initial_conditions[0]}', y=0.99, fontsize=10)
	plt.title('y(t) vs x(t) with t as a parameter')
	plt.legend()
	plt.grid(True)
    # Only show plot when all trajectories have been plotted
	if show_plot:
		plt.show()

def plot_equilibirum(Y,t, parameters_list, initial_conditions):
	fig, axs = plt.subplots(1, 2, figsize=(14, 5))

	# Plot Time Evolution of x(t) and y(t)
	axs[0].plot(t, Y[0], label=f'x(t)')
	axs[0].plot(t, Y[1], label=f'y(t)')
	axs[0].set_xlabel('Time t [seconds]')
	axs[0].set_ylabel('Population')
	axs[0].set_title('Time Evolution of x(t) and y(t)')
	axs[0].legend()
	axs[0].grid(True)

	# Plot y(t) vs x(t) for different initial conditions of y
	axs[1].plot(Y[0], Y[1], label=f'y(x(t))')
	axs[1].set_xlabel('x(t)')
	axs[1].set_ylabel('y(t)')
	axs[1].set_title('y(t) vs x(t)')
	axs[1].legend()
	axs[1].grid(True)

	# Add a shared supertitle for the figure
	fig.suptitle(f'(x*,y*)=({parameters_list[2]/parameters_list[3]:.2f},{parameters_list[0]/parameters_list[1]:.2f}), alpha={parameters_list[0]}, beta={parameters_list[1]}, gamma={parameters_list[2]}, delta={parameters_list[3]}', fontsize=12)

	# Show the plots
	plt.tight_layout()
	plt.show()


def plot_side_by_side(Y, t, parameters_list, initial_conditions):
    # Create a figure with two subplots side by side
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))

    # Plot Time Evolution of x(t) and y(t)
    axs[0].plot(t, Y[0], label=f'x(t)')
    axs[0].plot(t, Y[1], label=f'y(t)')
    axs[0].set_xlabel('Time t [seconds]')
    axs[0].set_ylabel('Population')
    axs[0].set_title('Time Evolution of x(t) and y(t)')
    axs[0].legend()
    axs[0].grid(True)

    # Plot y(t) vs x(t) for different initial conditions of y
    axs[1].plot(Y[0], Y[1], label=f'y(x(t))')
    axs[1].set_xlabel('x(t)')
    axs[1].set_ylabel('y(t)')
    axs[1].set_title('y(t) vs x(t)')
    axs[1].legend()
    axs[1].grid(True)

    # Add a shared supertitle for the figure
    fig.suptitle(f'alpha={parameters_list[0]}, beta={parameters_list[1]}, gamma={parameters_list[2]}, delta={parameters_list[3]}', fontsize=12)

    # Show the plots
    plt.tight_layout()
    plt.show()

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