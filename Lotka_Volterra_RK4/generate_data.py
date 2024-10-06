from rk4 import *
import numpy as np
import itertools
import pandas as pd
from tqdm import tqdm

import os

def generate_time_data(parameter_combinations_list, initial_conditions_list, time_boundary=50, dt=0.0001, output_folder='simulation_data'):
    '''
    Generate time series data for combinations of parameters and initial conditions,
    and save the data to separate CSV files per simulation.

    IN:
    - parameter_combinations_list: List of parameter lists [alpha, beta, gamma, delta].
    - initial_conditions_list: List of initial condition lists [x0, y0].
    - time_boundary: Time boundary for the simulation.
    - dt: Time step size.
    - output_folder: Folder to save the simulation data.

    The function saves a summary CSV file and individual CSV files for time series data.
    '''
    records = []  # List to hold summary data

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Generate all combinations of parameters and initial conditions
    all_combinations = list(itertools.product(parameter_combinations_list, initial_conditions_list))
    for idx, (parameters, initial_conditions) in tqdm(enumerate(all_combinations), total=len(all_combinations)):
        print(f"passed parameters, initial:  {parameters}, and {initial_conditions}")
        # Run the simulation
        Y, t = simulation(
            method=rk4_step,
            system=lotka_volterra_system,
            dt=dt,
            boundary=time_boundary,
            initial_conditions=initial_conditions,
            plotting_function=None,
            parameters_list=parameters
        )

        # Extract end values
        x_end = Y[0, -1]
        y_end = Y[1, -1]

        # Prepare the record dictionary
        simulation_id = f'simulation_{idx}'
        record = {
            'simulation_id': simulation_id,
            'alpha': parameters[0],
            'beta': parameters[1],
            'gamma': parameters[2],
            'delta': parameters[3],
            'x0': initial_conditions[0],
            'y0': initial_conditions[1],
            'time_boundary': time_boundary,
            'dt': dt,
            'x_end': x_end,
            'y_end': y_end,
            'timeseries_file': f'{simulation_id}.csv'
        }

        records.append(record)

        # Save time series data to a CSV file
        timeseries_df = pd.DataFrame({
            't': t,
            'x_t': Y[0, :],
            'y_t': Y[1, :]
        })
        timeseries_df.to_csv(os.path.join(output_folder, record['timeseries_file']), index=False)

    # Create a summary DataFrame and save to CSV
    summary_df = pd.DataFrame(records)
    summary_df.to_csv(os.path.join(output_folder, 'simulation_summary.csv'), index=False)

    print(f"Data has been saved to '{output_folder}' folder.")


def generate_parameter_combinations(alpha_range, beta_range, gamma_range, delta_range):
    '''
    Generates a list of lists containing all combinations of parameters alpha, beta, gamma, delta.
    IN:
    - alpha_range: List or array of alpha values
    - beta_range: List or array of beta values
    - gamma_range: List or array of gamma values
    - delta_range: List or array of delta values
    OUT:
    - parameter_combinations_list: List of parameter lists [alpha, beta, gamma, delta]
    '''
    # Generate all combinations of parameters
    parameter_values = itertools.product(alpha_range, beta_range, gamma_range, delta_range)
    parameter_combinations_list = []
    for alpha, beta, gamma, delta in parameter_values:
        params = [alpha, beta, gamma, delta]
        parameter_combinations_list.append(params)
    return parameter_combinations_list



def generate_initial_conditions(x0_range, y0_range):
    '''
    Generates a list of lists containing all combinations of initial conditions x0, y0.
    IN:
    - x0_range: List or array of x0 values
    - y0_range: List or array of y0 values
    OUT:
    - initial_conditions_list: List of initial condition lists [x0, y0]
    '''
    # Generate all combinations of initial conditions
    initial_values = itertools.product(x0_range, y0_range)
    initial_conditions_list = []
    for x0, y0 in initial_values:
        initial_conditions_list.append([x0, y0])
    return initial_conditions_list

