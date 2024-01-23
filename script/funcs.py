import subprocess
import pkg_resources
required_packages = ["numpy", "scipy", "matplotlib"] 
installed_packages = {pkg.key for pkg in pkg_resources.working_set}
missing_packages = [pkg for pkg in required_packages if pkg not in installed_packages]
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if missing_packages:
    subprocess.check_call(["python", "-m", "pip", "install", *missing_packages])

def four_parameter_logistic(x, ec50, hill_slope, top, bottom):
    # return bottom + (top - bottom) / (1 + (x / ec50)**(-hill_slope))
    return bottom + (top - bottom) / (1 + 10**((np.log10(ec50) - np.log10(dose)) * hill_slope))

def fit_4pl_model_py(x_data, y_data, initial_guess=(50, 1, 100, 0)):
    popt, _ = curve_fit(four_parameter_logistic, x_data, y_data, p0=initial_guess, maxfev=10000)
    return popt
  
def calculate_fitted_response(dose, popt):
    ec50, hill_slope, top, bottom = popt
    fitted_response = four_parameter_logistic(dose, ec50, hill_slope, top, bottom)
    return fitted_response

def plot_dose_response_with_fit_py(dose, response, fitted_response, title="Dose-Response Curve", color = "blue"):
    plt.figure(figsize=(8, 6))
    plt.scatter(dose, response, color=color, label='Observed Data')  # Plot observed data points
    plt.plot(dose, fitted_response, color=color, label='Fitted Curve')  # Plot fitted curve
    plt.title(title)
    plt.xlabel("Dose")
    plt.ylabel("Response")
    plt.legend()
    plt.grid(True)
    # plt.show()
    
def simulate_dose_response(
  n=100, response_type='inhibition', ec50=50, hill_coeff=1, 
  max_response=100, min_response=0, noise_level=0.05):
    dose = np.linspace(min_response, max_response, n)

    if response_type == 'inhibition':
        response = max_response / (1 + (dose/ec50)**hill_coeff)
    elif response_type == 'stimulation':
        response = min_response + (max_response - min_response) * (dose**hill_coeff / (ec50**hill_coeff + dose**hill_coeff))
    else:
        raise ValueError("response_type must be 'inhibition' or 'stimulation'")

    # Add random noise to simulate experimental variability
    response += np.random.normal(0, noise_level * max_response, n)

    return dose, response
  
  
