import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Load data from text file
data = np.loadtxt('path_to_your_file.txt')  # Load data into numpy array, adjust as necessary

# Separate data into x and y arrays
x_data = data[:, 0]  # Assuming the first column contains x values
y_data = data[:, 1]  # Assuming the second column contains y values

# Define a model function
#def model_function(x, *params):
    # Define the model equation here, for example, a simple polynomial
    #return params[0] * x**2 + params[1] * x + params[2]

# Initial guess for the parameters
#initial_guess = [1.0, 1.0, 1.0]  # Adjust as per your model

# Fit the model to the data
#params, covariance = curve_fit(model_function, x_data, y_data, p0=initial_guess)

# Extracting fit parameters
#fit_params = params

# Generate fitted curve
#fitted_curve = model_function(x_data, *fit_params)

# Plot original data and fitted curve
plt.scatter(x_data, y_data, label='Original Data')
#plt.plot(x_data, fitted_curve, color='red', label='Fitted Curve')
plt.xlabel('x')
plt.ylabel('y')
#plt.title('Curve Fitting')
plt.legend()
plt.show()

#print("Fit Parameters:", fit_params)import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Load data from text file
data = np.loadtxt('W:\Users\Marta\gvg.txt')  # Load data into numpy array, adjust as necessary

# Separate data into x and y arrays
x_data = data[:, 0]  # Assuming the first column contains x values
y_data = data[:, 1]  # Assuming the second column contains y values

# Define a model function
#def model_function(x, *params):
    # Define the model equation here, for example, a simple polynomial
    #return params[0] * x**2 + params[1] * x + params[2]

# Initial guess for the parameters
#initial_guess = [1.0, 1.0, 1.0]  # Adjust as per your model

# Fit the model to the data
#params, covariance = curve_fit(model_function, x_data, y_data, p0=initial_guess)

# Extracting fit parameters
#fit_params = params

# Generate fitted curve
#fitted_curve = model_function(x_data, *fit_params)

# Plot original data and fitted curve
plt.scatter(x_data, y_data, label='Original Data')
#plt.plot(x_data, fitted_curve, color='red', label='Fitted Curve')
plt.xlabel('x')
plt.ylabel('y')
#plt.title('Curve Fitting')
plt.legend()
plt.show()

#print("Fit Parameters:", fit_params)