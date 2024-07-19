from KDEpy import FFTKDE, TreeKDE
import numpy as np
from sklearn.model_selection import cross_val_score
from sklearn.base import BaseEstimator
import pdb

class KDEpyWrapper(BaseEstimator):
    def __init__(self, bandwidth):
        self.bandwidth = bandwidth
        self.kde = None
    
    def fit(self, X, y=None, sample_weight=None):
        self.kde = TreeKDE(kernel='gaussian', bw=self.bandwidth, norm=2)
        self.kde.fit(X, weights=sample_weight)
        return self
    
    def score(self, X, y=None):
        # Compute the density estimates for the input data
        density_estimates = self.kde.evaluate(X, eps=0.001)
        # Compute the log likelihood of the density estimates
        log_likelihood = np.log(density_estimates).sum()
        # Return the mean log likelihood
        return log_likelihood / len(density_estimates)
    
    def get_params(self, deep=True):
        return {"bandwidth": self.bandwidth}

    def set_params(self, **params):
        if "bandwidth" in params:
            self.bandwidth = params["bandwidth"]
        return self


# Create 2D data of shape (obs, dims)
data = np.random.randn(2**4, 2)
weights = data[:, 0] ** 2

xmin, xmax = 0, 9
ymin, ymax = 0, 7
data[:, 0] = np.abs(data[:, 0]) * (xmax - xmin) / (2 * np.max(np.abs(data[:, 0]))) + (xmax + xmin) / 2
data[:, 1] = np.abs(data[:, 1]) * (ymax - ymin) / (2 * np.max(np.abs(data[:, 1]))) + (ymax + ymin) / 2

# New data point
new_data_point = np.array([10, 2])

# Append the new data point to the array
data = np.append(data, [new_data_point], axis=0)
weights = data[:, 0] ** 2

print(data)
print(weights)

# Define the bandwidth grid
bandwidths = np.logspace(-1, 1, 19)

print(bandwidths)

#obtain the best bandwidth value to calculate the KDE
# Initialize variables to store results
#best_bandwidth = None
#best_score = float('-inf')  # Update this based on your evaluation metric
#
## Perform grid search
#for bandwidth in bandwidths:
#    # Create KDE instance with the current bandwidth
#    kde = KDEpyWrapper(bandwidth)
#    
#    # Compute scores using cross-validation
#    scores = cross_val_score(kde, data, cv=5, fit_params={'sample_weight': weights})
#    
#    # Compute mean score (you can use other aggregation functions)
#    mean_score = np.mean(scores)
#    
#    # Update best bandwidth if current bandwidth yields better performance
#    if mean_score > best_score:
#        best_score = mean_score
#        best_bandwidth = bandwidth
#
#print("Best bandwidth:", best_bandwidth)

# Assuming you want to define a grid with x axis between 0 and n, and y axis between 0 and m
n = 10  # Define the range for x axis
m = 8   # Define the range for y axis

# Calculate the number of grid points needed based on the range
num_points_x = n + 1
num_points_y = m + 1

# Generate equidistant grid points within the specified ranges
x_grid = np.arange(0, num_points_x)  # Generate grid points for x axis
y_grid = np.arange(0, num_points_y)  # Generate grid points for y axis

# Create a meshgrid from the x and y grid points
X, Y = np.meshgrid(x_grid, y_grid)

# Flatten the meshgrid to obtain grid points in the format required by FFTKDE
#grid_points = np.column_stack([X.ravel(), Y.ravel()])
#print(grid_points)

xmin, xmax, ymin, ymax = 0., n, 0., m
xx, yy = np.mgrid[xmin:xmax:num_points_x*1j, ymin:ymax:num_points_y*1j]
positions = np.vstack([xx.ravel(), yy.ravel()]).T

print(positions)

best_bandwidth = 1.0
# Compute the kernel density estimate
kde = FFTKDE(kernel='gaussian', bw = best_bandwidth, norm=2)
y = kde.fit(data, weights=weights).evaluate(positions)
print(y)
pdb.set_trace()
#print(len(positions))
#print(len(y))

#print(y)

#a = np.array([[1,4], [2,3]])
#b = np.array([[0,1], [0,2],[2, 3], [0,3], [1,4], [5,7]])
#
#print(np.all(a[:, np.newaxis, :] == b, axis=2))
#print(np.where(np.all(a[:, np.newaxis, :] == b, axis=2)))
#
#data_indices_all = np.where(np.all(a[:, np.newaxis, :] == b, axis=2))[1]
#
#print(data_indices_all)