import numpy as np
from scipy.special import sph_harm
from matplotlib import pyplot as plt
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def generate_field_on_sphere(radius, characteristic_scale, num_points=100):
    """
    Generate a field on a sphere with a given characteristic length scale.
    
    Parameters:
    radius (float): Radius of the sphere.
    characteristic_scale (float): Characteristic length scale on the sphere.
    num_points (int): Number of points to discretize the sphere (for theta and phi).
    
    Returns:
    theta (2D array): Array of theta values.
    phi (2D array): Array of phi values.
    field (2D array): Real part of the generated field on the sphere.
    """
    
    # Calculate the corresponding l for the characteristic scale
    # Assuming characteristic_scale is given in radians and is the arc length
    l_max = int(np.round(np.pi / characteristic_scale) + 10)
    l_min = int(np.round(np.pi / characteristic_scale) - 10)  # Narrow range around l_max
    
    # Generate grid for theta and phi
    theta, phi = np.linspace(0, np.pi, num_points), np.linspace(0, 2 * np.pi, num_points)
    theta, phi = np.meshgrid(theta, phi)
    
    # Initialize field
    field = np.zeros_like(theta, dtype=complex)
    
    # Generate field using spherical harmonics
    for l in tqdm(range(l_min, l_max + 1), desc='Generating field, l starts ' + str(l_min) + ':'):
        for m in tqdm(range(-l, l+1)):
            # Random coefficients
            coefficient = np.random.normal(size=2).view(complex)[0]
            
            # Add the spherical harmonic contribution to the field
            field += coefficient * sph_harm(m, l, phi, theta)
            if np.isnan(field).any():
                print('NaN encountered at l = ' + str(l) + ', m = ' + str(m))
                print('coefficient = ' + str(coefficient))
                print('sph_harm(m, l, phi, theta) = ' + str(sph_harm(m, l, phi, theta)))
                print('field = ' + str(field))
                return theta, phi, field.real
    
    # Return the real part of the field
    return theta, phi, field.real


def plot_field_on_sphere(theta, phi, field, radius=1):
    """
    Plot the field (theta, phi, field) on a sphere.
    
    Parameters:
    theta (2D array): Array of theta values.
    phi (2D array): Array of phi values.
    field (2D array): Real part of the generated field on the sphere.
    radius (float): Radius of the sphere.
    """
    # Convert spherical coordinates to Cartesian coordinates for plotting
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    
    # Normalize the field for better color representation
    field_norm = (field - field.min()) / (field.max() - field.min())
    
    # Plot the sphere with the field
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.viridis(field_norm),
                    linewidth=0, antialiased=False)
    
    # Hide the axes
    ax.set_axis_off()
    
    # Show the plot
    plt.show()

# Example usage:
# For a sphere of radius 1 and a characteristic scale of 0.1 radians
# theta, phi, field = generate_field_on_sphere(1, 0.1, 1000)
# np.save('field_2.npy', field)

# for rerunning
field = np.load('field.npy')
field_2 = np.load('field_2.npy')
field = field + field_2
num_points = 1000
theta, phi = np.linspace(0, np.pi, num_points), np.linspace(0, 2 * np.pi, num_points)
theta, phi = np.meshgrid(theta, phi)
# plt.contourf(theta, phi, field)
plt.hist(field.flatten(), bins=100)
print(np.mean(field.flatten()*field.flatten()))
plt.show()
# plot_field_on_sphere(theta, phi, field, 1)

# Output to text file
# with open('field.dat', 'wb') as out:
#     for i in range(field.shape[0]):
#         out.write(bytes(' '.join(str(num) for num in field[i,:]),'utf-8'))
#         out.write(bytes('\r\n','utf-8'))
