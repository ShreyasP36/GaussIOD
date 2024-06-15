# Gauss Method for Preliminary Orbit Determination
This repository contains MATLAB code for implementing the Gauss method to determine the preliminary orbit of a celestial object based on tracking station observations. 
The provided script includes the following features:
- **Tracking Station Position Calculation**: Determines the inertial position vector of the tracking station at multiple observation times considering Earth's flattening factor.
- **Directional Cosine Vectors**: Computes unit vectors in the direction of the slant range vector.
- **Scalar Quantities Computation**: Calculates cross products and dot products essential for further orbit determination steps.
- **Polynomial Coefficients and Roots Calculation**: Forms and solves an eighth-order polynomial to find the range at the second observation.
- **Initial Position and Velocity Vectors**: Derives the position and velocity vectors at the second observation time.
- **Iterative Refinement**: Applies an iterative algorithm for improved accuracy based on a tolerance level.
- **Orbital Elements Calculation**: Computes specific angular momentum, inclination, ascending node, eccentricity, argument of perigee, true anomaly, and semi-major axis.
- **Orbit Plotting**: Visualizes the calculated orbit in a 3D space along with the equatorial plane and Earth's position.

## Key Concepts
- **Gauss Method**:  An analytical technique for determining the orbit of an object from three observations, following *Algorithm 5.5* from [Orbital Mechanics for Engineering Students, 2nd Edition](https://www.sciencedirect.com/book/9780080977478/orbital-mechanics-for-engineering-students)
- **Lagrange Coefficients**: Used to relate the position and velocity vectors at different times.
- **Newton-Raphson Method**: An iterative numerical technique used for finding successively better approximations to the roots (or zeroes) of a real-valued function.
- **Universal Kepler Equation**: A generalized form of Kepler's equation that applies to all types of orbits (elliptical, parabolic, hyperbolic) and is used to compute the position of an orbiting body over time.
- **Stumpff Functions**: Special functions used in solving Kepler's equation for orbital motion.
- **Iterative Improvement**: Enhances the accuracy of the determined orbit by iterating on the calculated ranges, following *Algorithm 5.6* from [Orbital Mechanics for Engineering Students, 2nd Edition](https://www.sciencedirect.com/book/9780080977478/orbital-mechanics-for-engineering-students)
- **Orbital Elements Calculation**: Extracts orbital elements from the state vector using *Algorithm 4.2* from [Orbital Mechanics for Engineering Students, 2nd Edition](https://www.sciencedirect.com/book/9780080977478/orbital-mechanics-for-engineering-students)

## Getting Started
Clone the repository and run the MATLAB script to execute the orbit determination process. Ensure MATLAB is installed on your system.
```
git clone https://github.com/ShreyasP36/GaussIOD.git
cd GaussIOD
```

## Required 
To use this script, you need the following data:
- *Right Ascension* at three different time points.
- *Declination* at three different time points.
- *Local Sidereal Time* at three different time points.
- *Altitude* of the tracking station.
- *Latitude* of the tracking station.

## Usage
Open the [gauss_orbit_determination.m](./gauss_orbit_determination.m) script in MATLAB and assign data to the input variables as per your needs. Execute the script to display intermediate results and plot the determined orbit.

## Requirements
- MATLAB
- Basic understanding of orbital mechanics and celestial navigation

## Repository Structure 
- [gauss_orbit_determination.m](./gauss_orbit_determination.m): Main MATLAB script for orbit determination.
- [README.md](./README.md): Description and usage instructions.

## Example Output
Below is an example plot of the determined orbit:
![ExamplePlot](https://github.com/ShreyasP36/GaussIOD/blob/main/Example_plot.png)

## License
This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.

