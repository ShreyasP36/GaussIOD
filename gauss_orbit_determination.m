%% Gauss Method of preliminary orbit determination

%% Tracking station position

H = 1; % Altitude in km

phi = 40.28; % Latitude (in degrees)

%% Tracking Station Observations

% Time at three different timepoints (in s)
t1 = 0; 
t2 = 118.10;
t3 = 237.58;

% Right ascension at three different time points (in degrees)
alpha_1 = 43.537;
alpha_2 = 54.420;
alpha_3 = 64.318;

% Declination at three different time points (in degrees)
d1 = -8.7833;
d2 = -12.074;
d3 = -15.105;

% Local sidereal time at three different time points (in degrees)
theta_1 = 44.506;
theta_2 = 45;
theta_3 = 45.499;

%% Required Constants

Re = 6378; % Equatorial radius of earth (in km)
flat = 0.003353; % Flattening factor
u = 398600; % Graviatational Parameter

%% Convert angles from degrees to radians

phi = deg2rad(phi);
theta_1 = deg2rad(theta_1);
theta_2 = deg2rad(theta_2);
theta_3 = deg2rad(theta_3);
d1 = deg2rad(d1);
d2 = deg2rad(d2);
d3 = deg2rad(d3);
alpha_1 = deg2rad(alpha_1);
alpha_2 = deg2rad(alpha_2);
alpha_3 = deg2rad(alpha_3);

%% Inertial Position Vector of tracking station at 3 observation times

factor = sqrt(1 - (2*flat - flat^2) * (sin(phi))^2);
Re_factor = Re / factor;
R_x = @(theta) (Re_factor + H) * cos(phi) * cos(theta);
R_y = @(theta) (Re_factor + H) * cos(phi) * sin(theta);
R_z = @(theta) ((Re_factor*(1 - flat)^2)+ H) * sin(phi);

R1 = [R_x(theta_1), R_y(theta_1), R_z(theta_1)];  
R2 = [R_x(theta_2), R_y(theta_2), R_z(theta_2)]; 
R3 = [R_x(theta_3), R_y(theta_3), R_z(theta_3)]; 

display(R1)
display(R2)
display(R3)

%% Computing directional cosine vectors

unit_rho_x = @(d,alpha) (cos(d)*cos(alpha));
unit_rho_y = @(d,alpha) (cos(d)*sin(alpha));
unit_rho_z = @(d) (sin(d));

% Calculation of unit vectors in the direction of the slant range vector
unit_rho_1 = [unit_rho_x(d1,alpha_1) unit_rho_y(d1,alpha_1) unit_rho_z(d1)]; 
unit_rho_2 = [unit_rho_x(d2,alpha_2) unit_rho_y(d2,alpha_2) unit_rho_z(d2)];
unit_rho_3 = [unit_rho_x(d3,alpha_3) unit_rho_y(d3,alpha_3) unit_rho_z(d3)];

display(unit_rho_1)
display(unit_rho_2)
display(unit_rho_3)

%% Computation of scalar quantities

% Calculation of time intervals
T_1 = t1-t2;
T_3 = t3-t2;
T = t3-t1;

% Calculation of cross products
p1 = cross(unit_rho_2,unit_rho_3);
p2 = cross(unit_rho_1,unit_rho_3);
p3 = cross(unit_rho_1,unit_rho_2);
display(p1)
display(p2)
display(p3)

% Calculation of dot products
D0 = dot(unit_rho_1,p1);
D11 = dot(R1,p1);
D12 = dot(R1,p2);
D13 = dot(R1,p3);
D21 = dot(R2,p1);
D22 = dot(R2,p2);
D23 = dot(R2,p3);
D31 = dot(R3,p1);
D32 = dot(R3,p2);
D33 = dot(R3,p3);

%% Calculation of parameters A, B and E

A = (1/D0)*((-D12*(T_3/T))+D22+(D32*(T_1/T)));
B = (1/(6*D0))*((D12*(T_3^2-T^2)*(T_3/T))+(D32*(T^2-T_1^2)*(T_1/T)));
E = dot(R2,unit_rho_2);
R_sq = dot(R2,R2);
display(A)
display(B)
display(E)

%% Calculation of coefficients a,b and c of F(x) 

a = -(A^2+(2*A*E)+R_sq);
b = -2*u*B*(A+E);
c = -(u^2)*(B^2);
fprintf('The eigth order polynomial: x^8 + (%d)x^6 + (%d)x^3 + %d\n',a,b,c)

%% Finding roots of F(x) 

% Define the coefficients of the polynomial
coefficients = [1, 0, a, 0, 0, b, 0, 0, c];

% Solve for the roots
roots_of_polynomial = roots(coefficients);

% Display all roots
disp('All roots of the polynomial are:');
disp(roots_of_polynomial);

% Filter out the real roots
real_roots = roots_of_polynomial(imag(roots_of_polynomial) == 0);

% Select only the positive real roots
r_2= real_roots(real_roots > 0);
display(r_2)

%% Computing initial slant ranges 

rho_1_numerator = 6 * r_2^3 * (D31 * (T_1 / T_3) + D21 * (T / T_3)) + u * D31 * (T^2 - T_1^2) * (T_1 / T_3);
rho_1_denominator = 6 * r_2^3 + u * (T^2 - T_3^2);
rho_1 = (1 / D0) * (rho_1_numerator / rho_1_denominator - D11);
display(rho_1)
rho_2 = A + u * (B/r_2^3);
display(rho_2)
rho_3_numerator = 6 * r_2^3 * (D13 * (T_3 / T_1) - D23 * (T / T_1)) + u * D13 * (T^2 - T_3^2) * (T_3 / T_1);
rho_3_denominator = 6 * r_2^3 + u * (T^2 - T_1^2);
rho_3 = (1 / D0) * (rho_3_numerator / rho_3_denominator - D33);
display(rho_3)

%% Computing initial position vectors

r1 = R1 + rho_1*unit_rho_1;
r2 = R2 + rho_2*unit_rho_2;
r3 = R3 + rho_3*unit_rho_3;

display(r1)
display(r2)
display(r3)

%% Calculation of Lagrange coefficients

f1 = 1 - 0.5*(u/r_2^3)*T_1^2;
f3 = 1 - 0.5*(u/r_2^3)*T_3^2;
g1 = T_1 - (1/6)*(u/r_2^3)*T_1^3;
g3 = T_3 - (1/6)*(u/r_2^3)*T_3^3;

display(f1)
display(f3)
display(g1)
display(g3)

%% Computation of velocity v2

v2 = (1/(f1 * g3 - f3 * g1)) * (-f3 * r1 + f1 * r3);
display(v2)

%% Iterative algorithm for improvement based on tolerance

tol = 1e-3; % Define the tolerance
max_iter = 10; % Maximum number of iterations
iter_count = 0;

while true
    % Previous values of rho
    prev_rho_1 = rho_1;
    prev_rho_2 = rho_2;
    prev_rho_3 = rho_3;
    
    r_2 = sqrt(dot(r2, r2));
    v_2 = sqrt(dot(v2, v2));
    alpha_0 = (2 / r_2) - ((v_2^2) / u);
    vr2 = dot(v2, r2) / r_2;
    
    C = @(x) arrayfun(@(x) (x > 0) * (1 - cos(sqrt(x))) / x + (x == 0) * 0.5 + (x < 0) * (cosh(sqrt(-x)) - 1) / -x, x);
    S = @(x) arrayfun(@(x) (x > 0) * (sqrt(x) - sin(sqrt(x))) / (sqrt(x)^3) + (x == 0) * (1 / 6) + (x < 0) * (sinh(sqrt(-x)) - sqrt(-x)) / (sqrt(-x)^3), x);
    
    xi1 = sqrt(u) * T_1 / r_2;
    xi3 = sqrt(u) * T_3 / r_2;
    
    % Newton-Raphson method for xi1
    for iter = 1:max_iter
        F1 = sqrt(u) * T_1 - (r_2 * vr2 / sqrt(u)) * xi1^2 * C(alpha_0 * xi1^2) - (1 - alpha_0 * r_2) * xi1^3 * S(alpha_0 * xi1^2) - r_2 * xi1;
        F1_prime = -2 * (r_2 * vr2 / sqrt(u)) * xi1 * C(alpha_0 * xi1^2) - (3 * (1 - alpha_0 * r_2) * xi1^2 * S(alpha_0 * xi1^2) + (r_2 * vr2 / sqrt(u)) * xi1^2 * alpha_0 * xi1 * C(alpha_0 * xi1^2) - (1 - alpha_0 * r_2) * xi1^3 * alpha_0 * xi1 * S(alpha_0 * xi1^2)) - r_2;
        delta_xi1 = -F1 / F1_prime;
        xi1 = xi1 + delta_xi1;
        if abs(delta_xi1) < tol
            break;
        end
    end
    
    % Newton-Raphson method for xi3
    for iter = 1:max_iter
        F3 = sqrt(u) * T_3 - (r_2 * vr2 / sqrt(u)) * xi3^2 * C(alpha_0 * xi3^2) - (1 - alpha_0 * r_2) * xi3^3 * S(alpha_0 * xi3^2) - r_2 * xi3;
        F3_prime = -2 * (r_2 * vr2 / sqrt(u)) * xi3 * C(alpha_0 * xi3^2) - (3 * (1 - alpha_0 * r_2) * xi3^2 * S(alpha_0 * xi3^2) + (r_2 * vr2 / sqrt(u)) * xi3^2 * alpha_0 * xi3 * C(alpha_0 * xi3^2) - (1 - alpha_0 * r_2) * xi3^3 * alpha_0 * xi3 * S(alpha_0 * xi3^2)) - r_2;
        delta_xi3 = -F3 / F3_prime;
        xi3 = xi3 + delta_xi3;
        if abs(delta_xi3) < tol
            break;
        end
    end
    
    % Update Lagrange coefficients
    f1 = 1 - (xi1^2/r_2) * C(alpha_0 * xi1^2);
    f3 = 1 - (xi3^2/r_2) * C(alpha_0 * xi3^2);
    g1 = T_1 - (xi1^3/u) * S(alpha_0 * xi1^2);
    g3 = T_3 - (xi3^3/u) * S(alpha_0 * xi3^2);
    
    % Update c1 and c3
    c1 = g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);
    
    % Update rho values
    rho_1 = ((-D11) + (D21/c1) - ((c3 / c1) * D31))/D0;
    rho_2 = (-(c1 * D12) + (D22) + (- c3 * D32))/D0;
    rho_3 = (((-c1 / c3 )* D13)+ (D23 / c3) - D33)/D0;
    
    % Check if rho values are valid
    if any(isnan([rho_1, rho_2, rho_3]))
        error('NaN values detected in rho calculations.');
    end
    
    % Update position vectors
    r1 = R1 + rho_1*unit_rho_1;
    r2 = R2 + rho_2*unit_rho_2;
    r3 = R3 + rho_3*unit_rho_3;
    
    % Update velocity vector
    v2 = (1/(f1 * g3 - f3 * g1)) * (-f3 * r1 + f1 * r3);
    
    % Check for convergence
    if abs(rho_1 - prev_rho_1) < tol && abs(rho_2 - prev_rho_2) < tol && abs(rho_3 - prev_rho_3) < tol
        break;
    end
    
    iter_count = iter_count + 1;
    if iter_count >= max_iter
        break;
    end
    
    % Display intermediate results
    disp(['Iteration ', num2str(iter_count), ':']);
    disp(['rho_1: ', num2str(rho_1)]);
    disp(['rho_2: ', num2str(rho_2)]);
    disp(['rho_3: ', num2str(rho_3)]);
end

display(r2)
display(v2)

%% Obtaining Orbital elements from the state vector at t2

% Calculation of specific angular momentum
h = cross(r2,v2);
h_1 = sqrt(dot(h,h));
display(h_1)

% Calculation of inclination
incl = acosd(h(1,3)/h_1);
display(incl)

%Calculation of ascending node
z = [0 0 1];
N = cross(z,h);
N_1 = sqrt(dot(N,N));
if N(1,2)>=0
    ascen = acosd(N(1,1)/N_1);
else
    ascen = 360 - acosd(N(1,1)/N_1);
end
display(ascen)

% Calculation of eccentricity vector and eccentricity
ecc = (1/u) * ((v_2^2 - (u/r_2))*r2 - r_2*vr2*v2 );
e = sqrt(dot(ecc,ecc));
display(e)

% Calculation of argument of perigee
if ecc(1,3)>=0
    argu = acosd(dot(N,ecc)/(N_1*e));
else
    argu = 360 - acosd(dot(N,e)/(N_1*e));
end
display(argu)
    
% Calculation of true anomaly
if vr2>=0
    theta = acosd(dot(ecc,r2)/(e*r_2));
else
    theta = 360 - acosd(dot(ecc,r2)/(e*r_2));
end
display(theta)

% Semi-major axis
a = (1 / (2 / r_2 - v_2^2 / u)); 
display(a)

% Number of points for plotting the orbit
num_points = 1000;

% True anomaly range for a complete orbit (0 to 360 degrees)
theta_vals = linspace(0, 2 * pi, num_points);

% Pre-allocate position vectors
orbit_positions = zeros(3, num_points);

for i = 1:num_points
    theta_i = theta_vals(i);

    % Compute the radius at each true anomaly
    r_i = (a * (1 - e^2)) / (1 + e * cos(theta_i));

    % Compute the position in the perifocal coordinate system
    r_perifocal = [r_i * cos(theta_i); r_i * sin(theta_i); 0];

    % Rotation matrices
    R3_w = [cosd(argu), sind(argu), 0; -sind(argu), cosd(argu), 0; 0, 0, 1];
    R1_i = [1, 0, 0; 0, cosd(incl), sind(incl); 0, -sind(incl), cosd(incl)];
    R3_O = [cosd(ascen), sind(ascen), 0; -sind(ascen), cosd(ascen), 0; 0, 0, 1];

    % Combined rotation matrix from perifocal to geocentric equatorial frame
    Q_pX = (R3_O * R1_i * R3_w)';

    % Compute the position in the geocentric equatorial frame
    r_geocentric = Q_pX * r_perifocal;

    % Store the position
    orbit_positions(:, i) = r_geocentric;
end

% Plotting the orbit
figure;
plot3(orbit_positions(1, :), orbit_positions(2, :), orbit_positions(3, :), 'b', 'LineWidth', 1.5);
hold on;

% Adding the equatorial plane
theta_plane = linspace(0, 2*pi, 100);
r_plane = max(max(abs(orbit_positions))) * 1.1; % Adjust the radius of the plane based on the orbit
x_plane = r_plane * cos(theta_plane);
y_plane = r_plane * sin(theta_plane);
z_plane = zeros(size(theta_plane));

fill3(x_plane, y_plane, z_plane, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the Earth at the origin
plot3(0, 0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Earth');

% Enhance the plot
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Initial Orbit Plot');
legend('Orbit','Equatorial Plane','Origin', 'Location', 'best');
grid on;
axis equal;
view(3);






