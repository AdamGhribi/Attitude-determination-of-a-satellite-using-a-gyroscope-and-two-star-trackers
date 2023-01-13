%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% École Polytechnique de Montréal
% ELE6209A W2022 Final Project
% 
%             Attitude determination of a satellite
%            using a gyroscope and two star trackers
%
% Jacques Desfossés 69102
% Adam Ghribi     2161271
%
% This script implements the MEKF 
%
% History
% 01-Apr-2022 : Initial creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing/Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation time (s)
hours = 1.5;
simTime = hours * 3600;

% Time step (0.1 s)
dt = 0.1;

% Number of iterations
N = simTime/dt;

% Discrete time values
times = 0:dt:simTime;
timesMins = times/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gyroscope - Continuous Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gyro noise (rad/s)/sqrt(Hz)
sigma_g = sqrt(10)*(10^-6);

% Bias mean (in rad/s)
mu_bg = 0;

% Bias standard deviation: 0.0005°/hr for each axis (in rad/s)
stddev_bg = (0.0005/3600)*(pi/180);

% Bias noise (rad/s^2)/sqrt(Hz)
sigma_bg = sqrt(10)*(10^-10);

% Bias time constant (1/s)
lambda_bg = 1/60;

% Matrix Hg and Fg
Hg = eye(3);
Fg = -lambda_bg * eye(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Star Trackers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Noise for a single star tracker: 0.2 arcsec (in rad)
% Noise for combined star trackers is reduced by 1/2
sigma_s = (0.5) * (0.2/3600) * (pi/180);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gyroscope initial bias
bg_true = mu_bg*ones(3,1) + stddev_bg*randn(3,1);

% True initial attitude
qib_true = sqrt(2)/2 * [1;  0;  0;  1];

% True angular velocity (rad/s)
omega_b_bi_true = zeros(N+1, 3);
omega_b_bi_true(:, 1) = (pi/180)*0.1*sin(0.01*times);
omega_b_bi_true(:, 2) = (pi/180)*0.1*sin(0.0085*times);
omega_b_bi_true(:, 3) = (pi/180)*0.1*cos(0.0085*times);

% Initial state estimate
qib_hat = rpy2quat(quat2rpy(qib_true) + 0.0001*randn(3,1));
bg_hat = zeros(3,1);

% Initial covariance
Sigma11 = ((6/3600)*(pi/180))^2 * eye(3);   % Attitude error (rad^2)
Sigma22 = ((0.2/3600)*(pi/180))^2 * eye(3); % Gyro bias (rad/s)^2
Sigma12 = zeros(3,3);
Sigma21 = zeros(3,3);
Sigma = vertcat(horzcat(Sigma11, Sigma12), horzcat(Sigma21, Sigma22));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data to Record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rpyTrue = zeros(N+1, 3);
rpyTraj = zeros(N+1, 3);

rpyTrue(1, :) = dcm2rpy(quat2dcm(qib_true)'*Rin(0))';
rpyTraj(1, :) = dcm2rpy(quat2dcm(qib_hat)'*Rin(0))';

biasTrue = zeros(N+1, 3);
biasTraj = zeros(N+1, 3);

biasTrue(1,:) = bg_true';
biasTraj(1,:) = bg_hat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                      %%%
%%%%        MEKF          %%%
%%%%                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each k is 0.01 sec
for k=1:1:N
   % Gyroscope measurements 
   omega_b_ib_meas = omega_b_bi_true(k,:)' + Hg*bg_true + sigma_g*sqrt(dt)*randn(3,1);
 
   %
   % Kinematic state estimation (integration of mechanization)
   %
   omega_b_ib_hat = omega_b_ib_meas - Hg*bg_hat;
   qib_hat = integrateQuat(qib_hat, omega_b_ib_hat, dt);

   %
   % Covariance propagation : Integrate (7.16) using DT Lyapunov
   %
   A11 = antisym(-omega_b_ib_hat);
   A12 = -Hg;
   A21 =  zeros(3,3);
   A22 =  Fg;
   A = vertcat(horzcat(A11, A12), horzcat(A21, A22));

   B11 = -eye(3);
   B12 = zeros(3,3);
   B21 = zeros(3,3);
   B22 = eye(3);
   B = vertcat(horzcat(B11, B12), horzcat(B21, B22));

   Q11 = sigma_g^2 * eye(3);
   Q12 = zeros(3,3);
   Q21 = zeros(3,3);
   Q22 = sigma_bg^2 * eye(3);
   Q = vertcat(horzcat(Q11, Q12), horzcat(Q21, Q22));
   
   % Smaller time steps for covariance propagation
   Nsub = 5;
   for i=1:Nsub
      h = dt/Nsub;
      Phi = expm(A*h);
      Sigma = Phi*Sigma*Phi' + h*B*Q*B';
   end

   %
   %  True state propagation
   %
   qib_true = integrateQuat(qib_true, omega_b_bi_true(k,:)', dt);
   bg_true = expm(Fg*dt) * bg_true + sigma_bg*sqrt(dt)*randn(3,1);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Measurement update : Every 0.5 sec (every 5 gyro updates)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if mod(k,5) == 0

      % Measurement sensitivity matrix
      C = horzcat(eye(3), zeros(3,3));

      % Kalman Gain
      K = Sigma*C'/(C*Sigma*C'+ (sigma_s^2)*eye(3));

      % Covariance update
      Sigma = (eye(6) - K*C)*Sigma;

      % Star tracker measurement provides Euler angles
      %    1) Convert true attitude quaternion to Euler
      %    2) Add noise
      %    3) Convert back to quaternion for noisy attitude measurement
      qbi_meas = rpy2quat(quat2rpy(qib_true) + eye(3)*sigma_s*randn(3));

      % Measurement error (leverage Matlab quaternions)
      % Crassidis & Markley, 6.22a
      qk = quaternion(qbi_meas(1), qbi_meas(2), qbi_meas(3), qbi_meas(4));
      qk_hat = quaternion(qib_hat(1), qib_hat(2), qib_hat(3), qib_hat(4));

      % Multiplicative error for measurement
      dq = qk*(qk_hat^-1);
      [q0, q1, q2, q3] = parts(dq);
      dphimeas = 2*[q1/q0; q2/q0; q3/q0];

      % Error estimate
      dx_hat_plus = K*(dphimeas);
      dphi_plus = dx_hat_plus(1:3);
      dbg_plus = dx_hat_plus(4:6);

      % Update quaternion and perform reset
      XI1 = -qib_hat(2:4)';
      XI2 = qib_hat(1)*eye(3) - antisym(qib_hat(2:4));
      XI = vertcat(XI1, XI2);
      qib_hat_star = qib_hat + 0.5 * XI * dphi_plus;
      qib_hat = qib_hat_star/norm(qib_hat_star);

      % Update bias estimate
      bg_hat = bg_hat + dbg_plus;
   end

   %%%%%%%%%%%%%
   % Record data
   %%%%%%%%%%%%%
   Rin_k = Rin(times(k+1));
   rpyTrue(k+1, :) = dcm2rpy(quat2dcm(qib_true)'*Rin_k)';
   rpyTraj(k+1, :) = dcm2rpy(quat2dcm(qib_hat)'*Rin_k)';

   biasTrue(k+1,:) = bg_true';
   biasTraj(k+1,:) = bg_hat';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert roll/pitch/yaw (from rad to deg)
rpyTraj = rpyTraj*(180/pi);
rpyTrue = rpyTrue*(180/pi);

% Compute roll/pitch/yaw errors (from deg to arcsec)
rpyError = (rpyTraj-rpyTrue)*3600;

% Create 3-sigma lines for attitude (arcsec)
s3 = 3*sigma_s * (180/pi) * 3600;
rpy3sigplus  = [0  s3; timesMins(N+1)  s3];
rpy3sigminus = [0 -s3; timesMins(N+1) -s3];

% Convert bias (from rad/s to deg/hr)
biasTraj = biasTraj*(180/pi)*(3600);
biasTrue = biasTrue*(180/pi)*(3600);

% Compute bias error (deg/hr)
biasError = (biasTraj-biasTrue);
%%
% Create figures
figure

trueline='r';
trajline='b--';
sigline='r';

subplot(3,2,1);
sgtitle('Attitude Trajectory','Interpreter','latex');

plot(timesMins, rpyTrue(:,1), trueline, timesMins, rpyTraj(:,1), trajline);
legend('$\phi_{true}$', '$\phi$','Interpreter','latex');
title('Roll','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Orientation (deg)','Interpreter','latex');
grid 
subplot(3,2,2);
plot(timesMins, rpyError(:,1),  rpy3sigplus(:,1), rpy3sigplus(:,2), sigline, ...
     rpy3sigminus(:,1), rpy3sigminus(:,2), sigline);
axis([0 timesMins(N+1) -5*s3 5*s3])
legend('Error', '3$\sigma_{STR}$','Interpreter','latex');
title('Roll Error','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Error (arcsec)','Interpreter','latex');
grid
subplot(3,2,3);
plot(timesMins, rpyTrue(:,2), trueline, timesMins, rpyTraj(:,2), trajline);
legend('$\theta_{true}$', '$\theta$','Interpreter','latex');
title('Pitch','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Orientation (deg)','Interpreter','latex');
grid
subplot(3,2,4);
plot(timesMins, rpyError(:,2), rpy3sigplus(:,1), rpy3sigplus(:,2), sigline, ...
     rpy3sigminus(:,1), rpy3sigminus(:,2), sigline);
axis([0 timesMins(N+1) -5*s3 5*s3])
legend('Error', '3$\sigma_{STR}$','Interpreter','latex');
title('Pitch Error','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Error (arcsec)','Interpreter','latex');
grid
subplot(3,2,5);
plot(timesMins, rpyTrue(:,3), trueline, timesMins, rpyTraj(:,3), trajline);
legend('$\psi_{true}$', '$\psi$','Interpreter','latex');
title('Yaw','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Orientation (deg)','Interpreter','latex');
grid
subplot(3,2,6);
plot(timesMins, rpyError(:,3), rpy3sigplus(:,1), rpy3sigplus(:,2), sigline, ...
     rpy3sigminus(:,1), rpy3sigminus(:,2), sigline);
axis([0 timesMins(N+1) -5*s3 5*s3])
legend('Error', '3$\sigma_{STR}$','Interpreter','latex');
title('Yaw Error','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex')
ylabel('Error (arcsec)','Interpreter','latex');
grid
figure

subplot(3,2,1);
sgtitle('Bias Trajectory','Interpreter','latex');

plot(timesMins, biasTrue(:,1), trueline, timesMins, biasTraj(:,1), trajline);
legend('Bias$_{true}$', '$\hat{Bias}$','Interpreter','latex');
title('Bias x-axis','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Bias (deg/hr)','Interpreter','latex');
grid
subplot(3,2,2);
plot(timesMins, biasError(:,1));
title('Bias Error x-axis','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Bias (deg/hr)','Interpreter','latex');
grid
subplot(3,2,3);
plot(timesMins, biasTrue(:,2), trueline, timesMins, biasTraj(:,2), trajline);
legend('Bias$_{true}$', '$\hat{Bias}$','Interpreter','latex');
title('Bias y-axis','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Bias (deg/hr)','Interpreter','latex');
grid
subplot(3,2,4);
plot(timesMins, biasError(:,2));
title('Bias Error y-axis','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Bias (deg/hr)','Interpreter','latex');
grid
subplot(3,2,5);
plot(timesMins, biasTrue(:,3), trueline, timesMins, biasTraj(:,3), trajline);
legend('Bias$_{true}$', '$\hat{Bias}$','Interpreter','latex');
title('Bias z-axis','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Bias (deg/hr)','Interpreter','latex');
grid
subplot(3,2,6);
plot(timesMins, biasError(:,3));
title('Bias Error z-axis','Interpreter','latex');
xlabel('Time (min)','Interpreter','latex');
ylabel('Bias (deg/hr)','Interpreter','latex');
grid
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotation from {i} to {n}
function R=Rin(t)

   % Gravitational parameter of the Earth : 3.986*10^5 km^3/s^2
   mu = 3.986*10^5;

   % Radius : Distance of satellite w.r.t center of the Earth
   %          R = Radius of Earth + Altitude (km)
   rE = 6378;
   alt = 600;
   R = rE + alt;

   theta = sqrt(mu/R^3)*t;
   R = [ sin(theta) 0  cos(theta)
        -cos(theta) 0  sin(theta)
             0      -1     0     ];
   % Uncomment if attitude is desired as Rbi instead of Rbn
   % R=eye(3);
end

% Integration Quaternion Kinematic Equations (6.27)
function qi=integrateQuat(q, w, h)
   n = norm(w);
   Q = [0    -w(1) -w(2) -w(3);
        w(1)  0     w(3) -w(2);
        w(2) -w(3)  0     w(1);
        w(3)  w(2) -w(1)  0 ];
   qi = (cos(0.5*h*n)*eye(4) + (sin(0.5*h*n)/n)*Q)*q;
end

% Rotation matrix from quaternion (4.31)
function R=quat2dcm(q)
   R=zeros(3,3);
   q0=q(1); q1=q(2); q2=q(3); q3=q(4);
   R(1,1) = q0^2 + q1^2 - q2^2 - q3^2;
   R(1,2) = 2*(q1*q2 - q0*q3);
   R(1,3) = 2*(q0*q2 + q1*q3);
   R(2,1) = 2*(q1*q2 + q0*q3);
   R(2,2) = q0^2 - q1^2 + q2^2 - q3^2;
   R(2,3) = 2*(q2*q3 - q0*q1);
   R(3,1) = 2*(q1*q3 - q0*q2);
   R(3,2) = 2*(q0*q1 + q2*q3);
   R(3,3) = q0^2 - q1^2 - q2^2 + q3^2;
end

% Skew-symmetric matrix
function w_x = antisym(w)
%   assert(isvec(w,3));
   w_x = [ 0    -w(3)  w(2)
           w(3)  0    -w(1)
          -w(2)  w(1)  0   ];
end

% Quaternion to Euler (RPY), Equations (4.32, 4.33, 4.34)
function rpy = quat2rpy(q)
   q0=q(1); q1=q(2); q2=q(3); q3=q(4);
   
   roll  =  atan2(2*(q0*q1+q2*q3), 1-2*(q1^2+q2^2));
   pitch = -asin(2*(q1*q3-q0*q2));
   yaw   =  atan2(2*(q1*q2+q0*q3), 1-2*(q2^2+q3^2));

   rpy = [roll; pitch; yaw];
end

% DCM to RPY, Equations (4.18, 4.19, 4.20)
function rpy = dcm2rpy(R)
   roll  = atan2(R(3,2), R(3,3));
   pitch = -asin(R(3,1));
   yaw   = atan2(R(2,1), R(1,1));

   rpy = [roll; pitch; yaw];
end

% Euler (RPY) to Direction Cosine Matrix (DCM), Equation (4.17)
function R=rpy2dcm(rpy)
   phi   = rpy(1);
   theta = rpy(2);
   psi   = rpy(3);
   R = zeros(3,3);
   R(1,1) = cos(theta)*cos(psi);
   R(1,2) = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
   R(1,3) = sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
   R(2,1) = cos(theta)*sin(psi);
   R(2,2) = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi);
   R(2,3) = -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
   R(3,1) = -sin(theta);
   R(3,2) = sin(phi)*cos(theta);
   R(3,3) = cos(phi)*cos(theta);
end

% Euler (RPY) to Quaternion, Equation (p. 42).
function q=rpy2quat(rpy)
   phi   = rpy(1);
   theta = rpy(2);
   psi   = rpy(3);
   q = zeros(4,1);
   q(1) = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2);
   q(2) = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2);
   q(3) = cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2);
   q(4) = cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2);
end

% Quaternion to DCM (4.37)
function q = dcm2quat(R)
   q0 = sqrt(1 + R(1,1) + R(2,2) + R(3,3));
   q1 = sign(R(3,2) - R(2,3)) * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
   q2 = sign(R(1,3) - R(3,1)) * sqrt(1 - R(1,1) + R(2,2) - R(3,3));
   q3 = sign(R(2,1) - R(1,2)) * sqrt(1 - R(1,1) - R(2,2) + R(3,3));
   q = (1/2) * [q0; q1; q2; q3];
end