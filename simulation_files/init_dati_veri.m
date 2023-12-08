% Navigation and Guide Systems Projet
% Master Degree in Robotics and Automation Engineering, Dipartimento di Ingegneria, Universita´ di Pisa
% Authors: Rachele Nebbia Colomba, Chiara Sammarco
%    copyright            : (C) 2022 Technical University of Munich // Universitä di pisa    
%    email                : rachelenebbia <at> gmail <dot> com

%/ ---------------------------- EKF 3 DOF for bias estimation ------------------ //

%%Init FIle for EKF 3dof with real data 
%here you can find the parameter needed to run the simulation file KF_3_Dati_Veri.slx

%clear and close all workspace
clear all;close all;clc
%Sampling time
ST = 1/200;
T_50 = 1/50;

display('Loading data..')

%Load sensors dataset from file 
load('LOG00079_parsed_seg2.mat')
display('Loading plot_data..')

%run script to visualize the sensors dataset
run plot_data.m

%%Parameters
g_ned = [0; 0; -9.81]; %gravity vector 

% Constant accelerometer biases on 3 axes
b_x_vero = -0.3;
b_y_vero = 0.76;
b_z_vero = 0.47;

% Noise parameters
Ts_noise = 1/50; 
noise_gain = 0.01;

%noise variances on 3-axes position data
var_p_x = 2;
var_p_y = 2;
var_p_z = 0.30; 

%noise variances on 3-axes velocity data
var_v_x = 0.60;
var_v_y = 0.32;
var_v_z = 0.60;

%Gryscope variance vector
var_giroscopio  = [(deg2rad(0.16))^2;(deg2rad(0.16))^2;(deg2rad(0.16))^2];          % wx,wy,wz

% Low-pass filter for measurements
filter_acc = tf(1,[0.5 1]); 
filter_vel = tf(1,[1 1]);
filter_pos = tf(1,[1 1]);
filter_pos_Z = tf(1,[0.5 1]);

%%New EKF Parameters
px = 1e-4;
vx = 1e-6;
ax = 1e-3;

py = 1e-4;
vy = 1e-6;
ay = 1e-3;

pz = 1e-8; 
vz = 1e-9;
az = 1e-2;

bx = 1e-20;
by = 1e-20;
bz = 1e-20;

Q = single(diag([px,vx,ax,py,vy,ay,pz,vz,az,bx,by,bz]))*2;

% sensors covariances
% GPS cov px|py|vx|vy|vz
GPS_cov  = single([0.1, 0.1, 200*5, 200*5, 0.36])/2; %GPS covariance
IMU_cov = single([15*[0.07*5, 0.07*5] 15*50*0.09])/2; %IMU covariance
BARO_cov = single(1.7)/5; %baro covariance

% Standard filter weights
Q_old = single(diag([1e-8, 1e-1, 1e-2, 1e-8, 1e-1, 1e-2, 1e-5, 1e-4, 0.010]))*2;
GPS_cov_old  = single([0.1, 0.1, 200, 200, 0.36]/2); 
IMU_cov_old = single([15*[0.07, 0.07] 10*0.09]/2); %IMU covariance
BARO_cov_old = single(1.7)/2; %baro covariance