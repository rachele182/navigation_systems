% Navigation and Guide Systems Projet
% Master Degree in Robotics and Automation Engineering, Dipartimento di Ingegneria, Universita´ di Pisa
% Authors: Rachele Nebbia Colomba, Chiara Sammarco
%    copyright            : (C) 2021 Dipartimento di Ingegneria dell´Informazione (DII) // Universita´ di Pisa    
%    email                : rachelenebbia <at> gmail <dot> com

%/ ---------------------------- EKF 1 DOF for bias estimation ------------------ //

%%Init FIle for EKF 3dof
%here you can find the parameter needed to run the simulation file KF_3_gdl.slx

%clear and close all workspace
clear;close;clc
%Sampling time and frequency
ST = 1/200;
T_50 = 1/50;

display('Loading data..')

%Load sensors dataset from file 
load('LOG00079_parsed_seg2.mat')
display('Loading plot_data..')

%run script to visualize the sensors dataset
run plot_data.m

% Measurement noises: sampling time, variances
Ts_noise = 1/50; 
noise_gain = 0.001;

%noise variances on 3-axes position data
var_p_x = 1; 
var_p_y = 1;
var_p_z = 0.15; 

%noise variances on 3-axes velocity data
var_v_x = 0.30; 
var_v_y = 0.16;
var_v_z = 0.30;

%Measurement noises transfer function 
filter_pos = tf(1,[1 1]);
filter_pos_Z = tf(1,[0.5 1]);
filter_vel = tf(1,[1 1]);
filter_acc = tf(1,[0.5 1]);

% Constant accelerometer biases on 3 axes
b_x = -0.03;
b_y = 0.02;
b_z = 0.05;

% EKF filter parameters

%standard filter filter
Q_old = single(diag([1e-8, 1e-1, 1e-2, 1e-8, 1e-1, 1e-2, 1e-5, 1e-4, 0.010]))*2;
% GPS cov px|py|vx|vy|vz
GPS_cov_old  = single([0.1, 0.1, 200, 200, 0.36])/2; %GPS covariance
IMU_cov_old = single([15*[0.07, 0.07] 10*0.09])/2; %IMU covariance
BARO_cov_old = single(1.7)/2; %baro covariance

%new filter
px = 1e-2;
vx = 1e-6;
ax = 1e-4;

py = 1e-2;
vy = 1e-6;
ay = 1e-4;

pz = 1e-2; 
vz = 1e-6;
az = 1e-4;

bx = 1e-12;
by = 1e-12;
bz = 1e-12;

Q = single(diag([px,vx,ax,py,vy,ay,pz,vz,az,bx,by,bz]))*2;
 
% Sensors covariances 
GPS_cov = single([1*2, 1*2, 1, 1, 0.36])/2;
IMU_cov = single([15*[0.07, 0.07] 10*0.09])/2; %IMU covariance
BARO_cov = single(1.7)/2; %baro covariance
