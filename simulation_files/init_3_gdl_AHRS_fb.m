% Navigation and Guide Systems Projet
% Master Degree in Robotics and Automation Engineering, Dipartimento di Ingegneria, Universita´ di Pisa
% Authors: Rachele Nebbia Colomba, Chiara Sammarco
%    copyright            : (C) 2021 Dipartimento di Ingegneria dell´Informazione (DII) // Universita´ di Pisa    
%    email                : rachelenebbia <at> gmail <dot> com

%/ ---------------------------- EKF 3 DOF for bias estimation with feddback correction ------------------ //

%%Init FIle for EKF 3dof
%here you can find the parameter needed to run the simulation file KF_3_gdl_AHRS_Dataset1.slx + KF_3_gdl_AHRS_Dataset2.slx

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

%% Parameters
g_ned = [0; 0; -9.81]; %gravity vector 

%%Accelerometer Biases fixed values 
% Dataset 1
b_x = -0.3;
b_y = 0.2;
b_z = 0.5;
% Dataset 2 
bj_x = 0.4;
bj_y = 0.9; 
bj_z = 0.5;

% Noise parameters
Ts_noise = 1/50; 
noise_gain = 0.01;

%noise variances on 3-axes position data
var_p_x = 2;
var_p_y = 2;
var_p_z = 0.30; 

%noise variances on 3-axes velocitydata
var_v_x = 0.60;
var_v_y = 0.32;
var_v_z = 0.60;

%gyroscope variance vector
var_giroscopio  = [(deg2rad(0.16))^2;(deg2rad(0.16))^2;(deg2rad(0.16))^2];          % wx,wy,wz

% Low-pass filter for measurements
filter_acc = tf(1,[0.5 1]); 
filter_vel = tf(1,[1 1]);
filter_pos = tf(1,[1 1]);
filter_pos_Z = tf(1,[0.5 1]);

%%New EKF Parameters on states p=postion,v=velocity,a=acceleration,b=bias
px = 1e-8;
vx = 1e-9;
ax = 1e-5;

py = 1e-8;
vy = 1e-9;
ay = 1e-5;

pz = 1e-8; 
vz = 1e-9;
az = 1e-6;

bx = 1e-12;
by = 1e-12;
bz = 1e-12;

Q = single(diag([px,vx,ax,py,vy,ay,pz,vz,az,bx,by,bz]))*2;

% covariances
% GPS cov px|py|vx|vy|vz
GPS_cov  = single([0.1*50, 0.1*50, 200*5, 200*5, 0.36])/2; %GPS covariance
% GPS_cov  = single([1, 1, 1, 1, 0.36])/2;
IMU_cov = single([15*[0.07*50, 0.07*50] 15*50*0.09])/2; %IMU covariance
BARO_cov = single(1.7)/5; %baro covariance

% Standard filter weights
Q_old = single(diag([1e-8, 1e-1, 1e-2, 1e-8, 1e-1, 1e-2, 1e-5, 1e-4, 0.010]))*2;

% Sensors covariances 
GPS_cov_old  = single([0.1, 0.1, 200, 200, 0.36]/2); 
IMU_cov_old = single([15*[0.07, 0.07] 10*0.09]/2); %IMU covariance
BARO_cov_old = single(1.7)/2; %baro covariance



















bjy = bj_y - 0.1; 
