% Navigation and Guide Systems Projet
% Master Degree in Robotics and Automation Engineering, Dipartimento di Ingegneria, Universita´ di Pisa
% Authors: Rachele Nebbia Colomba, Chiara Sammarco
%    copyright            : (C) 2021 Dipartimento di Ingegneria dell´Informazione (DII) // Universita´ di Pisa    
%    email                : rachelenebbia <at> gmail <dot> com

%/ ---------------------------- EKF 1 DOF for bias estimation ------------------ //

%% Init FIle for EKF 1dof
%% here you can find the parameter needed to run the simulation file KF_1_gdl.slx

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

%sampling time white noise added to the measurements
Ts_noise = 1/50; 

% Constant accelemoter bias 
bias_body = 0.5;

% measurement noises transfer function
filter_pos = tf(1,[1 1]) %noise on fictious position data
filter_vel = tf(1,[1 1]); %noise on fictious velocity data
filter_acc = tf(1,[0.5 1]); %noise on fictious acceleration data
noise_gain = 0.01;

%EKF weights
Q = single(diag([1e-2, 1e-6, 1e-4, 1e-12]))*2; 

%Covariance Values for sensor
GPS_cov = single([1*2, 1])/2; %GPS covariance no overshoot
IMU_cov = single(15*0.07)/2; %IMU covariance
BARO_cov = single(1.7)/2; %baro covariance

% standar filter, load for comparison 
GPS_cov_old  = single([0.1, 0.1, 200, 200, 0.36])/2; %GPS covariance
IMU_cov_old = single([15*[0.07, 0.07] 10*0.09])/2; %IMU covariance
BARO_cov_old = single(1.7)/2; %baro covariance
Q_old = single(diag([1e-8, 1e-1, 1e-2, 1e-8, 1e-1, 1e-2, 1e-5, 1e-4, 0.010]))*2; 
