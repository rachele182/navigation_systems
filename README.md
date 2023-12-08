### <font color="green"> <span style="font-size:larger;"> Navigation Systems Project </font> </span>
<font color="green">**Authors**:</font>  Rachele Nebbia Colomba, Chiara Sammarco  
<font color="green">**Title**: Extended Kalaman Filter for Accelerometer bias estimation </font> 

This repositoty contains the code and results of simulation of an Extended Kalman Filter designed for a drone system to estimate the accelerometer bias (added here as a state of the system).  
The Project was part of *Guide and Navigation System* master degree course at the Department of Robotics and Automation Engineering, Universita´di Pisa.  
The EFK was designed and tested using navigation toolbox in **Matlab** and **Simulink**.  
The work include the simulation of the proposed filter with a real-dataset collected at the Universita´di Pisa during different fly experiments of a drone. 

Below a


The work is organized as follows:  

&#x1F539; **EKF 1 dof**, **EKF 3 dof**: the extended Kalman filter is designed first on simplified 1 degree of freedom system and three degree of freedom system, where the dataset was created in the same simulation environment;

&#x1F539; **EKF 3 dof with feedback correction**: the extended Kalman filter is tested and the estimated bias is given as a correction input to the asset computation (AHRS system);

&#x1F539; **EKF with feedback correction tested on real dataset**: the extended Kalman filter is tested on the real_dataset to evaluate the perfomance;

  - from scratch CNN,
  - pretrained CNN,
  - two-channel CNN using baseline patches (to pair in the NN bening tissue with corrispondend malign one )


For each test a simulink simulation was created and the results of proposed EKF was compared to a standard kalman filter without the bias estimation.

To have an explanation of the matlab files and how to run the simulink simulation please refere to the _contents.md_ file inside the folder. 

