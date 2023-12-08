% Master Degree in Robotics and Automation Engineering, Dipartimento di Ingegneria, Universita´ di Pisa
% Authors: Rachele Nebbia Colomba, Chiara Sammarco
%    copyright            : (C) 2022 Technical University of Munich // Universitä di pisa    
%    email                : rachelenebbia <at> gmail <dot> com

%%Description: Function to calculate the tranformation matrix to rotate from body frame to NED frame
%Inputs: phi, teta, pis [rad] --> body asset euler angles
%Outpit: C_bn = 3x3 matrix 

function C_bn = body2ned(phi,theta,psi)
    phi = single(phi);
    theta = single(theta);
    psi = single(psi);
    % matrici di rotazione elementari trasposte     
    C_x_t = [1,         0,         0;
             0,  cos(phi),  -sin(phi);
             0,  sin(phi),  cos(phi)];
    
    C_y_t = [ cos(theta),   0,   sin(theta);
                       0,   1,            0;
             -sin(theta),   0,  cos(theta)];   
    
    C_z_t = [ cos(psi), -sin(psi),        0;
              sin(psi),  cos(psi),        0;
                     0,         0,       1];
    
    C_bn = C_z_t*C_y_t*C_x_t;
end