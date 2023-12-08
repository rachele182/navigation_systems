### <font color="green"> <span style="font-size:larger;"> Contents of Four_Classes_Classification: </font> </span>

Here afther a brief description of the contents of the folder.   


- **Simulation Files**:
    - KF_1gdl.slx and KF_3gdl.slx contains the preliminar tests on the simulated   dataset dseigned for simplifies 1 degree and 3 degree of freedom respectively;
    - KF_3_gdl_AHRS_Datase1.slx and KF_3_gdl_AHRS_Dataset2.slx contains the simulation files with two different created dataset, with the feedback of the estimated bias used to correct asset AHRS estimation;
    - KF_AHRS_Dati_Veri.slx refers to the simulation run with the real dataset collected on a fly-trial run in February 2020. The data are loaded fro, *LOG00079_parsed_seg2 in test_volo_20200218*.

- **How to run**:
    1. run first matlab script init_1gdl.m --> open KF_1_gdl.slx --> start simulink simulation;
    2. run first matlab script init_3gdl.m --> open KF_3_gdl.slx --> start simulink simulation;
    3. run first matlab script init_3_gdl_AHRS_fb.m --> open KF_3_gdl_AHRS_Datase1.slx, KF_3_gdl_AHRS_Dataset2.slx --> start simulink simulation. 
	4. run first matlab script --> open init_dati_veri.nm --> KF_AHRS_Dati_Veri.slx --> start simulink simulation.

Please note that all the files are written in Matlab and simulink.  
The Kalman filter as well all the other functions are written as *MatlabFunctions* inside the simulink blocks.
In each cell you can visualize the obtained result. Please seutp workspace and collab notebook to run locally the codes. 