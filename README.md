<div align="center">

# Navigation Systems Project — EKF for Accelerometer Bias Estimation ✈️

**Extended Kalman Filter to estimate accelerometer bias for a drone navigation system**

**Authors:** Rachele Nebbia Colomba · Chiara Sammarco

*Project — MSc course "Guidance and Navigation Systems" · Department of Robotics and Automation Engineering, Università di Pisa*

[![MATLAB](https://img.shields.io/badge/MATLAB-Simulink-orange)](https://www.mathworks.com/products/simulink.html)
[![EKF](https://img.shields.io/badge/Estimation-Extended%20Kalman%20Filter-blueviolet)](https://en.wikipedia.org/wiki/Extended_Kalman_filter)
[![Nav](https://img.shields.io/badge/Use%20Case-Drone%20Navigation-green)](#)

</div>

---

## Overview

This repository contains the code and simulation results of an **Extended Kalman Filter (EKF)** designed for a **drone system** to estimate the **accelerometer bias**, which is added to the system as an **extra state**.

Built and tested with the **Navigation Toolbox** in **MATLAB / Simulink**, and validated against a **real dataset** collected at **Università di Pisa** during different drone flight experiments.

## What's Implemented

The work is organized as follows:

🔹 **EKF 1-DoF, EKF 3-DoF** — the EKF is first designed on a simplified **1- and 3-degree-of-freedom** system, with the dataset generated in the same simulation environment;

🔹 **EKF 3-DoF with feedback correction** — the EKF is tested and the estimated bias is fed back as a **correction input** to the attitude computation (**AHRS**);

🔹 **EKF with feedback correction on real dataset** — the filter is evaluated on the **real flight dataset** to assess performance.

For each test, a Simulink simulation was created and the proposed EKF was **compared against a standard Kalman filter without bias estimation**.

> ▶️ See **`contents.md`** inside each folder for file descriptions and how to run the Simulink simulations.

<div align="center">
  <img src="https://github.com/rachele182/navigation_systems/assets/75611841/b8a41c51-eb79-4c9b-b673-6eede590bbc5" width="365"/>
  <br/>
  <em>Schematic of the proposed filter, with the estimated accelerometer bias integrated as a correction into the attitude estimator (AHRS).</em>
</div>

## Key Skills & Tools

| Area             | What it demonstrates                          |
|------------------|-----------------------------------------------|
| State estimation | Extended Kalman Filter, bias estimation       |
| Navigation       | AHRS, attitude estimation, sensor correction  |
| Modeling         | 1-DoF → 3-DoF progressive filter design       |
| Validation       | Simulated **and** real drone flight dataset   |
| Tooling          | MATLAB / Simulink, Navigation Toolbox         |

## About

**Extended Kalman Filter** for **accelerometer bias estimation** on a drone, with **AHRS feedback correction** — designed in MATLAB/Simulink and validated on **real flight data** (Università di Pisa).
