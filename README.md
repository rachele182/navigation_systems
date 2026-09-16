<div align="center">

# Underwater Systems — AUV Simulation, Control & Navigation 🌊

**Complete modeling, control and Kalman-filter navigation of an Autonomous Underwater Vehicle (AUV) performing a lawn-mower survey in an unknown area**

**Authors:** Rachele Nebbia Colomba · Chiara Sammarco · Francesco Vezzi · Matteo Paiano *

*Final project — MSc course "Underwater Systems" · Robotics and Automation Engineering, Università di Pisa*

[![MATLAB](https://img.shields.io/badge/MATLAB-Simulink-orange)](https://www.mathworks.com/products/simulink.html)
[![Kalman](https://img.shields.io/badge/State%20Estimation-Kalman%20Filter-blueviolet)](https://en.wikipedia.org/wiki/Kalman_filter)
[![Robotics](https://img.shields.io/badge/Robotics-AUV%20Control-green)](https://en.wikipedia.org/wiki/Autonomous_underwater_vehicle)

</div>

---

## Overview

Final project for the **Underwater Systems** course (MSc Robotics and Automation Engineering, **Università di Pisa**).

The work delivers the **complete model of an AUV** — environment, sensors, controller and navigation filter — required to autonomously run a **lawn-mower survey** over an **unknown** seabed. Everything is built and validated in **MATLAB / Simulink**.

The integrated system is split into **five models**:

🔸 **Trajectory Generator** — computes the survey trajectory over the unknown area;  
🔸 **Vehicle Model** — AUV geometric parameters, dynamics and **thruster** positions;  
🔸 **Control** — **PID controllers** for the AUV thrusters;  
🔸 **Sensor + Environment** — full seabed model and the chosen sensor models;  
🔸 **Navigation** — the **Kalman navigation filter**.

<div align="center">
<img src="https://github.com/rachele182/navigation_systems/assets/75611841/16b22289-f5a4-4cf3-a26e-ecd3426b7a5f" width="375">
  <br/>
  <em>The modeled AUV, "Pasqualo".</em>
</div>

## Repository Structure

```
.
├── sensor_model/   # MATLAB/Simulink models of the underwater environment and sensors
├── mission/        # Simulink model + MATLAB scripts to run the full mission
├── animation/      # Script to animate the executed mission
└── README.md
```

> ▶️ See **`mission/contents.md`** for a description of each file and instructions to run the simulation and the animation.

## What's Inside

This project shows a complete, end-to-end AUV pipeline rather than a single block:

- **Trajectory planning** for area-coverage survey missions;
- **Hydrodynamic vehicle modeling** incl. thruster configuration;
- **PID-based thruster control**;
- **Sensor + environment simulation** (seabed reconstruction);
- **Kalman-filter state estimation / sensor fusion** for navigation.

## About This Repo (info note)

This project was built as **part of a team effort** with multiple students. For copyright reasons, the detailed scripts are disclosed in full only in **`sensor_model`** (the sensor/environment models developed by the authors marked with *). The `mission` and `animation` folders contain the **final integrated mission**, which is the result of merging all five modules.

<div align="center">
  <img src="https://github.com/rachele182/navigation_systems/assets/75611841/39082569-4841-47a7-8545-c70805ac7949" width="425">
</div>

## Key Skills & Tools

| Area                | What it demonstrates                              |
|---------------------|---------------------------------------------------|
| State estimation    | Kalman filter, sensor fusion, navigation          |
| Control             | PID thrust allocation, trajectory tracking        |
| Modeling            | Hydrodynamics, thruster configuration, seabed env |
| Simulation          | MATLAB / Simulink, full mission + animation       |
| Robotics            | AUV/underwater robotics, autonomous navigation    |

## About

Complete **simulation, control and navigation-filter design** for an AUV performing a lawn-mower (**area-coverage**) survey in an **unknown** area — built in MATLAB/Simulink, with **Kalman-filter** state estimation.



