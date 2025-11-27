# Cascaded-PID-based-EKF-for-Quadrotor-UAVs

This repository presents a lightweight cascaded PID-based Extended Kalman Filter that uses optical flow and IMU data for full quadrotor UAV stabilization. The controller is tuned in MATLAB and deployed as embedded code on a Teensy 4.1.

## MATLAB Tuning

The MATLAB folder contains the mathematical models for roll, pitch, yaw, and altitude dynamics together with the PID controllers. Users can use the MATLAB PID Tuner or manually adjust the controller gains to match their specific quadrotor system.

## Arduino Code

The Arduino implementation is provided in `pidekfquad.ino`. The code uses several libraries, including the SBUS library for the receiver input. The controller parameters and filtering settings follow the values obtained from the tuned MATLAB model.

## Hardware Overview

- Flight controller: Teensy 4.1  
- Sensors: IMU and optical flow sensor MTF 02P  
- Receiver: SBUS-compatible RC receiver  

## Usage

1. Tune the PID controllers in MATLAB for your specific quadrotor model.
2. Export or copy the tuned gains into `pidekfquad.ino`.
3. Compile and upload the Arduino code to the Teensy 4.1.
4. Power the quadrotor and verify stabilization performance in a safe test environment.


