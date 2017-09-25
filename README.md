# Unscented Kalman Filter Project 
Self-Driving Car Engineer Nanodegree Program

---

This code was developed to track a vehicle using to simulated measurement inputs, LIDAR and Radar. The LIDAR unit provides x,y position of the object being tracked and the RADAR provides position, velocity and angle. The structure of this project could be used to track different systems provided the measurement model and system model equations are know. 

The Unscented Kalman Filter (UKF) provides the advantage of the the Extended Kalman Filter (EKF) in that it handles the non-linearities of the system better and the Jacobians of the system do not have to be calculated. The computation time of the UKF is slightly longer than the EKF but it is of the same order of magnitude. In almost all cases the UKF performs better at tacking a nonlinear system. 

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../dat
a/output2.txt` 

## Udacity Self Driving Car Original Project Code
The base code for the project was provided by Udacity and can be found [here](https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project)
The project has slightly changed since I did it to now run with the Udacity Term 2 Simulator.


