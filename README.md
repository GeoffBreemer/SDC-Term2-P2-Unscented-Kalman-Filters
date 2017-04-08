# Project 2: Unscented Kalman Filters
This document describes the submission for **Project 2: Unscented Kalman Filters**. Boilerplate code provided by Udacity was improved by:

1. removing repetitive code

2. detecting zeros for certain values and setting them to a small value instead, e.g. to avoid dividing by zero when preparing the Jacobian

3. encapsulating code instead of having public class variables, e.g. as keen be seen in the `KalmanFilter` class.

The plots shown below were generated using the `ekf-visualization.ipynb` notebook, which was provided by Udacity.

## Running the Project
Perform the `Basic Build Instructions` provided by Udacity. Then run the project from the `build` directory on the two data files using:

`./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt ../output/output1.txt`

and

`./UnscentedKF ../data/sample-laser-radar-measurement-data-2.txt ../output/output2.txt`

## Results
The results for input file `sample-laser-radar-measurement-data-1.txt` are:


The output is available in file `output\output1.txt`. The results are visualised in the plot below:

![image1]

The results for input file `sample-laser-radar-measurement-data-2.txt` are:

The output is available in file `output\output2.txt`. The results are visualised in the plot below:

![image2]

[//]: # (Image References)

[image1]: output/output1.png "Output1.txt visualisation"
[image2]: output/output2.png "Output2.txt visualisation"