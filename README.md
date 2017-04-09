# Project 2: Unscented Kalman Filters
This document describes the submission for **Project 2: Unscented Kalman Filters**. Boilerplate code provided by Udacity was improved by removing repetitive code, removing unused variables, reusing P1 code and updating `main.cpp` to ensure it works with the visualization notebook. 

The plots shown below were generated using the `ukf-visualization-extended.ipynb` notebook, which was provided by Udacity though did require minor changes.

## Running the Project
Perform the `Basic Build Instructions` provided by Udacity. Then run the project from the `build` directory on the two data files using:

`./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt ../output/output1.txt`

and

`./UnscentedKF ../data/sample-laser-radar-measurement-data-2.txt ../output/output2.txt`

## Results (using both radar and laser)
The results for input file `sample-laser-radar-measurement-data-1.txt` are:

~~~~
Accuracy - RMSE:
0.0672946
0.0756271
0.558942
0.569706
Done processing ../data/sample-laser-radar-measurement-data-1.txt
~~~~

The output is available in file `output\output1.txt`. The results are visualised in the plot below:

![image1]

The results for input file `sample-laser-radar-measurement-data-2.txt` are:

~~~~
Accuracy - RMSE:
0.181708
0.180736
0.261946
0.283641
Done processing ../data/sample-laser-radar-measurement-data-2.txt
~~~~

The output is available in file `output\output2.txt`. The results are visualised in the plot below:

![image2]

The plot below shows the ground truth versus the estimates velocity:

![image3]

## Radar only
The results for `sample-laser-radar-measurement-data-1.txt` using only the **radar** measurements are shown below:

~~~~
Accuracy - RMSE:
0.128342
0.14145
0.62974
0.646617
Done processing ../data/sample-laser-radar-measurement-data-1.txt
~~~~

The results for `sample-laser-radar-measurement-data-2.txt` using only the **radar** measurements are shown below:

~~~~
Accuracy - RMSE:
0.249143
0.677432
0.820647
0.421344
Done processing ../data/sample-laser-radar-measurement-data-2.txt
~~~~

For both data files the accuracy drops when only radar measurements are used, as expected. 

## Laser only
The results for `sample-laser-radar-measurement-data-1.txt` using only the **laser** measurements are shown below:

~~~~
Accuracy - RMSE:
0.0932132
0.0813661
0.653073
0.614186
Done processing ../data/sample-laser-radar-measurement-data-1.txt
~~~~

The results for `sample-laser-radar-measurement-data-2.txt` using only the **laser** measurements are shown below:

~~~~
Accuracy - RMSE:
0.207084
0.187352
0.363471
0.303702
Done processing ../data/sample-laser-radar-measurement-data-2.txt
~~~~

Again, for both data files the accuracy drops when only laser measurements are used, as expected. However, for both files the **position accuracy** is still very similar to when using both measurement types, and better than when using radar measurements only.

As expected, **velocity accuracy** is not affected too badly when only using laser measurements compared to when both measurement types are used. Velocity accuracy for the *first* data set is only slightly worse compared to using radar measurements only. Interestingly, velocity accuracy for the *second* data set *improves* a lot when using laser measurements only. This may be caused by the small number of radar measurements availabe in the second data set.

[//]: # (Image References)

[image1]: output/output1.png "Output1.txt visualisation"
[image2]: output/output2.png "Output2.txt visualisation"
[image3]: output/output2_v.png "Output2.txt ground truth and estimated velocity"