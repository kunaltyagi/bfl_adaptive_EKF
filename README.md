# bfl_adaptive_EKF
Extension on Orocos BFL library to include a version of EKF which adapts Q and R matrices from an initial guess

## TODO
The algorithm requires complete change of internal computational algorithm, so either
* change the visibility of the struct in ExtendedKalmanFilter and KalmanFilter to protected, *OR*
* create it from srratch, copying everything from KalmanFilter and ExtendedKalmanFilter
