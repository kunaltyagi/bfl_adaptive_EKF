# bfl_adaptive_EKF
Extension on Orocos BFL library to include a version of EKF which adapts Q and R matrices from an initial guess

## TODO
The algorithm requires complete change of internal computational algorithm, so either
* change the visibility of the struct in ExtendedKalmanFilter and KalmanFilter to protected, *OR*
* create it from scratch, copying everything from KalmanFilter and ExtendedKalmanFilter

Right now, I have adopted the second strategy.
I have yet to test it, so use it at your own peril :P However, I'm sure that it works, it just required certain optimizations.
