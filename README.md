# series-elastic-actuation-impedance-control

This repo addresses impedance controller design problem of a class of SEA cascaded control architectures, which is composed of outer-impedance and inner-torque feedback loops. Our study proposes gain design criterion solving optimal controller gains by maximizing phase-margin-based stability. Meanwhile, we observe a trade-off between impedance and torque controller gains and analyze their interdependence in terms of closed-loop stability and overall impedance performance. With the proposed controller design criterion, we adopt frequency-domain methods to thoroughly analyze the effects of time delays, filtering and load inertia on SEA impedance performance. A novel impedance performance metric, defined as "Z-region", is proposed to simultaneously quantify achievable impedance magnitude range (i.e., Z-width) and frequency range (i.e., Z-depth). Maximizing the Z-region enables SEA-equipped robots to achieve a wide variety of impedance tasks without alternating the control structure.

This is a MATLAB algorithm implementation of the demos in the following paper:

Y. Zhao, S. Jorgensen, N. Paine, L. Sentis, [Impedance Control and Performance Measure of Series Elastic Actuators](http://sites.utexas.edu/hcrl/files/2016/01/08016601.pdf), IEEE Transactions on Industrial Electronics, 2018, 65(3), 2817-2827.

A video of experimental validation is attached [here](https://youtu.be/biIdlcAMPyE) 

Please feel free to try the demo scripts. Your open issues/questions/comments are welcome (yezhao@utexas.edu, lsentis@austin.utexas.edu).