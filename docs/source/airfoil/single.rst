**************************
Running single analysis
**************************

It is highly recommended to go through sample generation process with CST or FFD before proceeding ahead.
Typically, surrogate-based optimization methods involve creating a surrogate model and finding a new infill point
based on some criteria. Once infill point is found, expensive analysis needs to be performed to evaluate the 
objective function value.

Blackbox provides a method to run these expensive single analysis using the surrogate model. This is useful when you want to run
a single analysis and do not want to go through the process of creating a sample and surrogate model.
