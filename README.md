### Simulation Results for ENPM667 Final Project

In this project, we have derived the equations of motion for a nonlinear system. The system in question is a cart with two hanging pendulums.
We have linearized the system around an equilibrium point. This linearized system will assist in designing the LQR/LQG controllers.

### Part C

In this section we explore the conditions in which the linearized system is controllable. This can be solved both symbolically or numerically.

### Part D

In this section we check that the parameters given are in fact controllable for the linearized system and then obtain an LQR controller. We Then tune the controller until it gives a reasonable output.
This tuned LQR controller is then applied to the nonlinear system. We also check for stability using Lyapunov's Indirect method.

### Part E

In this section we explore different ouput vector combinations, and determine which of them are observable on the linearized system

### Part F

In this section, we take the observable results from part E and obtain the best Lunenberger Observer for each of them. This observer is then applied to both the linearized and nonlinear system

### Part G

In this section, we choose the smallest observable output vector and combine the previously derived LQR controller with a Kalman Filter to obtain an LQG controller. 
