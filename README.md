# C program for solving the second order ordinary differential equation with a linear boundary condition
The second order ordinary differential equation is written in form $y''+p(x)y'+q(x)y=f(x)$

The boundary condition:

$$\begin{cases}
  \alpha_{0}y(x_0)+\alpha_{1}y'(x_0)=A\\
  \beta_{0}y(x_n)+\beta_{1}y'(x_n)=B
\end{cases}$$

## Finite difference method for solving the boundary problem
1. Partition the interval $\[x_0,x_n\]$ into subintervals of length $h$
2. Approximate the derivatives using finite differences
3. Solve the system of linear equations:

$$\begin{cases}
   \frac{y_{i+1}-2y_{i}+y_{i-1}}{h^2}+p_{i}\frac{y_{i+1}-y_{i-1}}{2h}+q_{i}y_{i}=f_{i}\\
   \alpha_{0}y(x_0)+\alpha_{1}y'(x_0)=A & \quad(1)\\
   \beta_{0}y(x_n)+\beta_{1}y'(x_n)=B
\end{cases}$$

## Double sweep method for solving the system (1)
1. Write the system in the tri-diagonal form
2. Calculate double sweep coefficients
3. Find unknown $y_{i}$ in the reverse order
