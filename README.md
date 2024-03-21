# Post-Modern and Non-Linear Control

This is the repo for course projects of ENGG 149 at Dartmouth College.

### For system implementation

```Matlab
system = ss(Ac,Bc,C,D);
sysd = c2d(system,dt,'zoh');
[a,b,c,d] = ssdata(sysd)
```

### For pole placement

```Matlab
e = eig(a)*0.95;
M = place(a', c', e);
```

### for system simulation

```Matlab
F1 = dlqr(a,b,Q1,R1);
```

### for LQR controller

```Matlab
[y,x] = dlsim(a,b,c,d,u,x_k);
```

## Project 1: Review of Modern Control Theory
## Project 2:The Extended Kalman Filter (EKF)


