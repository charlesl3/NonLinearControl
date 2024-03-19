% p1
Ac = [0 1; -100 0];
Bc = [0;1];
C = [-100  0];
D = 1;
dt = 0.01;
%%
% p2
A =  eye(2) + Ac*dt + 1/factorial(2)*(Ac*dt)^2 + 1/factorial(3)*(Ac*dt)^3 + 1/factorial(4)*(Ac*dt)^4+1/factorial(5)*(Ac*dt)^5
B = (eye(2)*dt + 1/factorial(2)*Ac*dt^2 + 1/factorial(3)*Ac^2*dt^3 + 1/factorial(4)*Ac^3*dt^4+1/factorial(5)*Ac^4*dt^5)*Bc
C = C
D = D
system = ss(Ac,Bc,C,D);
sysd = c2d(system,dt,'zoh');
[a,b,c,d] = ssdata(sysd)

u = zeros(500,1); 
u(1) = 1;% u(0) = 1
y = [];

x_k = [0;0]; % x(0)
for i=1:500
    y_k = c*x_k+d*u(i);
    x_k = a*x_k+b*u(i);
    y = [y;y_k];
end

figure(1)
plot([1:500]*0.01,y)
xlabel('time(sec)');title('Unit Pulse Response (Output)');

%%
% p3
t2 =0:0.01:5;
u = zeros(1, length(t2));   
X0_lqr = [1;-1];

Q1 = 10*eye(2);
R1 = 1;
F1 = dlqr(a,b,Q1,R1);
y2 = [];
x_k = [1;-1]; % x(0)
x = x_k;
u2 = [];
J_lrq = 0;
for i=1:500
    u = -F1*x_k;
    J_lrq = J_lrq + x_k'*Q1*x_k + u'*R1*u;
    y_k = (c-d*F1)*x_k;
    x_k = (a-b*F1)*x_k;
    y2 = [y2,y_k];
    u2 = [u2;u];
    x = [x,x_k];
end
figure(2)
subplot(3,1,1)
plot([1:500]*0.01,x(1,1:500)')
title('LQR Closed-loop Position');
subplot(3,1,2)
plot([1:500]*0.01,x(2,1:500)')
title('LQR Closed-loop Velocity');
subplot(3,1,3)
plot([1:500]*0.01,u2)
xlabel('time(sec)');title('LQR Closed-loop Control Input');

cost = [];
 
for j = [0.95:0.01:0.99];
    x_k = [1;-1]; % x(0) 
    e = eig(a)*j;
    F = place(a, b, e);
    J = 0;
    for i=1:500
        u = -F*x_k;
        J = J + x_k'*Q1*x_k + u'*R1*u;
        x_k = (a-b*F)*x_k;
    end
    cost = [cost;J];
end

figure(3)
hold on
plot(cost,'-*');
yline(J_lrq);
xlabel('poles');title('cost of different controllers');
legend('cost by other controllers','cost by LQR controller')
hold off 
%%
%p4

randn('seed',1);
u = 10*randn(500,1);  
x_k = [1;-1]; % x(0)
y = [];
x = [];

for i=1:500
    y_k = c*x_k+d*u(i);
    x_k = a*x_k+b*u(i);
    y = [y,y_k];
    x = [x,x_k];
end
figure(4)
subplot(4,1,1)
plot([1:500]*0.01,x(1,1:500)')
title('Position (State 1)');
subplot(4,1,2)
plot([1:500]*0.01,x(2,1:500)')
title('Velocity (State 2)');
subplot(4,1,3)
plot(y)
title('Acceleration (output)');
subplot(4,1,4)
plot(u)
xlabel('time(sec)');title('Random input excitation');


%lsim
x_k = [1;-1]; % x(0)
[y,x] = dlsim(a,b,c,d,u,x_k);
figure(5)
subplot(4,1,1)
plot([1:500]*0.01,x(:,1))
title('simulated Position (State 1)');
subplot(4,1,2)
plot([1:500]*0.01,x(:,2))
title('simulated Velocity (State 2)');
subplot(4,1,3)
plot(y)
title('simulated Acceleration (output)');
subplot(4,1,4)
plot(u)
xlabel('time(sec)');title('Random input excitation');
%%
%p5
e = eig(a)*0.95;
M = place(a', c', e);
M = -M';
x_ob = [];
x_k = [0;0];

for i=1:500
    x_k = (a+M*c)*x_k+(b+M*d)*u(i)-M*y(i);
    x_ob = [x_ob,x_k];
end

figure(6)
subplot(2,1,1)
hold on
plot([1:500]*0.01,x_ob(1,:))
plot([1:500]*0.01,x(:,1))
title('True vs. Estimated Position');
hold off

subplot(2,1,2)
hold on
plot([1:500]*0.01,x_ob(2,:))
plot([1:500]*0.01,x(:,2))
xlabel('time(sec)');title('True vs. Estimated Velocity');
hold off
%%
%p6
F = -F;
M = M;

x_k = [1;-1;0;0];
A_cl = [a, b*F; -M*c,a+M*c+b*F];
x_cl = [];
u_cl = [];
for i=1:500
    u_cl = [u_cl;F*x_k(1:2)];
    x_k = A_cl*x_k;
    x_cl = [x_cl,x_k];
end

figure(7)
subplot(3,1,1)
plot([1:500]*0.01,x_cl(1,:))
title('Observer-Based LQR Closed-Loop Position');
subplot(3,1,2)
plot([1:500]*0.01,x_cl(2,:))
title('Observer-Based LQR Closed-Loop Velocity');
subplot(3,1,3)
plot(u_cl)
xlabel('time(sec)');title('Observer-Based LQR Closed-Loop Control Input');