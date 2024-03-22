%p1

dt = 0.001;
x_k = 10;
X_k = [x_k];
U_k = [];

for k = 1:2000;
    u_k = x_k^3-cos(x_k)-10*x_k;
    x_k = (1-10*dt)*x_k;
    U_k = [U_k;u_k];
    X_k = [X_k;x_k];
end

figure;
subplot(2,1,1)
plot(0:0.001:1.999,X_k(1:2000))
title('Closed-Loop x(t)');
subplot(2,1,2)
plot(0:0.001:1.999,U_k)
xlabel('time(sec)');title('Control u(t)');
%%
%p2

A = [0 1 ; 0 0];
B = [0;1];
e = [-2;-3];
F = -place(A, B, e);
x_k = [1 ;-1];
X_k = [x_k];
U_k = [];

for k = 1:5000;
    x_1 = x_k(1);
    x_2 = x_k(2);
    u_k = x_1^2+(1/cos(x_2))*(F(1)*x_1+F(2)*sin(x_2));
    U_k = [U_k;u_k];
    x_k = dt*[sin(x_2);-x_1^2]+dt*[0;1]*u_k+x_k;
    X_k =[X_k,x_k];
end

figure;
subplot(2,1,1)
plot(0:0.001:5,X_k)
title('Closed-Loop x(t)');
subplot(2,1,2)
plot(0:0.001:4.999,U_k)
xlabel('time(sec)');title('Control u(t)');
%%
%p3

x_k = [1;-1];
X_k = [x_k];
U_k = [];

for k = 1:5000;
    x_1 = x_k(1);
    x_2 = x_k(2);
    u_k = x_1^2-2*x_2;
    U_k = [U_k;u_k];
    x_k = dt*[sin(x_2);-x_1^2]+dt*[0;1]*u_k+x_k;
    X_k =[X_k,x_k];
end

figure;
subplot(2,1,1)
plot(0:0.001:5,X_k)
title('Closed-Loop x(t)');
subplot(2,1,2)
plot(0:0.001:4.999,U_k)
xlabel('time(sec)');title('Control u(t)');
%%
%p4

A = [0 1 ; 0 0];
B = [0;1];
e = [-2;-3];
F = -place(A, B, e);
x_k = [1;-1];
X_k = [x_k];
U_k = [];

for k =1:5000;
    x_1 = x_k(1);
    x_2 = x_k(2);
    u_k = x_1-2*(1-x_1^2)*x_2+F(1)*x_1+F(2)*x_2;
    U_k = [U_k;u_k];
    x_k = dt*[x_2;-x_1+2*(1-x_1^2)*x_2]+dt*[0;1]*u_k+x_k;
    X_k =[X_k,x_k]; 
end

figure;
subplot(2,1,1)
plot(0:0.001:5,X_k)
title('Closed-Loop x(t)');
subplot(2,1,2)
plot(0:0.001:4.999,U_k)
xlabel('time(sec)');title('Control u(t)');