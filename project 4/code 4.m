%p1
warning('off')
for k=1:50, ystar(k,1)=1-cos(pi*k/50);end
ydes=[zeros(50,1);ystar;2*ones(100,1);ystar(50:-1:1);zeros(50,1)];


M = eye(5);
Cd_real = [8 -4 0 0 0; -4 8 -4 0 0; 0 -4 8 -4 0 ;0 0 -4 8 -4;0 0 0 -4 8];
K_real = [240 -120 0 0 0; -120 240 -120 0 0; 0 -120 240 -120 0 ;0 0 -120 240 -120;0 0 0 -120 240];
Bf_real = [1 0; 0 0;0 1;0 0;0 0];

Cd = 0.5*[8 -4 0 0 0; -4 8 -4 0 0; 0 -4 8 -4 0 ;0 0 -4 8 -4;0 0 0 -4 8];
K = [240 -120 0 0 0; -120 240 -120 0 0; 0 -120 240 -120 0 ;0 0 -120 240 -120;0 0 0 -120 240]/1.2;
Bf = [1 0; 0 0;0 1;0 0;0 0];

% step 2
Ac_real = [zeros(5),eye(5);-inv(M)*K_real,-inv(M)*Cd_real]
Bc_real = [zeros(5,2);inv(M)*Bf_real]
C_real = [0 0 0 0 1 0 0 0 0 0]
D_real = [0 0]
dt = 0.01;
system = ss(Ac_real,Bc_real,C_real,D_real);
sysd = c2d(system,dt,'zoh');
[a_real,b_real,c_real,d_real] = ssdata(sysd);

Ac = [zeros(5),eye(5);-inv(M)*K,-inv(M)*Cd]
Bc = [zeros(5,2);inv(M)*Bf]
C = [0 0 0 0 1 0 0 0 0 0]
D = [0 0]
dt = 0.01;
system = ss(Ac,Bc,C,D);
sysd = c2d(system,dt,'zoh');
[a,b,c,d] = ssdata(sysd);

r =[];
for p = 0:299;
    r = [r;c*a^p*b(:,1)];
end

P = tril(toeplitz(r));
Q = (10^8)*eye(300);
S = (10^2)*eye(300);
alpha = 0.7;
L = alpha*inv(P'*Q*P+S)*P'*Q;


% no disturbances
error = [];
u0 = zeros(300,1);
x0 = zeros(10,1);
Y = [];
for i  = 1:100;
    Yi = [];
    for k = 0:299;
        x0 = a_real*x0 + b_real*[u0(k+1);0]; %x(1)
        y = c_real*x0; %y(1)
        Yi = [Yi;y];
    end
    e = ydes - Yi;
    u0 = L*e+u0;
    Y = [Y,Yi];
    error = [error;norm(e,2)];
end
  

figure;
subplot(2,1,1)
plot(0.01:0.01:3,Y);ylim([-0.5 2.5])
xlabel('time(sec)');ylabel('Output');title('Output Trajectories at Different Learning Iterations');
subplot(2,1,2)
plot(0:0.01:2.999,zeros(300,1));ylim([-1 1])
xlabel('time(sec)');ylabel('Dist Input');title('Disturbance Input');

figure;
subplot(2,1,1)
hold on
plot(0.01:0.01:2.999,Yi(2:end));ylim([-0.5 2.5])
plot(0:0.01:2.999,ydes)
hold off
xlabel('time(sec)');ylabel('Desired and Output Trajectories');title('Desired Output Trajectory vs. Output Trajectory at Last Iteration');
subplot(2,1,2)
plot(0:0.01:2.999,u0);ylim([-500 2500])
xlabel('time(sec)');ylabel('Control Input');title('Control Input at Last Iteration');

figure;
semilogy(error,"-*");ylim([10^(-1) 10^2])
xlabel('Iteration Number');ylabel('Norm of Output Tracking Error');title('Norm of Output Tracking Error vs. Iteration Number');



% with disturbances
dist = @(t) -1000*exp(-t)*sin(2*t);
error = [];
u0 = zeros(300,1);
x0 = zeros(10,1);
Y = [];
for i  = 1:100;
    Yi = [];
    for k = 0:299;
        x0 = a_real*x0 + b_real*[u0(k+1);-1000*exp(-k*dt)*sin(2*k*dt)]; %x(1)
        y = c_real*x0; %y(1)
        Yi = [Yi;y];
    end
    e = ydes - Yi;
    u0 = L*e+u0;
    Y = [Y,Yi];
    error = [error;norm(e,2)];
end


figure;
subplot(2,1,1)
plot(Y(2:end,:));ylim([-4 4])
xlabel('time(sec)');ylabel('Output');title('Output Trajectories at Different Learning Iterations');
subplot(2,1,2)
fplot(dist,[0 3]);ylim([-600 200])
xlabel('time(sec)');ylabel('Disturbance Input');title('Disturbance Input');

figure;
subplot(2,1,1)
hold on
plot(0.01:0.01:2.999,Yi(2:end));ylim([-0.5 2.5])
plot(0:0.01:2.999,ydes)
hold off
xlabel('time(sec)');ylabel('Desired and Output Trajectories');title('Desired Output Trajectory vs. Output Trajectory at Last Iteration');
subplot(2,1,2)
plot(0:0.01:2.999,u0);ylim([-500 2500])
xlabel('time(sec)');ylabel('Control Input');title('Control Input at Last Iteration');

figure;
semilogy(error,"-*");ylim([10^(-1) 10^2])
xlabel('Iteration Number');ylabel('Norm of Output Tracking Error');title('Norm of Output Tracking Error vs. Iteration Number');




%%
%p2
format short e
n = 10;
randn('seed',1), A = randn(2*n,n)
randn('seed',2), b = randn(2*n,1)
format long e
XLS = pinv(A)*b

mu_init = zeros(10,1);
var_init = 100*eye(10);
J = @(A,x,b) (A*x-b)'*(A*x-b);

for i = 1:2000;
    pop = mvnrnd(mu_init,var_init,1000)';
    CE = [];
    for c = 1:size(pop,2);
        z = pop(:,c);
        CE = [CE,J(A,z,b)];
    end
    [~, idx] = sort(CE, 'ascend');

    pop_sorted = pop(:,idx);
    pop = pop_sorted(:,1:400);

%     C = zeros(size(pop,1));
%     for c = 1:size(pop,2);
%         z = pop(:,c);
%         C = C+(((z-mean(z))*(z-mean(z))')./length(z));
%     end
%     
%     mu_init = mean(pop,2);
%     var_init = C./size(pop,2);
    mu_init = mean(pop,2);
    var_init = cov(pop');

    XCE = pop(:,1);

    if sum(abs(XCE-XLS)<=7e-10) == 10;
        break;
    end   

end

XCE

XLS-XCE