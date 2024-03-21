 % p1
warning('off')
xt = @(k) 30 + 4*exp(k)*cos(10*k);
vxt = @(k) 4*exp(k)*cos(10*k)-40*exp(k)*sin(10*k);
ht = @(k) 10+exp(2*k);
vht = @(k) 2*exp(2*k);

k = 0:1:1999;

figure;
subplot(2,1,1)
fplot(xt,[0 2])
title('True x-position');ylim([0 60])
subplot(2,1,2)
fplot(ht,[0 2])
title('True h-position');ylim([0 100])
%%
figure;
subplot(2,1,1)
fplot(vxt,[0 2])
title('True x-velocity');ylim([-500 500])
subplot(2,1,2)
fplot(vht,[0 2])
xlabel('time(sec)');title('True h-velocity');ylim([0 150])
%%
%p2 and p3
randn('seed',1);v_r = randn(2000,1);
randn('seed',2);v_theta = 0.01*randn(2000,1);
X_naive = [];
H_naive = [];
X_noise = [];
H_noise = [];
r_noise = [];
theta_noise = [];
ind = 0;

for k = 0:0.001:1.999;
    ind = ind+1;
    x_k = 30 + 4*exp(k)*cos(10*k);
    h_k = 10+exp(2*k);
    r_k = sqrt(x_k^2+h_k^2);
    theta_k = atan(h_k/x_k);
    r_k_noise = sqrt(x_k^2+h_k^2)+v_r(ind);
    theta_k_noise = atan(h_k/x_k)+v_theta(ind);

    x_true_k = r_k*cos(theta_k);
    h_true_k = r_k*sin(theta_k);
    x_noise_k = r_k_noise*cos(theta_k_noise);
    h_noise_k = r_k_noise*sin(theta_k_noise);

    r_noise =[r_noise;r_k_noise];
    theta_noise = [theta_noise;theta_k_noise];
    X_naive = [X_naive;x_true_k];
    H_naive = [H_naive;h_true_k];
    X_noise = [X_noise;x_noise_k];
    H_noise = [H_noise;h_noise_k];
end

VX_naive = [];
VH_naive = [];
VX_noise = [];
VH_noise = [];
for k = 1:1999;
    VX_naive = [VX_naive;(X_naive(k+1)-X_naive(k))/0.001];
    VH_naive = [VH_naive;(H_naive(k+1)-H_naive(k))/0.001];
    VX_noise = [VX_noise;(X_noise(k+1)-X_noise(k))/0.001];
    VH_noise = [VH_noise;(H_noise(k+1)-H_noise(k))/0.001];
end

figure;
plot(X_naive,H_naive);ylim([0 70]);
xlabel('x');ylabel('h');title('True flight path of the fly');

figure;
subplot(2,1,1)
plot(0:0.001:1.999,r_noise);ylim([0 100]);
title('r w/ noise');
subplot(2,1,2)
plot(0:0.001:1.999,theta_noise)
xlabel('time(sec)');title('theta w/ noise');


figure;
subplot(2,1,1)
hold on
fplot(xt,[0 2])
plot(0:0.001:1.999,X_noise);
title('Estimated vs. True x-position by Naive Method');
hold off
subplot(2,1,2)
hold on
fplot(ht,[0 2])
plot(0:0.001:1.999,H_noise);
hold off
xlabel('time(sec)');title('Estimated vs. True h-position by Naive Method');


figure;
subplot(2,1,1)
hold on
fplot(vxt,[0 2])
plot(0:0.001:1.998,VX_noise)
title('Estimated vs. True x-velocity by Naive Method');
hold off
subplot(2,1,2)
hold on
fplot(vht,[0 2])
plot(0:0.001:1.998,VH_noise)
hold off
xlabel('time(sec)');title('Estimated vs. True h-velocity by Naive Method');

figure;
hold on
plot(X_naive,H_naive);
plot(X_noise,H_noise,'.');ylim([0 70]);
hold off
xlabel('x');ylabel('h');title('Estimated flight path (red dots) vs. truth (blue) by naive method');
legend('truth','estimated')

%%
% p5 
dt = 0.001;
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
G = [0 0;dt 0;0 0;0 dt];
z_kk = [1;0;0;0];
P_kk = 10^-2*eye(4);
Q = G*(10^7*eye(2))*G';
R = [1 0; 0 10^-4];
EKF_Z = [z_kk];

for k = 1:1999;
    z_k1k = A*z_kk;
    x_k = z_k1k(1);
    h_k = z_k1k(3);
    C_k1 = [x_k/(sqrt(h_k^2+x_k^2)) 0 h_k/(sqrt(h_k^2+x_k^2)) 0; -h_k/(h_k^2+x_k^2) 0 x_k/(h_k^2+x_k^2) 0];

    P_k1k = A*P_kk*A'+Q;
    Lk1 = P_k1k*C_k1'*pinv(R+C_k1*P_k1k*C_k1');
    P_kk = (eye(4)-Lk1*C_k1)*P_k1k;
    y_measure = [r_noise(k+1);theta_noise(k+1)];
    g_k1 = [sqrt(h_k^2+x_k^2);atan(h_k/x_k)];
    z_kk = z_k1k + Lk1*(y_measure - g_k1);
    EKF_Z = [EKF_Z,z_kk];
end

figure;
subplot(2,1,1)
hold on
fplot(xt,[0 2])
plot(0:0.001:1.999,EKF_Z(1,:))
title('Estimated vs. True x-position by EKF');
hold off
subplot(2,1,2)
hold on
fplot(ht,[0 2])
plot(0:0.001:1.999,EKF_Z(3,:))
hold off
xlabel('time(sec)');title('Estimated vs. True y-position by EKF');
legend('truth','estimated')

figure;
subplot(2,1,1)
hold on
fplot(vxt,[0 2])
plot(0:0.001:1.999,EKF_Z(2,:))
title('Estimated vs. True x-velocity by EKF');
hold off
subplot(2,1,2)
hold on
fplot(vht,[0 2])
plot(0:0.001:1.999,EKF_Z(4,:));ylim([-100 400])
hold off
xlabel('time(sec)');title('Estimated vs. True y-velocity by EKF');
legend('truth','estimated')

figure;
hold on
plot(X_naive,H_naive);
plot(EKF_Z(1,:),EKF_Z(3,:),'.');ylim([0 70]);
hold off
xlabel('x');ylabel('h');title('Estimated flight path (red) vs. truth (blue) by EKF');
legend('truth','estimated')
