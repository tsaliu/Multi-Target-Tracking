clc; clear; close all;

[data]=textread('rmsev1.txt','%f');
data=reshape(data,[2,100])';

[data2]=textread('rmsev2.txt','%f');
data2=reshape(data2,[2,100])';

[data3]=textread('rmsev3.txt','%f');
data3=reshape(data3,[2,100])';

k=0:1:length(data)-1;


figure(1) %300 400
subplot(2,1,1)
plot(k,data(:,1)), grid;
xlabel('Time(k)'); title('RMSE v_x 100 Runs');
subplot(2,1,2)
plot(k,data(:,2)), grid;
xlabel('Time(k)'); title('RMSE v_y 100 Runs');

figure(2) %700 300
subplot(2,1,1)
plot(k,data2(:,1)), grid;
xlabel('Time(k)'); title('RMSE v_x 100 Runs');
subplot(2,1,2)
plot(k,data2(:,2)), grid;
xlabel('Time(k)'); title('RMSE v_y 100 Runs');
figure(3) %700 700
subplot(2,1,1)
plot(k,data3(:,1)), grid;
xlabel('Time(k)'); title('RMSE v_x 100 Runs');
subplot(2,1,2)
plot(k,data3(:,2)), grid;
xlabel('Time(k)'); title('RMSE v_y 100 Runs');
