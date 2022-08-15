clear;clc;close all
T = 0:0.1:4-0.1;
T = T';
x = 15*ones(length(T),1);
y1 = linspace(-10,10,length(T)/2)';
y2 = linspace(10,-10,length(T)/2)';
y = [y1 ; y2];
z = 15*ones(length(T),1);
tta1 = zeros(length(T),1);
tta2 = zeros(length(T),1);
tta3 = zeros(length(T),1);
for i = 1:length(T)
    [tta1(i),tta2(i),tta3(i)] = inv_kin(x(i),y(i),z(i));
    
end
t1.time=[ ];
t1.signals.values=tta1;
t1.signals.dimensions=1;
t2.time=[ ];
t2.signals.values=tta2;
t2.signals.dimensions=1;
t3.time=[ ];
t3.signals.values=tta1;
t3.signals.dimensions=1;