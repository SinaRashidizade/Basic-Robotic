clear; clc; close all;
%% computing transformation matrix

syms tta1(t) tta2(t) tta3(t);

tta1dot = diff(tta1(t),t,1);
tta2dot = diff(tta2(t),t,1);
tta3dot = diff(tta3(t),t,1);

tta1dd = diff(tta1dot,t,1);
tta2dd = diff(tta2dot,t,1);
tta3dd = diff(tta3dot,t,1);

q = [tta1(t); tta2(t); tta3(t)];
qdot = [tta1dot; tta2dot; tta3dot];
qdd = [tta1dd; tta2dd; tta3dd];

PI = sym(pi);
T1 = compute_DH_mod_tf(0,0,335,tta1(t));
T2 = compute_DH_mod_tf(75,-PI/2,0,tta2(t));
T3 = compute_DH_mod_tf(270,0,0,tta3(t)-PI);

T = cat(3,T1,T2,T3);

T_ee = simplify(T1*T2*T3);
%% Some of Equations

T_eedot = diff(T_ee,t);
P_eedot = T_eedot(1:3,4);

R_ee = T_ee(1:3,1:3);    
R_eedot = T_eedot(1:3,1:3);
S = R_eedot*transpose(R_ee);

T0 = sym('T0',[4,4,3]);
T0(:,:,1) = T(:,:,1);
for n=2:3
    T0(:,:,n) = T0(:,:,n-1) * T(:,:,n);
end

R = sym('R',[3,3,3]);
p = sym('p',[3,1,3]);
for n=1:3
    R(:,:,n) = T(1:3,1:3,n);
    p(:,:,n) = T(1:3,4,n);
end

w = sym('w',[3,1,3]);
w0 = [0; 0; 0];
w(:,:,1) = transpose(R(:,:,1))*w0 + [0; 0; qdot(1)];
for n=2:3
    w(:,:,n) = transpose(R(:,:,n))*w(:,:,n-1) + [0; 0; qdot(n)];
end

v = sym('v',[3,1,3]);
v0 = [0; 0; 0];
v(:,:,1) = transpose(R(:,:,1)) * (v0 + cross(w0,p(:,:,1)));
for n=2:3
    v(:,:,n) = transpose(R(:,:,n)) * (v(:,:,n-1) + cross(w(:,:,n-1),p(:,:,n)));
end
%% Calculate Dynamic with Lagrange method

% Inertia Matrix
I = zeros(3,3,3);
I(:,:,1) = [384957.06 -2691.73 43626.55; -2691.73 180695.89 -19024.35; 43626.55 -19024.35 309313.62];
I(:,:,2) = [ 5886595.04 330234.57 368561.17; 330234.57 4128626.67 2824947.69; 368561.17 2824947.69 2826241.18];
I(:,:,3) = [179937.28 34003.18 -94339.01; 34003.18 207856.24 -34453.66; -94339.01 -34453.66 157860.07];

% center of mass from solid
cm = zeros(3,1,3);
cm(:,:,1) = [-2.47;32.05;-9.31];
cm(:,:,2) = [16.87;116.62;144.30];
cm(:,:,3) = [43.73;21.06;-45.00]; 

% Mass parametrs
m = zeros(3,1);
m(1) = 168.13;
m(2) = 167.87;
m(3) = 45.54;

R0 = sym('R0',[3,3,3]);
p0 = sym('p0',[3,1,3]);
zc = sym('zc',[3,1,3]);
for n=1:3
    R0(:,:,n) = T0(1:3,1:3,n);
    p0(:,:,n) = T0(1:3,4,n);
    zc(:,:,n) = T0(1:3,3,n);    % T_cito0 = T_citoi * T_ito0
end

c0 = sym('c0',[3,1,3]);
for n=1:3
    c0(:,:,n) = R0(:,:,n)*cm(:,:,n) + p0(:,:,n);
end

Jvc = sym('Jvc',[3,3,3]);
Jwc = sym('Jwc',[3,3,3]);
for i=1:3
    for j=1:3
        for k=1:3
            Jvc(i,j,k) = 0;
            Jwc(i,j,k) = 0;
        end
    end
end

n = 1;
while n<4
    for i=1:n
        Jvc(:,i,n) = cross(zc(:,:,i),(c0(:,:,n)-p0(:,:,i)));
    end
    n = n+1;
end

n = 1;
while n<4
    for i=1:n
        Jwc(:,i,n) = zc(:,:,i);
    end
    n = n+1;
end

Jvc = simplify(Jvc,"Steps",150);
Jwc = simplify(Jwc,"Steps",150);

D = sym('D',[3,3]);
for i=1:3
    for j=1:3
        D(i,j) = 0;
    end
end
for i=1:3
    D = D + (m(i) * (transpose(Jvc(:,:,i))*Jvc(:,:,i))) + (transpose(Jwc(:,:,i))*I(:,:,i)*Jwc(:,:,i));
end

Cf = sym('Cf',[3,3,3]);
for k=1:3
    for j=1:3
        for i=1:3
            Cf(i,j,k) = 0.5 * ( diff(D(i,k),q(j)) + diff(D(j,k),q(i)) - diff(D(i,j),q(k)) );
        end
    end
end

C = sym('C',[3,3]);
for i=1:3
    for j=1:3
        C(i,j) = 0;
    end
end
for k=1:3
    for i=1:3
        for j=1:3
            C(k,i) = C(k,i) + Cf(k,i,j)*qdot(j);
        end
    end
end

gravity = [0; 0; 9.780327];
syms P
P = 0;
for i=1:3
    P = P + transpose(gravity)*c0(:,:,i)*m(i);
end
g = sym('g',[3,1]);
for i=1:3
    g(i) = diff(P,q(i));
end

tau_L = D*qdd + C*qdot + g;
% tau_L = simplify(D*qdd + C*qdot + g, "steps", 50);