close all
Ix = 0.0028;                % rotational moment of inertia
Iy = 0.0501;                % rotational moment of inertia
Iz = 0.0501;                % rotational moment of inertia
m = 10;                     % mass of the rocket in kg
rho = 1.17;                 % density in kg/m^3
S = 0.0799846*0.08155432;   % wing area in m^2
d = 0.1;                    % distance from roll axis to aerodynamic center
Cldf = 0.01;                % coefficient of lift relative to fin angle
p0 = 0.01;                  % initial roll rate rad/s
q0 = 0.01;                  % initial pitch rate rad/s
r0 = 0.01;                  % initial yaw rate rad/s
v = 70;                     % tunnel velocity m/s


syms p pdot q qdot r rdot df

xdot = [-q*r*(Iz-Iy)/Ix     %pdot
        -p*r*(Ix-Iz)/Iy     %qdot
        -q*p*(Iy-Ix)/Iz];   %rdot

Asyms = jacobian(xdot,[p,q,r])  % symbolic A

A = double(subs(subs(subs(Asyms,p,p0),q,q0),r,r0))  % numerical A

B = [3*0.5*rho*S*d*v^2*Cldf/Ix
        0
        0               ];

C = eye(3);

Q = (diag([1/0.5,1/0.5,1/0.5]).^2);

R = (diag(1/9)).^2;

[K,P,E] = lqr(A,B,Q,R);
N = 0;  % no H matrix for controlled variable model
        % z(t) = Gx +Hu     N = G'H

Acontrolled = (A-B*K)

% Uncontrolled Eigenvalues
eig(A)

% Controlled Eigenvalues
eigen = eig(Acontrolled)

sys = ss((Acontrolled),B,C,0,'OutputName',{'P' 'Q' 'R'},'StateName',{'P' 'Q' 'R'},...
    'InputName','Fin Angle');

% Impulse response to show stability
l = 0;
ts = [0,0,0;];
for a = [-9 0 9]
l = l+1;
figure
t = linspace(0,300,100000);
[y] = impulse(a*sys,t);
subplot(3,1,1)
plot(t,y(:,1))
subplot(3,1,2)
plot(t,y(:,2))
subplot(3,1,3)
plot(t,y(:,3))

for i = 1:length(y)
    if abs(y(length(y)-i+1,1))<=0.01   % Settle threshold
        ts(l) = t(length(y)-i+1);  % Settle time for each angle impulse
    end
end
end

% Controllabliity Test, must be full rank (3)
rank([B , (Acontrolled)*B, (Acontrolled)^2*B])

% Observability Test, must be full rank (3)
rank([C , (Acontrolled)*C, (Acontrolled)^2*C])

% Hamiltonian Matrix
H = [A-B*R^(-1)*N' , -B*R^(-1)*B'
     -Q+N*R^(-1)*N' , -(A-B*R^(-1)*N')']





