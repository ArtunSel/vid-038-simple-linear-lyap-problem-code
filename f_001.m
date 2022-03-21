
%% find equilibrium-point using "fimplicit"
clear all,close all,clc;
syms x1 x2 x3 real
f1=-x1-10*x2;
f2=-x2;

range1=.1;
fimplicit(@(x1,x2)-x1-10*x2,range1*[-1,1,-1,1]); hold on;
fimplicit(@(x1,x2)-x2,range1*[-1,1,-1,1]);
%% find equilibrium-point using "vpasolve"
x1e = []; x2e = [];
try
%     [x1e,x2e,x3e] = vpasolve([f1;f2;f3],[x1;x2;x3]);
    [x1e,x2e] = vpasolve([f1;f2],[x1;x2]);
    ind = imag(x1e)==0 & imag(x2e)==0; % Take only real solutions
    x1e = double(x1e(ind));
    x2e = double(x2e(ind));
catch
    x1e = []; x2e = [];
end
x1e
x2e

%% part-1 use "quiver" to visualize the vec-fields
clear all,close all,clc;
range1=10;
[x1,x2]=meshgrid(range1*[-1:0.1:1],range1*[-1:0.1:1]);
x1dot=-x1-10*x2;
x2dot=-x1-10*x2;
quiver(x1,x2,x1dot,x2dot) 
%% check the eigenvalues
A=[-1,-10;0,-1];
eig([-1,-10;0,-1])
P=lyap(A',1*eye(2))
eig(P)
eig([P*A]+[P*A]')
%% let us try the generic-lyap fcn
clear all,close all,clc;
syms x1 x2 real
f1=-x1-10*x2;
f2=-x2;
V=0.5*(x1^2+x2^2);

% Vdot=gradient(V)'*[f1;f2];
Vdot=jacobian(V,[x1;x2])*[f1;f2];
Vdot=simplify(Vdot)
Vdot=expand(Vdot)
R=1;
theta_vec=linspace(-pi,pi,100);
max_pt=-1e5;
for ii=1:1:length(theta_vec)
    temp1=subs(Vdot,[x1,x2],...
        [R*cos(theta_vec(ii)),R*sin(theta_vec(ii))]);
    temp1=double(temp1);
    if (temp1>max_pt) max_pt=temp1; end 
%     if ~(temp1<=0)
%         error('positive val detected->BAD');
%     end
end
max_pt


%% let us try the "x'Px" as the lyap-fcn
clear all,close all,clc;
syms x1 x2 real
f1=-x1-10*x2;
f2=-x2;
A=[-1,-10;0,-1];
P=lyap(A',1*eye(2))
V=[x1;x2]'*P*[x1;x2];
% Vdot=gradient(V)'*[f1;f2];
Vdot=jacobian(V,[x1;x2])*[f1;f2];
Vdot=simplify(Vdot)
R=1;
theta_vec=linspace(-pi,pi,100);
max_pt=-100;
for ii=1:1:length(theta_vec)
    temp1=subs(Vdot,[x1,x2],...
        [R*cos(theta_vec(ii)),R*sin(theta_vec(ii))]);
    temp1=double(temp1);
    if (temp1>max_pt) max_pt=temp1; end 
%     if ~(temp1<=0)
%         error('positive val detected->BAD');
%     end
end
max_pt

%% let us look at a single trajectory
clear all,close all,clc;
fig1=figure(1);fig1.Color=[1,1,1];
ax1=axes('Parent',fig1);
    set(0,'CurrentFigure',fig1);
    set(fig1,'currentaxes',ax1);
tspan=[0:0.01:30]; x0=randi([-10,10],2,1);
wt=tspan;
f=randi([1,10],1,1); w=sin(2*pi*f*tspan);
% w=square(tspan);
[t,x]=ode45(@(t,x) odefcn(t,x,wt,w),tspan,x0);

% plot(t,x(:,1),'r-','LineWidth',[2],"Parent",ax1);hold on;
% plot(t,x(:,2),'r-','LineWidth',[2],"Parent",ax1);hold on;
% yline(1);hold on;yline(-1);hold on;

A=[-1,-10;0,-1];
P=lyap(A',1*eye(2))
% V=[x1;x2]'*P*[x1;x2];
V=zeros(size(x(:,1)));
for ii=1:1:length(V)
V(ii)=[x(ii,1);x(ii,2)]'*P*[x(ii,1);x(ii,2)];
end
plot(t,V,'r-','LineWidth',[2],"Parent",ax1);

function xdot=odefcn(t,x,wt,w)
w=interp1(wt,w,t);
xdot=zeros(2,1);
x1=x(1);x2=x(2);
xdot(1)=-x1-10*x2;
xdot(2)=-x2;
end























%