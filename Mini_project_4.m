

%
% mini-project 4 example script

clear all; close all;

zz=zeros(3,1); ex = [1;0;0]; ey = [0;1;0]; ez = [0;0;1];

% symbolic ABB IRB 1200 robot

syms L1 L2 L3 L4 L5 positive
syms q1 q2 q3 q4 q5 q6 real

%% define ABB IRB 1200 robot symbolically


%Finds the Positions
P01=0*ex+L1*ez;
P12=zz;
P23=L2*ez;
P34=L3*ez+L4*ex;
P45=zz;
P56=zz;
P6T=L5*ex;

%H
h1=ez;
h2=ey;
h3=ey;
h4=ex;
h5=ey;
h6=ex;

%q
q=[q1;q2;q3;q4;q5;q6];

%R
R01=rot(h1,q1);
R12=rot(h2,q2);
R23=rot(h3,q3);
R34=rot(h4,q4);
R45=rot(h5,q5);
R56=rot(h6,q6);


% forward kinematics and Jacobian

irb1200_s.P=[P01 P12 P23 P34 P45 P56 P6T];
irb1200_s.H=[h1 h2 h3 h4 h5 h6];
irb1200_s.joint_type=[0 0 0 0 0 0];
irb1200_s.q=q;

irb1200_s=fwddiffkiniter(irb1200_s);

% analytical Jacobian in base frame


%% 1B
R02=R01*R12;
R03=R02*R23;
R04=R03*R34;
R05=R04*R45;
R06=R05*R56;

R31 = (R12*R23)';
R32 = R23';
R33 = eye(3);
R43 = R34;
R35 = R34*R45;
R36 = R34*R45*R56;

P_14 = P12 + R12*P23 + R23*P34;
P_24 = P23 + R23*P34;
P_34 = P34;
P_44 = [0; 0 ;0];
P_54 = [0; 0 ;0];
P_64 = [0; 0 ;0];

P_4T = P34+ R34*P45 + R45*P56 + R56*P6T;


JT_0=simplify(irb1200_s.J);

 phi_T4 = phi(eye(3),P_4T);
% 
 JT_0=irb1200_s.J;
 R0T=simplify(irb1200_s.T(1:3,1:3));
 p0T=simplify(irb1200_s.T(1:3,4));
 J4_0=phi(eye(3,3),-R0T*P6T)*JT_0;
% 
 R02=rot(h1,q1)*rot(h2,q2);
% 
 J4_2=[R02' zeros(3,3);zeros(3,3) R02']*J4_0;
 J4_2(4:6,4:6)=simplify(J4_2(4:6,4:6));
 J4_2(1:3,4:6)=simplify(J4_2(1:3,4:6));
 J4_2(4:6,1:3)=simplify(J4_2(4:6,1:3));
 J4_2(1:3,1:3)=simplify(J4_2(1:3,1:3));
% pretty(J4_2);
% 
J4_3=simplify([R23' zeros(3,3);zeros(3,3) R23']*J4_2);
R03 = R01*R12*R23;
J4_3=simplify(phi(R03',-R0T*P6T)*JT_0);
% pretty(J4_3);
% 
% JT_0_n = simplify(phi(R03,R03'*R0T*P6T)*J4_3);
% pretty(JT_0_n);
% 
% c = simplify(JT_0-JT_0_n);


%% init

% define ABB IRB 1200 robot

%L
L1=399.1;
L2=448; %L2=350;
L3=42;
L4=451; %L4=351;
L5=82;

% P
p01=0*ex+L1*ez;p12=zz;p23=L2*ez;p34=L3*ez+L4*ex;p45=zz;p56=zz;p6T=L5*ex;

% H
h1=ez;h2=ey;h3=ey;h4=ex;h5=ey;h6=ex;

% 
irb1200.P=[p01 p12 p23 p34 p45 p56 p6T];
irb1200.H=[h1 h2 h3 h4 h5 h6];
irb1200.joint_type=[0 0 0 0 0 0];

% % define collision body for abb 1200
radius=.01;
[irb1200_rbt,colLink]=defineRobot(irb1200,radius);
 
% S-shaped curve

load S_sphere_path_uniform l pS
% 
% find the end effector frame
r=.5;
pc=r*ez;
N=length(pS);
xT=zeros(3,N);zT=zeros(3,N);yT=zeros(3,N);
quat=zeros(4,N); % allocate space for unit quaternion representation of R_{0T}

for i=1:N   
    xT(:,i)=(pS(:,i)-pc)/norm(pS(:,i)-pc);
    if i<N
        yT(:,i)=(pS(:,i+1)-pS(:,i));
    else
        yT(:,i)=yT(:,i-1);
    end
    yT(:,i)=yT(:,i)-yT(:,i)'*xT(:,i)*xT(:,i);
    yT(:,i)=yT(:,i)/norm(yT(:,i));
    zT(:,i)=cross(xT(:,i),yT(:,i));
    R=[xT(:,i) yT(:,i) zT(:,i)];
    quat(:,i)=R2q(R);
end


%% 1C
%comparing iterative Jacobian calculation in class vs. MATLAB
%

N=500;

for i=1:N
    
% random joint angles
q=(rand(6,1)-.5)*2*pi;

% iterative Jacobian calculation
irb1200.q=q;
tstart=tic;
irb1200=fwddiffkiniter(irb1200);
t(i)=toc(tstart);
J=irb1200.J;

% MATLAB function
tstart=tic;
J1=geometricJacobian(irb1200_rbt,q,'body7');
t1(i)=toc(tstart);
% check difference
diffJ1(i)=(norm(J-J1));

% using symbolic expression for Jacobian computation
R0T=irb1200.T(1:3,1:3);
tstart=tic;
J4_3_num=J4_3func(q);
R03=rot(irb1200.H(:,1),q(1))*rot(irb1200.H(:,2),q(2))*rot(irb1200.H(:,3),q(3));
R3T=rot(irb1200.H(:,4),q(4))*rot(irb1200.H(:,5),q(5))*rot(irb1200.H(:,6),q(6));
JT_0=phi(R03,R3T*irb1200.P(:,end))*J4_3_num;
t2(i)=toc(tstart);

diffJ2(i)=norm(J-JT_0);

end

figure(10);plot((1:N),t,'^',(1:N),t1,'o',(1:N),t2,'x','linewidth',2);
xlabel('run #');ylabel('sec');
title('Jacobian computation time');
legend('Iterative method','MATLAB','symbolic');
axis([1 N 0 .002]);

%% 2a

%Determinate of J4_3:
%L2*sin(q5)*(L4*cos(q3)+L3*sin(q3))*(L4*cos(q2+q3)+L3*sin(q2+q3)+L2*sin(q2));

Determinant = det(J4_3);
simplify(Determinant);

%wrist Singularity



%% PART 3

%forward  kin check

%
% proj4kincheck.m
%
% forward/inverse fkinematics checker for an elbow arm
% using MATLAB inversekinematics solver vs. subproblem decomposition
% 
% make sure fwdkin_example_rbt is in the path to generate the elbow_rbt
% rigid body tree
% 
% output plots compare EE solution accuracy and computation times
%

%clear all; close all;

% MATLAB robotics toolbox inverse kinematics object
ik = inverseKinematics('RigidBodyTree',irb1200_rbt);

% set storage
N=50;

n=length(irb1200.H);
q=zeros(n,N);
qsol=zeros(n,N);
qsol1=zeros(n,8,N);
qsol2=zeros(n,N);
T=zeros(4,4,N);
Tsol=zeros(4,4,N);
Tsol1=zeros(4,4,8,N);
Jsol2=zeros(6,6,N);
errT=zeros(N,1);
errT1=zeros(N,8);
errT2=zeros(N,1);
telapsed=zeros(1,N);
telapsed1=zeros(1,N);
telapsed2=zeros(1,N);

for i=1:N
    % random arm configuration
    q(:,i)=(rand(6,1)-.5)*pi;
    irb1200.q=q(:,i);
    % forward kinematics
    irb1200=fwddiffkiniter(irb1200);
    % nominal end effector pose
    T(:,:,i)=irb1200.T;
    % MATLAB's inverse kinematics
    % keep track of the computation time
    tstart=tic;
    [qsol(:,i),solnInfo]=...
        ik('body7',irb1200.T,ones(1,6),q(:,i)+q(:,i)*.1*randn);
    telapsed(i) = toc(tstart);    
    % forward kinematics again to compare with EE pose
    irb1200.q=qsol(:,i);
    irb1200=fwddiffkiniter(irb1200);
    Tsol(:,:,i)=irb1200.T;
    % find EE pose error
    errT(i)=norm(T(:,:,i)-Tsol(:,:,i),'fro');   
    % now use our own exact inverse kinematics
    irb1200.q=q(:,i);
    irb1200.T=T(:,:,i);
    % keep track of the computation time
    tstart=tic;
    irb1200=invkinelbow(irb1200);
    telapsed1(i) = toc(tstart);
    % there are 8 solutions
    qsol1(:,:,i)=irb1200.q;
    % find all the EE pose
    for j=1:8;
        irb1200.q=qsol1(:,j,i);
        irb1200=fwddiffkiniter(irb1200);
        Tsol1(:,:,j,i)=irb1200.T;
        Jsol1(:,:,j,i)=irb1200.J;
        errT1(i,j)=norm(T(:,:,i)-Tsol1(:,:,j,i),'fro');
    end
    % our own iterative inverse kinematics
    irb1200.MaxIter=100;
    irb1200.StepSize=.4;
    irb1200.Weights=[50;50;50;1;1;1];
    irb1200.q=q(:,i)+q(:,i).*.1.*randn(6,1);
    irb1200.T=T(:,:,i);
    tstart=tic;
    irb1200=invkin_iterJ(irb1200);
    telapsed2(i) = toc(tstart);
    qsol2(:,i)=irb1200.q;
    irb1200=fwddiffkiniter(irb1200);
    Tsol2(:,:,i)=irb1200.T;
    Jsol2(:,:,i)=irb1200.J;
    errT2(i)=norm(T(:,:,i)-Tsol2(:,:,i),'fro');
end 

figure(9);plot((1:N),errT,'bx',(1:N),errT2,'linewidth',2); hold on;
for j=1:8
    plot((1:N),errT1(:,j),'ro','linewidth',2);
end
hold off;
xlabel('random test number');ylabel('end effector error');
title('End effector pose error || T - T_{solve} ||_F');
legend('MATLAB','iterative Jacobian','exact by subproblem');

fprintf('max MATLAB EE error: %g \n',max(errT));
fprintf('max subproblem EE error: %g \n',max(max(errT1)));
fprintf('max iterative invkin EE error: %g \n',max(max(errT2)));

figure(20);plot((1:N),telapsed,'bx',(1:N),telapsed1,'ro',...,
    (1:N),telapsed2,'g^','linewidth',2);
xlabel('random test number');ylabel('elapsed time (sec)');
title('inverse kinematics computation time');
legend('MATLAB','exact by subproblem','iterative Jacobian');

fprintf('max and average MATLAB invkin computation time: %g, %g \n',...
    max(telapsed),mean(telapsed));
fprintf('max and average subproblem invkin computation time: %g, %g \n',...
    max(telapsed1),mean(telapsed1));
fprintf('max and average iterative invkin computation time: %g, %g \n',...
    max(telapsed2),mean(telapsed2));



%% Part 4

% plot the spherical S
load S_sphere_path
figure(4);plot3(p_S(1,:),p_S(2,:),p_S(3,:),'rx','linewidth',3);
%{
p_s(1,:) = x 
p_s(2,:) = y
p_s(3,:) = z
%}
clear T;

xlabel('x');ylabel('y');zlabel('z');
hold on;
% add a 3d sphere
surf(X,Y,Z)
% make it transparent
alpha .5
axis(r*[-1 1 -1 1 0 2]);axis('square');
view(120,10);

% convert to equal path length grid
diffS=vecnorm(diff(p_S')');
%path length
ls=[0 cumsum(diffS)];
%final path length
lf=sum(diffS);
N=100;
l=(0:lf/N:lf);

pS=interp1(ls,p_S',l,'spline')';
% plot it out again with equal path length
figure(2);plot3(pS(1,:),pS(2,:),pS(3,:),'rx','linewidth',3);
xlabel('x');ylabel('y');zlabel('z');
hold on;
% 3d sphere
surf(X,Y,Z)
% make it transparent
alpha 0.4
axis(r*[-1 1 -1 1 0 2]);axis('square');
view(120,10);

% check the path length is indeed equal
dlvec=vecnorm(diff(pS')');
figure(3);plot(dlvec,'x')
dl=mean(dlvec);
disp(max(abs(dlvec-dl)));

% save it as the path file
save S_sphere_path_uniform l pS

% define abb 1200 robot using POE convention
irb1200.P=[p01 p12 p23 p34 p45 p56 p6T]/1000;
irb1200.H=[h1 h2 h3 h4 h5 h6];
irb1200.joint_type=[0 0 0 0 0 0];
%irb1200.R6T=R6T;

% define collision body for abb 1200
radius=.01;
[irb1200_rbt,colLink]=defineRobot(irb1200,radius);
 
% 
% choose the inverse kinematics solution
%

for i=1:N
    % specify end effector SE(3) frame
    %Td{i}=[[xT(:,i) yT(:,i) zT(:,i)]*R6T' pS(:,i);[0 0 0 1]];
    Td{i}=[[xT(:,i) yT(:,i) zT(:,i)] pS(:,i);[0 0 0 1]];
    irb1200.T=Td{i};
    %
    irb1200=invkinelbow(irb1200); % << you need to supply this!!!
    %
    for k=1:8
        q(:,i,k)=irb1200.q(:,k);
    end
        % check forward kinematics to make sure the IK solution is correct
    for k=1:8
        irb1200.q=q(:,i,k);
        irb1200=fwddiffkiniter(irb1200);
        T{i,k}=irb1200.T;
    end

end

% choose the pose to visualize
ksol=1

for i=1:N
    % show robot pose (ever 5 frames)
    if mod(i,5)==0
        disp(norm(T{i,ksol}-Td{i}));
        figure(2);show(irb1200_rbt,q(:,i,ksol),'collision','on');
        view(150,10);
    end
end
% compute end effector linear and angular velcoity 

lsdot=.01;
for i=1:N-1
    dt(i) = (ls(i+1)-ls(i))/lsdot;
    for k=1:8
        qdot(:,i,k)=(q(:,i+1,k)-q(:,i,k))/dt(i);
        Ri1=T{i+1,k}(1:3,1:3);
        Ri=T{i,k}(1:3,1:3);
        w(:,i,k)=vee(Ri1*Ri'-eye(3,3))/dt(i);
        pi1=T{i+1,k}(1:3,4);
        pi=T{i,k}(1:3,4);
        v(:,i,k)=(pi1-pi)/dt(i);
    end
end

for k=1:8
   for i=1:6
       maxqdot(i,k)=max(qdot(i,:,k));
   end
   fprintf('maximum qdot for pose %d \n', k);
   disp(maxqdot(:,k)');   
end



%% functions
%save IRB1200Jacobian JT_0 J4_2 J4_3

%
% phi.m
% 
% propagation of spatial velocity
%
function phimat=phi(R,p)

    phimat=[R zeros(3,3);-R*hat(p) R];
    
end

%
% hat.m (converting a vector into a skew-symmetric cross-product matrix
%
% khat = hat(k)
%

function khat = hat(k)
  
  khat=[0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];

end

function robot = fwddiffkiniter(robot)
    n = length(robot.joint_type);
    joint_type = robot.joint_type;
    
    q = robot.q;
    H = robot.H;
    P = robot.P;
    T = eye(4,4);
    J = zeros(6,n);
    
    for i = 1:n
        h = H(1:3,i);
        if joint_type(i)==0
            R = expm(hat(h)*q(i));
            p = P(1:3,i);zeros(1,3);
        else
            R = eye(3,3);
            p = p(1:3,i)+q(i)*h;
        end
        J = phi(eye(3,3),T(1:3,1:3)*p)*J;
        
        if joint_type(i)==0
            J(:,i)=[T(1:3,1:3)*h;zeros(3,1)];
        else
            J(:,i)=[zeros(3,1);T(1:3,1:3)*h];
        end
        T = T*[R p;zeros(1,3) 1];
    end
    
    robot.J = phi(eye(3,3),T(1:3,1:3)*P(:,n+1))*J;
    robot.T = T*[eye(3,3) P(:,n+1);0 0 0 1];
end

function R = rot(k,theta)
    k=k/norm(k);
    R=eye(3,3)+sin(theta)*hat(k)+(1-cos(theta))*hat(k)*hat(k);
end 

function J=J4_3func(q)
L2=448;L3=42;L4=451;

J=[-sin(q(2) + q(3)),0, 0, 1, 0, cos(q(5));
    0, 1,   1, 0, cos(q(4)),  sin(q(4))*sin(q(5));
   cos(q(2) + q(3)), 0,   0, 0, sin(q(4)), -cos(q(4))*sin(q(5));
   0, L3 + L2*cos(q(3)),  L3, 0,       0, 0;
L4*cos(q(2) + q(3)) + L3*sin(q(2) + q(3)) + L2*sin(q(2)), 0, 0, 0,  0, 0;
   0, L2*sin(q(3)) - L4, -L4, 0, 0, 0];
     
end


function q=R2q(R)
  
  q=zeros(4,1);
  q(1)=.5*sqrt(trace(R)+1);
  if abs(q(1))<1e-5
    [k,theta]=R2kth(R);
    q(2:4)=k;
  else
    q(2:4)=vee(R-R')/4/q(1);
  end
end


function k = vee(K)

k=[-K(2,3);K(1,3);-K(1,2)];

end


% solve for q subtended between p1 and p2
%    k determines the sign of q
%
% input: k,p1,p2 as R^3 column vectors
% output: q (scalar)
%

function q=subprob0(k,p1,p2)

if ((k'*p1)>sqrt(eps)|(k'*p2)>sqrt(eps))
  error('k must be perpendicular to p and q');
end

p1=p1/norm(p1);
p2=p2/norm(p2);

q=2*atan2(norm(p1-p2),norm(p1+p2));

if k'*(cross(p1,p2))<0
  q=-q;
end 
end

%
% q=subprob1(k,p1,p2)
%
% solve for q from
%
% exp(k x q) p1 = p2
%
% input: k,p1,p2 as R^3 column vectors
% output: q (scalar)
%

function q=subprob1(k,p1,p2)

p2=p2/norm(p2)*norm(p1);

if norm(p1-p2)<sqrt(eps);q=0;return;end
  
k=k/norm(k);
pp1=p1-(p1'*k)*k;
pp2=p2-(p2'*k)*k;

epp1=pp1/norm(pp1);
epp2=pp2/norm(pp2);

q=subprob0(k,epp1,epp2);
%q=atan2(k'*(cross(epp1,epp2)),epp1'*epp2);
end

%[q1,q2]=subprob2(k1,k2,p1,p2)
%
% solve for theta1 and theta2 from
%
% exp(k1 x q1) p1 = exp(k2 x q2) p2 
%  
% input: k1,k2,p1,p2 as R^3 column vectors
%
% output: q1 and q2 as 2x1 columns corresponding to the two solutions
%

function [q1,q2]=subprob2(k1,k2,p1,p2)

p2=p2/norm(p2)*norm(p1);
k12=k1'*k2;
pk1=p1'*k1;
pk2=p2'*k2;

% check if solution exists

if abs(k12^2-1)<eps;theta1=[];theta2=[];
    q1=[NaN;NaN];q2=[NaN;NaN];
    disp('no solution (k1 and k2 are collinear)');
    return;
end

a=[1 -k12; -k12 1]*[pk1;pk2]/(1-k12^2);

% 
% check if solution exists
%
cond=(norm(p1)^2-norm(a)^2-2*a(1)*a(2)*k12);

% special case: 1 solution
if abs(cond)<eps;
  v=[k1 k2]*a;
  q1a=subprob1(k1,p1,v);
  q2a=subprob1(k2,p2,v);
  q1=[q1a;q1a];
  q2=[q2a;q2a];
end

% special case: no solution
if cond<0
    q1=[NaN NaN];q2=[NaN NaN];
    disp('no solution (two cones do not intersect)');
    return;
end

gamma=sqrt(cond)/norm(cross(k1,k2));

% general case: 2 solutions

q1=zeros(2,1);
q2=zeros(2,1);

v1=[k1 k2 cross(k1,k2)]*[a;gamma];
v2=[k1 k2 cross(k1,k2)]*[a;-gamma];
q1(1)=subprob1(k1,p1,v1);
q1(2)=subprob1(k1,p1,v2);

q2(1)=subprob1(k2,p2,v1);
q2(2)=subprob1(k2,p2,v2);

end

%
% q=subprob3(k,p1,p2,d)
%
% solve for theta from
%
% norm(p2-exp(k x q) p1) = d
%
% input: k,p1,p2 as R^3 column vectors, delta: scalar
% output: q (2x1 vector, 2 solutions)
%

function q=subprob3(k,p1,p2,d)

pp1=p1-k'*p1*k;
pp2=p2-k'*p2*k;
dpsq=d^2-(k'*(p1-p2))^2;

if dpsq<0;theta=[NaN;NaN];return;end

if dpsq==0;theta=subprob1(k,pp1/norm(pp1),pp2/norm(pp2));return;end
  
bb=(norm(pp1)^2+norm(pp2)^2-dpsq)/(2*norm(pp1)*norm(pp2));
if abs(bb)>1; theta=[NaN;NaN];return;end

phi=acos(bb);

q0=subprob1(k,pp1/norm(pp1),pp2/norm(pp2));
q=zeros(2,1);

q(1)=q0+phi;
q(2)=q0-phi;

end

function robot = invkinelbow(robot)
%     ex = [1;0;0]; ey = [0;1;0]; ez = [0;0;1];
%     T=robot.T;
%     R0T = T(1:3,1:3);
%     p0T = T(1:3,4);
%     
%     h1=ez;
%     h2=ey;
%     h3=ey;
%     h4=ex;
%     h5=ey;
%     h6=ex;
% 
%     
%     p01 = robot.P(:,1);
%     p12 = robot.P(:,2);
%     p23 = robot.P(:,3);
%     p34 = robot.P(:,4);
%     p45 = robot.P(:,5);
%     p56 = robot.P(:,6);
%     p6T = robot.P(:,7);
%     
%     
%    p16_0 = p0T-p01-R0T*p6T;
%    
%    %solves for q3_1 and q3_2
%    q3 = subprob3(h3,-p34,p23,norm(p16_0));
%    %q3_1 = q3(1);
%    %q3_2 = q3(2);
%    
%    %         sol 1      sol 2       sol 1      sol 2
%    q3 = [q3(1) q3(1) q3(2) q3(2) q3(1) q3(1) q3(2) q3(2)];
%    
%    %solves Rotataion y matrix for q3_1
%    %Ry_q3_1 = rot(ey,q3_1);
%    %solves Rotation y matrix for q3_2
%    %Ry_q3_2 = rot(ey,q3_2);
%    %p2 for subproblem
%    %p2_1 = (p23+Ry_q3_1*p34);
%    %p2_2 = (p23+Ry_q3_2*p34);
%    
%    %solves for solution 1 of q1 and q2
%    [q1_1,q2_1] = subprob2(-h1,h2,p16_0,(p23+(rot(h3,q3(1))*p34)));
%    %solves for solution 2 of q1 and q2
%    [q1_2,q2_2] = subprob2(-h1,h2,p16_0,(p23+(rot(h3,q3(3))*p34)));
%    
%    %      sol 1   sol 2  sol 3   sol 4        
%    q1 = [q1_1(1) q1_1(2) q1_2(1) q1_2(2) q1_1(1) q1_1(2) q1_2(1)  q1_2(2)];
%    q2 = [q2_1(1) q2_1(2) q2_2(1) q2_2(2) q2_1(1) q2_1(2) q2_2(1)  q2_2(2)];
%    
%    R03_1 = rot(h1,q1(1))* rot(h2,q2(1))*rot(h3, q3(1));
%    R03_2 = rot(h1,q1(2))* rot(h2,q2(2))*rot(h3, q3(2));
%    R03_3 = rot(h1,q1(3))* rot(h2,q2(3))*rot(h3, q3(3));
%    R03_4 = rot(h1,q1(4))* rot(h2,q2(4))*rot(h3, q3(4));
%    
%    
%    %q4 and q5 using solution 1 
%    [q4_1,q5_1] = subprob2(-h4,h5,(R03_1')* R0T * h6,h6);
%    [q4_2,q5_2] = subprob2(-h4,h5,(R03_2')* R0T * h6,h6);
%    [q4_3,q5_3] = subprob2(-h4,h5,(R03_3')* R0T * h6,h6);
%    [q4_4,q5_4] = subprob2(-h4,h5,(R03_4')* R0T * h6,h6);
%    
%    q4 = [q4_1(1) q4_2(1) q4_3(1) q4_4(1) q4_1(2) q4_2(2) q4_3(2) q4_4(2)];
%    q5 = [q5_1(1) q5_2(1) q5_3(1) q5_4(1) q5_1(2) q5_2(2) q5_3(2) q5_4(2)];
%    
%    
%    P2_1 = R03_1*rot(h4,q4(1))*rot(h5,(q5(1)));
%    P2_2 = R03_2*rot(h4,q4(2))*rot(h5,(q5(2))); 
%    P2_3 = R03_3*rot(h4,q4(3))*rot(h5,(q5(3)));
%    P2_4 = R03_4*rot(h4,q4(4))*rot(h5,(q5(4)));
%    P2_5 = R03_1*rot(h4,q4(5))*rot(h5,(q5(5)));
%    P2_6 = R03_2*rot(h4,q4(6))*rot(h5,(q5(6)));
%    P2_7 = R03_3*rot(h4,q4(7))*rot(h5,(q5(7)));
%    P2_8 = R03_4*rot(h4,q4(8))*rot(h5,(q5(8)));
%    
%    q6_1 = subprob1(h6,h5,P2_1'*R0T*h5);
%    q6_2 = subprob1(h6,h5,P2_2'*R0T*h5);
%    q6_3 = subprob1(h6,h5,P2_3'*R0T*h5);
%    q6_4 = subprob1(h6,h5,P2_4'*R0T*h5);
%    q6_5 = subprob1(h6,h5,P2_5'*R0T*h5);
%    q6_6 = subprob1(h6,h5,P2_6'*R0T*h5);
%    q6_7 = subprob1(h6,h5,P2_7'*R0T*h5);
%    q6_8 = subprob1(h6,h5,P2_8'*R0T*h5);
%    
%    q6 = [q6_2 q6_1 q6_4 q6_3 q6_6 q6_5 q6_8 q6_7];
%    
%    q_sol = [q1; q2; q3; q4; q5; q6];
%     
%    robot.q = q_sol
%     

%create shorthand for all the robot variables that will be used for
%calculations
p01 = robot.P(:,1);
p12 = robot.P(:,2);
p23 = robot.P(:,3);
p34 = robot.P(:,4);
p45 = robot.P(:,5);
p56 = robot.P(:,6);
p6T = robot.P(:,7);

h1=robot.H(:,1);
h2=robot.H(:,2);
h3=robot.H(:,3);
h4=robot.H(:,4);
h5=robot.H(:,5);
h6=robot.H(:,6);

R = robot.T(1:3, 1:3); p=robot.T(1:3,4);

%solve for q3 using subproblem 3
q3sol = subprob3(h3, -p34, p23, norm(p-p01-R*p6T));

%solve for the two distinct solutions of q1 q2 using subproblem 2
    for i=1:2
        [q1sol(:,i),q2sol(:,i)]=subprob2(-h1,h2,p-p01-R*p6T,p23+rot(h3,q3sol(i))*p34);
    end

    qsol = zeros(6,8);

    qsol(1:3,1:4) = [q1sol(:,1)' q1sol(:,2)'; ...
    q2sol(:,1)' q2sol(:,2)' ; ...
    q3sol(1) q3sol(1) q3sol(2) q3sol(2)];
    qsol(1:3,5:8) = qsol(1:3,1:4);

%solve for the four possible wrist angles
    for i=1:4
        %solve for q4 and q5 using subproblem 2
        R03 = rot(h1,qsol(1,i)) * rot(h2,qsol(2,i)) * rot(h3,qsol(3,i));
        [q4vec, q5vec] = subprob2(-h4, h5,R03'*R*h6, h6);
        q4a = q4vec(1); q4b = q4vec(2);
        q5a = q5vec(1); q5b = q5vec(2);
        qsol(4:5,i)= [q4a;q5a];
        qsol(4:5,i+4)= [q4b;q5b];

        R05a=R03*rot(h4, q4a)*rot(h5,q5a);
        R05b=R03*rot(h4, q4b)*rot(h5,q5b);

        %solve for the q6 values for each respective q4/q5 pair using
        %subproblem 1
        qsol(6,i) = subprob1(h6, h5, R05a'*R*h5);
        qsol(6,i+4) = subprob1(h6, h5, R05b'*R*h5);
    end

    %overwrite the solutions for q into the given robot
    robot.q=qsol;

end


% invkin_iterJ.m
%
% inverse kinematics using Jacobian iteration (planar)
%

function robot=invkin_iterJ(robot,N,alpha)

    % parameters
    N=robot.MaxIter;
    alpha=robot.StepSize;
    weights =robot.Weights;
    % target (R,p)
    p0Td=robot.T(1:3,4);    
    R0Td = robot.T(1:3,1:3);
        
    % set up storage space
    q0=robot.q;
    n=length(q0); % # of joints
    q=zeros(n,N+1); 
    q(:,1)=q0; % output joint displacements
    p0T=zeros(3,N+1); % output p
    quatT=zeros(4,N+1); % output quaternion

    % iterative update
    for i=1:N
        % forward kinematics
        robot.q=q(:,i);
        robot=fwddiffkiniter(robot);
        quatT(:,i)=R2q(robot.T(1:3,1:3));
        p0T(:,i)=robot.T(1:3,4);  
        % task space error: angular error and position error
        dX=[R2qv(robot.T(1:3,1:3)*R0Td');p0T(:,i)-p0Td];      
        % Jacobian update - note 10 times the gain for orientation        
        % qq=q(:,i)-alpha*pinv(robot.J)*dX;
        qq=q(:,i)-alpha*robot.J'*inv(robot.J*robot.J'+.01*diag(1./weights))*dX;
        q(:,i+1)=(qq>pi).*(-2*pi+qq)+(qq<-pi).*(2*pi+qq)+(qq<pi).*(qq>-pi).*qq;
    end
    % final iteration
    robot.q=q(:,N+1);
    robot=fwddiffkiniter(robot);
end

%
% R2qv.m
%
% converts R in SO(3) to vector quaternion q_vec
%

function qv=R2qv(R)

q01=.5*sqrt(trace(R)+1);
q02=-.5*sqrt(trace(R)+1);
    if abs(q01)<1e-5
    [k,theta]=R2kth(R);
    qv=k;
    else
    qv(:,1)=vee(R-R')/4/q01;
    qv(:,2)=vee(R-R')/4/q02;
    end
 qv = qv(:,1);
end
