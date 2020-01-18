clc;clear;close all;
% load track information
load('TestTrack2.mat');
load('LeadCarA_input_init.mat')
load('LeadCarB_input_init.mat')
load('s_distance.mat')
% read the track data 
range = 1:246;
centerx = TestTrack.cline(1,range); centery = TestTrack.cline(2,range);
blx = TestTrack.bl(1,range);        bly = TestTrack.bl(2,range);
brx = TestTrack.br(1,range);        bry = TestTrack.br(2,range);
rlx = TestTrack.leftlane(1,range);  rly = TestTrack.leftlane(2,range);
rrx = TestTrack.rightlane(1,range); rry = TestTrack.rightlane(2,range);
head = TestTrack.theta(1,range);
% poly-fit the track data 
n_poly = range(end)*5;
x_poly = spline(1:range(end),centerx,linspace(1,range(end),n_poly));        y_poly = spline(1:range(end),centery,linspace(1,range(end),n_poly));
blx_poly = spline(1:range(end),blx,linspace(1,range(end),n_poly));          bly_poly = spline(1:range(end),bly,linspace(1,range(end),n_poly));
brx_poly = spline(1:range(end),brx,linspace(1,range(end),n_poly));          bry_poly = spline(1:range(end),bry,linspace(1,range(end),n_poly));
rlx_poly = spline(1:range(end),rlx,linspace(1,range(end),n_poly));          rly_poly = spline(1:range(end),rly,linspace(1,range(end),n_poly));        
rrx_poly = spline(1:range(end),rrx,linspace(1,range(end),n_poly));          rry_poly = spline(1:range(end),rry,linspace(1,range(end),n_poly));        
head_poly = spline(1:range(end),head,linspace(1,range(end),n_poly));
% construct reference inputs and states 
U_ref = zeros(2,n_poly-1);
U_ref(2,:) = 1000;
Y_ref = zeros(4,n_poly);
Y_ref(1,:) = 10;
Y_ref(2,:) = rrx_poly(1:n_poly);
Y_ref(3,:) = rry_poly(1:n_poly);
Y_ref(4,:) = head_poly(1:n_poly);

% plot the track
figure(1)
subplot(1,2,1)
plot(blx_poly,bly_poly,'k',brx_poly,bry_poly,'k')
axis equal;grid on;hold on
plot(x_poly,y_poly,'k--');
axis([200,1500,-300,1000])
subplot(1,2,2)
plot(blx_poly,bly_poly,'k',brx_poly,bry_poly,'k')
axis equal;grid on;hold on
plot(x_poly,y_poly,'k--');


%% Switch
flag_change = 0;        % vehicle has not changed traj
flag_back = 0;          % vehicle has not moved back 
flag_ready_change = 1;  % vehicle is ready to change
flag_ready_back = 0;    % vehicle is not ready to move back

%% MPC
% MPC parameter
dt = 0.01;
T = 0:dt:25000*dt;   % estimate to finish the track within 250 seconds, actually 208.27 seconds 
npred = 10;
ndec = 4*(npred+1)+2*npred;
A = @(i) eye(4)+dt*...
    [                                                               -1/20, 0, 0,                                                                                   0;
    cos(U_ref(1,i))*cos(Y_ref(4,i)) - (sin(U_ref(1,i))*sin(Y_ref(4,i)))/2, 0, 0, -Y_ref(1,i)*(cos(U_ref(1,i))*sin(Y_ref(4,i)) + (cos(Y_ref(4,i))*sin(U_ref(1,i)))/2);
    cos(U_ref(1,i))*sin(Y_ref(4,i)) + (cos(Y_ref(4,i))*sin(U_ref(1,i)))/2, 0, 0, -Y_ref(1,i)*((sin(U_ref(1,i))*sin(Y_ref(4,i)))/2 - cos(U_ref(1,i))*cos(Y_ref(4,i)));
                                                    (5*sin(U_ref(1,i)))/8, 0, 0,                                                                                   0];

B = @(i) dt*...
    [                                                                                 0, 1/2000;
    -Y_ref(1,i)*((cos(U_ref(1,i))*sin(Y_ref(4,i)))/2 + cos(Y_ref(4,i))*sin(U_ref(1,i))),      0;
    -Y_ref(1,i)*(sin(U_ref(1,i))*sin(Y_ref(4,i)) - (cos(U_ref(1,i))*cos(Y_ref(4,i)))/2),      0;
                                                       (5*Y_ref(1,i)*cos(U_ref(1,i)))/8,      0];
% cost function
Q = diag([1500 10 10 1]);R = [1 0;0 1];
H = blkdiag(kron(eye(npred+1),Q),kron(eye(npred),R));
% actual trajectory 
Y_act = zeros(4,length(T));
Y_act(:,1) = [10;290.81;-177.42;1.96];
U_act = zeros(2,length(T));
% initialize other vehicle
Y_act_A = zeros(4,length(T));
Y_act_A(:,1) = LeadCarA.init;
Y_act_B = zeros(4,length(T));
Y_act_B(:,1) = LeadCarB.init;
CarX0 = plot(290.81,-177.42);       CarX1 = CarX0;
CarA0 = plot(LeadCarA.init(2:3));   CarA1 = CarA0;
CarB0 = plot(LeadCarB.init(2:3));   CarB1 = CarB0;

for i = 1:length(T)
    % Here we look to see what the closest y_ref is to where we currently
    % are. In general, MPC tracks behind the reference trajectory so we 
    % reset the reference for the error calc to be close to where i is.
    k = knnsearch([Y_ref(2,:);Y_ref(3,:)]',[Y_act(2,i),Y_act(3,i)]);
    if k == n_poly-10
        break
    end
    
    % implementing reference trajectory altering algorithm
    if flag_change == 1
        flag_change = 0;
        flag_ready_change = 0;   
        flag_ready_back = 1;
        disp('change traj')
        [rrx_new,rry_new] = change_traj_right_to_left(k,rlx_poly,rly_poly,rrx_poly,rry_poly);
        
        U_ref(2,:) = 1200;
        Y_ref(1,k+1:k+9) = 12;
        Y_ref(2,k+1:k+9) = rrx_new;
        Y_ref(3,k+1:k+9) = rry_new;
        
        Y_ref(1,k+10:end) = 12;
        Y_ref(2,k+10:end) = rlx_poly(k+10:end);
        Y_ref(3,k+10:end) = rly_poly(k+10:end);
        
        passing_moment = i;
    elseif flag_back == 1
        flag_back = 0;
        flag_ready_back = 0;
        disp('go back')
        [rlx_new,rly_new] = change_traj_left_to_right(k,rlx_poly,rly_poly,rrx_poly,rry_poly);
        
        U_ref(2,:) = 1000;
        Y_ref(1,k+1:k+9) = 10;
        Y_ref(2,k+1:k+9) = rlx_new;
        Y_ref(3,k+1:k+9) = rly_new;
        
        Y_ref(1,k+10:end) = 10;
        Y_ref(2,k+10:end) = rrx_poly(k+10:end);
        Y_ref(3,k+10:end) = rry_poly(k+10:end);
    end
    if ~exist('passing_moment','var')
    elseif i == passing_moment+2000
        flag_ready_change = 1;
    end
    % changing input based on new reference trajectory
    Y_error = Y_act(:,i)-Y_ref(:,k);
    [low,up] = bound_cons(k,U_ref,npred,ndec,[-0.5,0.5,-5000,2500]);
    [Aeq,beq] = eq_cons(k,A,B,npred,ndec,Y_error);
    options = optimoptions('quadprog','ConstraintTolerance',10^-4,'OptimalityTolerance',10^-3,'Display','off');
    z_hor = quadprog(H,[],[],[],Aeq,beq,low,up,[],options);
    U_act(:,i) = z_hor(45:46)+U_ref(:,k);
    % simulating main vehicles
    [~,Y_temp] = ode45(@(t,x) vehicle(x,U_act(:,i)),0:dt:dt,Y_act(:,i));
    Y_act(:,i+1) = Y_temp(end,:)';
    % simulating other vehicles
    [~,Y_temp_A] = ode45(@(t,x) vehicle(x,LeadCarA.input(:,i)),0:dt:dt,Y_act_A(:,i));
    Y_act_A(:,i+1) = Y_temp_A(end,:)';
    [~,Y_temp_B] = ode45(@(t,x) vehicle(x,LeadCarB.input(:,i)),0:dt:dt,Y_act_B(:,i));
    Y_act_B(:,i+1) = Y_temp_B(end,:)';
    % calculating distance along the path between vehicles
    S = s(k);
    S_A = s(finding_index(Y_act_A(2,i),Y_act_A(3,i),x_poly,y_poly));
    S_B = s(finding_index(Y_act_B(2,i),Y_act_B(3,i),x_poly,y_poly));
    % implementing switch system
    if 0 < S_A-S && S_A-S < 50 && flag_ready_change == 1
        flag_change = 1;
        disp('Approaching Car A')
    end
    if -10 < S_A-S && S_A-S < -5 && flag_ready_back == 1
        flag_back = 1;
        disp('Car A is passed')
    end
    if 0 < S_B-S && S_B-S < 50 && flag_ready_change == 1
        flag_change = 1;
        disp('Approaching Car B')
    end
    if -10 < S_B-S && S_B-S < -5 && flag_ready_back == 1
        flag_back = 1;
        disp('Car B is passed')
    end
    % plotting every 25 steps 
    if (mod(i,25) == 0)
        subplot(1,2,1)
        [CarX0,CarA0,CarB0] = plot_vehicle(Y_act,Y_act_A,Y_act_B,i,CarX0,CarA0,CarB0,0);
        subplot(1,2,2)
        [CarX1,CarA1,CarB1] = plot_vehicle(Y_act,Y_act_A,Y_act_B,i,CarX1,CarA1,CarB1,1);
        string2 = ['Speed is ',num2str(Y_act(1,i))];
        title(string2)
    end
end

%% Sub-function
function [rrx_new,rry_new] = change_traj_left_to_right(index,rlx_poly,rly_poly,rrx_poly,rry_poly)
% passing seven points
percentage = [0.95 0.9 0.8 0.65 0.5 0.32 0.2 0.1 0.05];
rrx_new = zeros(1,length(percentage));
rry_new = zeros(1,length(percentage));
for i = 1:length(percentage)
    percent = percentage(i);
    [x_new,y_new] = point_shift(percent,index+i,rlx_poly,rly_poly,rrx_poly,rry_poly);
    rrx_new(i) = x_new;
    rry_new(i) = y_new;
end
end

function [rrx_new,rry_new] = change_traj_right_to_left(index,rlx_poly,rly_poly,rrx_poly,rry_poly)
% passing seven points
percentage = [0.05 0.1 0.2 0.35 0.5 0.65 0.8 0.9 0.95];
rrx_new = zeros(1,length(percentage));
rry_new = zeros(1,length(percentage));
for i = 1:length(percentage)
    percent = percentage(i);
    [x_new,y_new] = point_shift(percent,index+i,rlx_poly,rly_poly,rrx_poly,rry_poly);
    rrx_new(i) = x_new;
    rry_new(i) = y_new;
end
end

function [x_new,y_new] = point_shift(percent,index,rlx_poly,rly_poly,rrx_poly,rry_poly)
lx = rlx_poly(index);
ly = rly_poly(index);
rx = rrx_poly(index);
ry = rry_poly(index);

x_new = (lx-rx)*percent+rx;
y_new = (ly-ry)*percent+ry;
end

function index = finding_index(x,y,x_poly,y_poly)
position = [x_poly;y_poly]';
Req_position = [x,y];
index = knnsearch(position,Req_position);
end

function [CarX,CarA,CarB] = plot_vehicle(Y_act,Y_act_A,Y_act_B,i,CarX,CarA,CarB,condition)
delete(CarX);delete(CarA);delete(CarB)
CarX = plot(Y_act(2,i),Y_act(3,i),'bh','linewidth',2);
CarA = plot(Y_act_A(2,i),Y_act_A(3,i),'ro','linewidth',2);
CarB = plot(Y_act_B(2,i),Y_act_B(3,i),'ro','linewidth',2);
if condition == 1
    axis([Y_act(2,i)-30, Y_act(2,i)+30, Y_act(3,i)-30, Y_act(3,i)+30]);
end
drawnow
end

function [Aeq,beq] = eq_cons(idx,A,B,npred,Ndec,Zini)
Aeq = -eye((npred+1)*4,Ndec);
beq = zeros((npred+1)*4,1);
Aeq(1:4,1:4) = eye(4);
beq(1:4) = Zini;
for i = 0:npred-1
    Aeq((i+1)*4+(1:4),i*4+(1:4)) = A(idx+i);
    Aeq((i+1)*4+(1:4),(npred+1)*4+i*2+(1:2)) = B(idx+i);
end
end

function [Lb,Ub] = bound_cons(idx,U_ref,npred,Ndec,input_range)
Lb = -inf(Ndec,1);
Ub = inf(Ndec,1);
for i = 0:npred-1
    Lb(44+i*2+1) = input_range(1)-U_ref(1,idx+i);
    Ub(44+i*2+1) = input_range(2)-U_ref(1,idx+i);
    Lb(44+i*2+2) = input_range(3)-U_ref(2,idx+i);
    Ub(44+i*2+2) = input_range(4)-U_ref(2,idx+i);
end
end

function state_dot = vehicle(state,input)
% Constants
L1 = 1.6;   L2 = 0.8;   R = 1/5;
Ak = 100;   m = 2000;   b = 100;
% Inputs
% delta = interp1(T,input(:,1),t,'previous','extrap');
% Fd = interp1(T,input(:,2),t,'previous','extrap');
delta = input(1);Fd = input(2);
% States
u = state(1);x = state(2);y = state(3);psi = state(4);
% Outputs
u_dot = (Fd-b*u)/m;
x_dot = (((-L2/L1)*sin(delta)*sin(psi)) + (cos(delta)*cos(psi)))*u;
y_dot = (((L2/L1)*sin(delta)*cos(psi)) + (cos(delta)*sin(psi)))*u;
psi_dot = (1/L1)*sin(delta)*u;
state_dot = [u_dot;x_dot;y_dot;psi_dot];
end


