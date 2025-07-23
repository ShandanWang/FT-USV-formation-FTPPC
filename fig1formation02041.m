clc
clear
close all
%% initialize of Graph
Nf = 6; % number of followers
Nl = 1; % number of leaders
W = [0, 1, 0, 0, 1, 0;
     1, 0, 0, 0, 1, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 1, 0, 1, 0;
     1, 1, 0, 1, 0, 1;
     0, 0, 0, 0, 1, 0;];
L = [2 -1  0  0  -1  0;
    -1  2  0  0  -1  0;
     0  0  1 -1   0  0;
     0  0 -1  2  -1  0;
    -1 -1  0 -1  4  -1;
     0  0  0  0  -1  1];

B = diag([1,0,1,0,0,0]);
H = L+B;
p = [1, 0, 1, 0, 0, 0];

%% initialize of FTPPC
rho_0_x = 10; rho_inf_x = 1.2; T_pre_x = 10; mu_x = 2;
rho_0_y = 10;  rho_inf_y = 1.2; T_pre_y = 10; mu_y = 2;
rho_0_z = 10;  rho_inf_z = 1.2; T_pre_z = 10; mu_z = 2;
rho_0=10;
Rho_t_x = 0;
Rho_t_y = 0;
Rho_t_z = 0;
rho_inf=1.2;
T_pre = 10;
mu=2;
rho_t=0;
%% initialize of parameters
e_normalized=0;
delta= 0.7 * ones(18, 1);
Gamma_i = 0;  % 实际 Gamma_i初值， 你可以根据具体情况调整
hat_omega_i0 = 0; % 观测的初值
hat_Gamma_i = 0;  % 观测的初值
omega_i0 = 0;
omega_im_previous = 0;
q = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
omega = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
domega = zeros(3,6);   
omega_im = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;];    
omega_i0 = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;];    
hat_omega_i0 = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; 
hat_Gamma_i= [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; 
a = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
hat_omega_i = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;];  
   dhat_a = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
   
   dhat_a1 = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
hat_a = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
hat_a1 = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
do = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
PI=[0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
 PI1=[0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; % states of the followers
   
 sumsliding = zeros(3, 6); 
constq = [0;0;0];
xilast_Eqi = zeros(3, 1);        
dxilast_Eqi = zeros(3, 1); 
xi_Eqi = zeros(3, 6); 
dxi_Eqi = zeros(3, 6);
ddxi_Eqi = zeros(3, 1);
I  = [1,0,0;
      0,1,0;
      0,0,1];
I6 =  ones(6, 6); 
HEQI = zeros(18,1);  
HDEQI = zeros(18,1);  
sumsliding = zeros(18,1);  
S1 = zeros(18,1);  
dS1 = zeros(18,1);  
S2 = zeros(18,1);  
TAOEQ = zeros(18,1);  
DTAOSW = zeros(18,1);  
TAOSW = zeros(18,1);  
TAO = zeros(18,1);  
col_Eai_ssum = zeros(18,1);  
col_Eai_sum = zeros(18,1);  
%% initialize of sliding mode
S1 = zeros(3, 6); % Nf 是 follower 数量
S2 = zeros(18, 1);
dS1 = zeros(3, 6);
ks = 10;
z1_all = zeros(3, 6);
z2_all = zeros(3, 6);
%% initialize of controller
taoeq = zeros(3, 6);
taosw = zeros(3, 6);
dtaosw = zeros(3, 6);
tao = zeros(3, Nf);

%% initialize of parameters
l1 = 0.1; % 等效控制增益
l2 = 0.1; % 等效控制增益
r1 = 1/3; % 滑模指数
r2 = 7/3; % 滑模指数
r3 = 2;
k3 = 1; % 切换控制增益

   gamma = 1000;
   integralterm = [0,0,0]';
   sigma = [0,0,0]';
   
 chi1 = [0,0,0]';
 chi2 = [0,0,0]';
 phi1 = [0,0,0]';
 phi2 = [0,0,0]';
   qqq=0.7;
   ppp=1.3;
   kkk1 = 1; % Set appropriate value
kkk2 = 1; % Set appropriate value
ggg1 = 1; % Set appropriate value
ggg2 = 1; % Set appropriate value
r0_1 = 2-qqq; % Set appropriate value
r0_2 = 1; % Set appropriate value
r_inf_1 = 2-ppp; % Set appropriate value
r_inf_2 = 1; % Set appropriate value
d0 = qqq-1; % Set appropriate value
d_inf = ppp-1; % Set appropriate value
epsilon = 0.4; % Set appropriate value

k0 = [10,0,0;0,5,0;0,0,3];      % 
g0 = [10,0,0;0,5,0;0,0,3];      

p0 = 0.75;      % can be ajusted according to the specific situation
q0 = 1.25;      % can be ajusted according to the specific situation

k1 = [16,0,0;0,10,0;0,0,6];      
g1 = [16,0,0;0,10,0;0,0,6];      

p1 = 0.75;      
q1 = 1.25;      

k2 = [20,0,0;0,10,0;0,0,10];      
g2 = [20,0,0;0,10,0;0,0,10];      

p2 = 0.75;      
q2 = 1.25;      
   
x0 = [0,0,0]';      % state of the leader
v0 = [0,0,0]';      % state of the leader
a0 = [0,0,0]';      % state of the leader

nu = [0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;
       0, 0, 0, 0, 0, 0;]; 
nu2_i = [0, 0, 0]'; 
tao = [0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;];

xi = [0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;];
Eai_sum =  [0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;];
Eai_ssum =  [0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;];
xi2_i = [0, 0, 0]';
flag1 = 0;
dt = 0.01;
tf = 60;
t_sample = [0:dt:tf];

t = 0;

%% parameters

eta11 = [5/4;0;0];
eta12 = [3/4;0;0];
eta13 = [1/4;0;0];
eta14 = [-1/4;0;0];
eta15 = [-3/4;0;0];
eta16 = [-5/4;0;0];

eta21 = [0;1/2;0];
eta22 = [1/2;0;0];
eta23 = [1;-1/2;0];
eta24 = [0;-0.5;0];
eta25 = [-1;-0.5;0];
eta26 = [-1/2;0;0];

eta31 = [-1;1/2;0];
eta32 = [0;1/2;0];
eta33 = [1;1/2;0];
eta34 = [1;-1/2;0];
eta35 = [0;-1/2;0];
eta36 = [-1;-1/2;0];

etas = [eta11, eta12, eta13, eta14, eta15, eta16;
         eta21, eta22, eta23, eta24, eta25, eta26;
         eta31, eta32, eta33, eta34, eta35, eta36
        ];

%model parameters
m11 = 25.8;
m22 = 24.661;    
m23 = 1.095;
m32 = 1;
m33 = 2.76;

%control parameters
alpha1 = 0.2;
alpha2 = 2*alpha1/(alpha1+1);
alpha3 = (4-3*alpha2)/(2-alpha2);
alpha4 = (4-3*alpha2)/(3-2*alpha2);

beta = 4;
c1 = 1;
c2 = 1;
c3 = 0.2;
c4 = 0.2;
c5 = 2;
qq=3;
 flag = 0;

aa1= [2,0,0;0,2,0;0,0,1];
aa2 = [1,0,0;0,2,0;0,0,1];

ro = [0,0,0]';
theta = 0;
miu = 8;
ro1 = 0;
sigma1 =0;
integralterm1 = 0;
 disturbance = [0,0,0]';
%% store the states for plot
Q = zeros(3,1,6);
Q(:,1,:) = q;
NU = zeros(3,1,6);
NU(:,1,:) = nu;
TAO = zeros(18,1);
Tao(:,1,:) = tao;
Do = zeros(3,1,6);
Do(:,1,:) = do; 
PII = zeros(3,1,6);
PII(:,1,:) = PI; 
Disturbance = [disturbance];
X0 = [x0];
formationx = [0];
formationy = [0];
formationtheta = [0];
formation = [0,0,0]';
omega_im_store = [0,0,0]'; 
hat_omega_i_store = [0,0,0]'; 

% tracking error
Eqi = [0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;]; 

% store the tracking error for plot
EQI = zeros(3,1,6);
EQI(:,1,:) = Eqi;
z_a = zeros(3,1);

XI_Eqi = zeros(18,1);
DXI_Eqi = zeros(18,1);
DDXI_Eqi = zeros(18,1);

R =  zeros(3,3,6);
%% main body
for k = 1:size(t_sample,2)-1 % time iteration

dxi_Eqi_last = dxi_Eqi;
xi_Eqi_last = xi_Eqi;
 dx0 = v0;
dv0 = a0;

%leader kinemics equation
    for i = 1:Nf % agent iteration
        z1=0; z2=0; zz1 = 0;z3=0;
        % follower kinemics equation parameters
   M = [m11, 0,   0; 
        0 ,  m22, m23;
        0,   m32, m33];

        C = [0,             0,       -m22*nu(2,i)-m23*nu(3,i);
              0,             0,       m11*nu(1,i);
       m22*nu(2,i)+m23*nu(3,i), -m11*nu(1,i),    0];

        d11 = 0.723 + 1.327 * abs(nu(1,i)) + 5.866*nu(1,i)^2;
        d22 = 0.861 + 36.282*abs(nu(2,i)) + 8.050*abs(nu(3,i));
        d23 = -0.108+ 0.845*abs(nu(2,i)) + 3.340*abs(nu(3,i));
        d32 = -0.105- 5.044*abs(nu(2,i)) - 0.130*abs(nu(3,i));
        d33 = 1.900 - 0.08*abs(nu(2,i)) + 0.75*abs(nu(3,i));
        
        D = [d11,  0,   0;
               0,  d22, d23;
               0,  d32, d33;];

        R(:,:,i) = [cos(q(3,i)), -sin(q(3,i)), 0;
             sin(q(3,i)), cos(q(3,i)), 0;
                   0,      0,        1;];  

        S = [0, -nu(3,i), 0;
            nu(3,i), 0, 0;
            0, 0, 0;];
 
disturbance = [3+3*sin(0.02*t)+2.5*sin(0.1*t);-4+2*sin(0.02*t)+0.5*sin(0.3*t);-3.3+2*sin(0.02*t)+2.5*sin(0.3*t)];
%         %follower kinemics equation 
         nu(:,i) = inv(R(:,:,i))*omega(:,i);
         xi(:,i) = R(:,:,i)*S*nu(:,i) - R(:,:,i)*inv(M)*(C*nu(:,i)+D*nu(:,i));
        domega(:,i) = R(:,:,i)*inv(M)*(tao(:,i))+xi(:,i)+disturbance;   
        omega(:,i) = omega(:,i) + domega(:,i)*dt;
        q(:,i) = q(:,i) + omega(:,i)*dt;

        if t<20
            kk=1;t_reset = t;
        elseif t>=20&&t<40
            kk=2;t_reset = t-20;
        else
            kk=3;t_reset = t-40;
        end
        
        % measurement error 
        x0 = [8*sin(0.1*t), -8*cos(0.1*t), 1*0.1*t]';
        v0 = [0.8*cos(0.1*t), 0.8*sin(0.1*t), 0.1]';
        a0 = [-0.08*sin(0.1*t), 0.08*cos(0.1*t), 0]';
        
        Eqi(:,i) = q(:,i)-etas((kk-1)*3+1:kk*3,i)-x0;
        omega_i = omega(:,i); 
        omega_im(:,i) = omega_i .* (1 + (rand(size(omega_i)) - 0.5) * 0.6);  %messurement noise 30%
        Eomegai(:,i) = omega_im(:,i) -v0;
         Eai(:,i) = a(:,i)-a0;
         Edomegai(:,i) = domega(:,i)-a0;
        EQI(:,k+1,:) = Eqi;

        for j = 1:Nf
            % auxiliary variables
            Eqij = q(:,i)-q(:,j)-etas((kk-1)*3+1:kk*3,i)+etas((kk-1)*3+1:kk*3,j); %i=6
            Eaij = a(:,i)-a(:,j);
            Eomegaij =  omega_im(:,i)- omega_im(:,j);
            Edomegaij =  domega(:,i)- domega(:,j);
            
            z1 = W(i,j)*Eqij+z1;
            z2 = W(i,j)*Eomegaij+z2;
            z3 = W(i,j)*Edomegaij+z3;
            
            zz1 = W(i,j)*Eaij + zz1;
        end
z1_all(:,i) = z1 + p(i) * Eqi(:,i);
z2_all(:,i) = z2 + p(i) * Eomegai(:,i);
z3_all(:,i) = z3 + p(i) *  Edomegai(:,i);

omega_i0(:,i) = omega_i0(:,i)+omega_im(:,i)*dt; 
e_omega0 = hat_omega_i0(:,i) - omega_i0(:,i);

dot_hat_omega_i0 = hat_omega_i(:,i) - epsilon * k0 * abs(e_omega0 / epsilon^2).^p0 .* sign(e_omega0 / epsilon^2) ...
                 - epsilon * g0 *abs(e_omega0 / epsilon^2).^q0 .* sign(e_omega0 / epsilon^2);
             
         nu2_i = inv(R(:,:,i))*hat_omega_i(:,i);
         xi2_i = R(:,:,i)*S*nu(:,i) - R(:,:,i)*inv(M)*(C*nu(:,i)+D*nu(:,i));    
         
dot_hat_omega_i = hat_Gamma_i(:,i) + R(:,:,i)*inv(M)*(tao(:,i)) +xi2_i - k1* abs(e_omega0 / epsilon^2).^p1.* sign(e_omega0 / epsilon^2) ...
                    - g1 * abs(e_omega0 / epsilon^2).^q1.* sign(e_omega0 / epsilon^2) ;
                
dot_hat_Gamma_i = -(k2 / epsilon) * abs(e_omega0 / epsilon^2).^p2 .* sign(e_omega0 / epsilon^2) ...
                    - (g2 / epsilon) * abs(e_omega0 / epsilon^2).^q2 .* sign(e_omega0 / epsilon^2);             

 hat_omega_i0(:,i) = hat_omega_i0(:,i) + dot_hat_omega_i0 * dt; 
hat_omega_i(:,i) = hat_omega_i(:,i) + dot_hat_omega_i * dt;
hat_Gamma_i(:,i) = hat_Gamma_i(:,i) + dot_hat_Gamma_i * dt;

do(:,i) = hat_Gamma_i(:,i);

       ddqri = a(:,i);      
      flag1 = R(:,:,i)*inv(M);

         dai = -beta * sign(zz1+p(i)*Eai(:,i)) - c5*sig( zz1+p(i)*Eai(:,i),qq);

         a(:,i) = a(:,i) + dai*dt;
omega_im_previous = omega(:,i);

        if mod(k,900)==0&&k>90
            formation = [formation, q(1:3,i)];
            if i==6
                formation = [formation, q(1:3,1)];
            end
            flag = 1; 
        end   
    end
%H, I  HEQI HDEQI sumsliding S1 dS1 S2 TAOEQ DTAOSW TAOSW TAO
HH = kron(H,I);
XI_Eqi_last = XI_Eqi;
DXI_Eqi_last = DXI_Eqi;

sig = @(z,alpha) sign(z).*abs(z).^alpha;
%xi(:,i),do(:,i)transformed to column
col_xi = xi(:);
col_do = do(:);
col_a = a(:);
% col_a = [a0;a0;a0;a0;a0;a0];
col_z1_all_x = z1_all(1,:);  % x 
col_z1_all_y = z1_all(2,:);  % y 
col_z1_all_z = z1_all(3,:);  % z 
n = length(col_z1_all_x);
col_z1_all = reshape([col_z1_all_x; col_z1_all_y; col_z1_all_z], [], n);
col_z1_all = col_z1_all(:);

col_z2_all_x = z2_all(1,:);  % x 
col_z2_all_y = z2_all(2,:);  % y 
col_z2_all_z = z2_all(3,:);  % z 
n = length(col_z2_all_x);
col_z2_all = reshape([col_z2_all_x; col_z2_all_y; col_z2_all_z], [], n);
col_z2_all = col_z2_all(:);

col_z3_all_x = z3_all(1,:);  % x 
col_z3_all_y = z3_all(2,:);  % y 
col_z3_all_z = z3_all(3,:);  % z 
col_z3_all = [col_z3_all_x(:);col_z3_all_y(:);col_z3_all_z(:)];
n = length(col_z3_all_x);
col_z3_all = reshape([col_z3_all_x; col_z3_all_y; col_z3_all_z], [], n);
col_z3_all = col_z3_all(:);
%FTPPC
if t_reset < T_pre_x
    rho_t_x = (rho_0_x - rho_inf_x) * (1 - t_reset / T_pre_x)^mu_x + rho_inf_x;
    rho_t_dot_x = -(rho_0_x - rho_inf_x)  * mu/T_pre * (1-t_reset/T_pre)^(mu-1);
    rho_t_ddot_x = (rho_0_x - rho_inf_x) * mu * (mu-1)/T_pre^2 * (1-t_reset/T_pre)^(mu-2);
else
    rho_t_x = rho_inf_x;
end

if t_reset < T_pre_y
    rho_t_y = (rho_0_y - rho_inf_y) * (1 - t_reset / T_pre_y)^mu_y + rho_inf_y;
    rho_t_dot_y = -(rho_0_y - rho_inf_y)  * mu/T_pre * (1-t_reset/T_pre)^(mu-1);
    rho_t_ddot_y = (rho_0_x - rho_inf_y) * mu * (mu-1)/T_pre^2 * (1-t_reset/T_pre)^(mu-2);
else
    rho_t_y = rho_inf_y;
end

if t_reset < T_pre_z
    rho_t_z = (rho_0_z - rho_inf_z) * (1 - t_reset / T_pre_z)^mu_z + rho_inf_z;
    rho_t_dot_z = -(rho_0_z - rho_inf_z)  * mu/T_pre * (1-t_reset/T_pre)^(mu-1);
    rho_t_ddot_z = (rho_0_z - rho_inf_z) * mu * (mu-1)/T_pre^2 * (1-t_reset/T_pre)^(mu-2);
else
    rho_t_z = rho_inf_z;
end

rho_t = concatenate_arrays(rho_t_x, rho_t_y, rho_t_z);
rho_t_dot = concatenate_arrays(rho_t_dot_x, rho_t_dot_y, rho_t_dot_z);
rho_t_ddot = concatenate_arrays(rho_t_ddot_x, rho_t_ddot_y, rho_t_ddot_z);

Rho_t_x =  [Rho_t_x, rho_t_x];       
Rho_t_y =  [Rho_t_y, rho_t_y];   
Rho_t_z =  [Rho_t_z, rho_t_z];   

e_normalized_x = col_z1_all_x / rho_t_x;
e_normalized_y = col_z1_all_y / rho_t_y;
e_normalized_z = col_z1_all_z / rho_t_z;
e_normalized = concatenate_arrays(e_normalized_x, e_normalized_y, e_normalized_z);
e_dot_normalized_x = col_z2_all_x./rho_t_x - col_z1_all_x.*rho_t_dot_x./rho_t_x.^2;  
e_dot_normalized_y = col_z2_all_y./rho_t_y - col_z1_all_y.*rho_t_dot_y./rho_t_y.^2; 
e_dot_normalized_z = col_z2_all_z./rho_t_z - col_z1_all_z.*rho_t_dot_z./rho_t_z.^2;  
e_dot_normalized = concatenate_arrays(e_dot_normalized_x, e_dot_normalized_y, e_dot_normalized_z);

e_ddot_normalized_x = (col_z3_all_x.*rho_t_x-2.*col_z2_all_x.*rho_t_dot_x-col_z1_all_x.*rho_t_ddot_x)/rho_t_x^2+2.*col_z1_all_x.*rho_t_dot_x.^2./rho_t_x^3;
e_ddot_normalized_y = (col_z3_all_y.*rho_t_y-2.*col_z2_all_y.*rho_t_dot_y-col_z1_all_y.*rho_t_ddot_y)/rho_t_y^2+2.*col_z1_all_y.*rho_t_dot_y.^2./rho_t_y^3;
e_ddot_normalized_z = (col_z3_all_z.*rho_t_z-2.*col_z2_all_z.*rho_t_dot_z-col_z1_all_z.*rho_t_ddot_z)/rho_t_z^2+2.*col_z1_all_z.*rho_t_dot_z.^2./rho_t_z^3;
e_ddot_normalized = concatenate_arrays(e_ddot_normalized_x, e_ddot_normalized_y, e_ddot_normalized_z);

%error transform
    H_function = 1./(delta+e_normalized)+1./(1-e_normalized);
    H_function_dot = -1./(delta+e_normalized).^2 + 1./(1-e_normalized).^2;
     XI_Eqi =  log((delta+e_normalized)./(delta.*(1-e_normalized)));

        DXI_Eqi = H_function.*e_dot_normalized;
        DDXI_Eqi =  H_function_dot.*e_dot_normalized.^2 + H_function.*e_ddot_normalized;

 DXI_Eqi_diff = (XI_Eqi-XI_Eqi_last)/dt;

 DDXI_Eqi_diff = (DXI_Eqi-DXI_Eqi_last)/dt;

sumsliding = sumsliding...
    -(-c1*sig(XI_Eqi,1/5)-c2*sig(DXI_Eqi,1/3)-c3*sig(XI_Eqi,9/5)-c4*sig(DXI_Eqi,9/7))*dt; %3N*1
    
S1 = DXI_Eqi+sumsliding;

dS1 = - l1 * abs(S1).^r1 .* sign(S1) ...
           - l2 * abs(S1).^r2 .* sign(S1);

S2 = dS1 + l1 .* abs(S1).^r1 .* sign(S1) ...
          + l2 .* abs(S1).^r2 .* sign(S1);
rho_t_vector = [rho_t_x * ones(size(col_z1_all_x(:))); 
                rho_t_y * ones(size(col_z1_all_y(:))); 
                rho_t_z * ones(size(col_z1_all_z(:)))];
rho_t_dot_vector = [rho_t_dot_x * ones(size(col_z1_all_x(:))); 
                rho_t_dot_y * ones(size(col_z1_all_y(:))); 
                rho_t_dot_z * ones(size(col_z1_all_z(:)))];
rho_t_ddot_vector = [rho_t_ddot_x * ones(size(col_z1_all_x(:))); 
                rho_t_ddot_y * ones(size(col_z1_all_y(:))); 
                rho_t_ddot_z * ones(size(col_z1_all_z(:)))];            
TAOEQ = col_a - col_xi - col_do ...
             +inv(HH)*(((- l1 .* abs(S1).^r1 .* sign(S1) ...
             - l2 .* abs(S1).^r2 .* sign(S1)...
            +-c1 .* sig(XI_Eqi,1/5)-c2 .* sig(DXI_Eqi,1/3)-c3 .* sig(XI_Eqi,9/5)-c4 .* sig(DXI_Eqi,9/7)).*rho_t_vector...
        -H_function_dot .* e_dot_normalized.^2.*rho_t_vector) .* 1./H_function + (2.*col_z2_all .* rho_t_dot_vector+col_z1_all .* rho_t_ddot_vector)./rho_t_vector...
            -2.*col_z1_all.*(rho_t_dot_vector).^2./rho_t_vector.^2);
 
DTAOSW = -k3 *abs(S2).^r3.* tanh(ks*S2);
TAOSW = TAOSW + DTAOSW * dt;
R_matrix = blkdiag(R(:,:,1),R(:,:,2),R(:,:,3),R(:,:,4),R(:,:,5),R(:,:,6));
M_matrix = blkdiag(M,M,M,M,M,M);
MR_inverse = M_matrix * inv(R_matrix);

TAO =  MR_inverse*(TAOEQ+TAOSW);
tao = reshape(TAO, 3 ,6);
 
    x0 = x0 + dx0*dt;
    v0 = v0 + dv0*dt;
    t = t + dt;
    Q(:,k+1,:) = q;
    Do(:,k+1,:) = do;
    PII(:,k+1,:) = PI;
    Tao(:,k+1,:) = tao;
    NU(:,k+1,:) = nu;
    Disturbance = [Disturbance,disturbance];
    X0 = [X0, x0];
    
    omega_im_store(:, k+1) = omega_im(:,i);
    hat_omega_i_store(:, k+1) = hat_omega_i(:,i);
    
end

save('Q.mat','Q');
save('EQI.mat','EQI'); 
save('do_fixed_time_data.mat', 'Do');
save('dis.mat', 'Disturbance');


figure (1);
subplot(3,1,1); 
plot(t_sample, omega_im_store(1,:),  'LineWidth', 2); hold on;
plot(t_sample, hat_omega_i_store(1,:,1), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Velocity', 'Interpreter', 'latex');
legend('$\omega_{im}$', '${\hat{\omega}}_i$', 'Interpreter', 'latex');

subplot(3,1,2);
plot(t_sample, omega_im_store(2,:), 'LineWidth', 2); hold on;
plot(t_sample, hat_omega_i_store(2,:,1),  'LineWidth', 2);
xlabel('Time (s)');
ylabel('Velocity', 'Interpreter', 'latex');
legend('$\omega_{im}$', '${\hat{\omega}}_i$', 'Interpreter', 'latex');

subplot(3,1,3);
plot(t_sample, omega_im_store(3,:),  'LineWidth', 2); hold on;
plot(t_sample, hat_omega_i_store(3,:),  'LineWidth', 2);
xlabel('Time (s)');
ylabel('Velocity', 'Interpreter', 'latex');
legend('$\omega_{im}$', '${\hat{\omega}}_i$', 'Interpreter', 'latex');

figure (2);

legend_width = 0.1;
legend_height = 0.6;

left_position = 0.1;
legend_position = left_position + 1 - legend_width;

pos1 = [left_position, 0.7, 0.6, 0.25];
pos2 = [left_position, 0.4, 0.6, 0.25];
pos3 = [left_position, 0.1, 0.6, 0.25];

ax1 = subplot('Position', pos1);
plot(t_sample, EQI(1,:,1), '--', t_sample, EQI(1,:,2), ':',...
    t_sample, EQI(1,:,3), '-.', t_sample, EQI(1,:,4), '-',...
    t_sample, EQI(1,:,5), '--', t_sample, EQI(1,:,6), '-.',...
    t_sample, Rho_t_x,'k', t_sample, -Rho_t_x,'k',...，
    'LineWidth', 2.5);
ylabel(ax1, {'x_{e}(m) '});
lgd1 = legend(ax1, 'i = 1', 'i = 2', 'i = 3', 'i = 4', 'i = 5', 'i = 6','boundary');
set(lgd1, 'Position', [legend_position, 0.1, legend_width, legend_height], 'NumColumns', 1);

ax2 = subplot('Position', pos2);
plot(t_sample, EQI(2,:,1), '--', t_sample, EQI(2,:,2), ':',...
    t_sample, EQI(2,:,3), '-.', t_sample, EQI(2,:,4), '-',...
    t_sample, EQI(2,:,5), '--', t_sample, EQI(2,:,6), '-.',...
    t_sample, Rho_t_y,'k',t_sample, -Rho_t_y,'k',...
    'LineWidth', 2.5);
ylabel(ax2, {'y_{e}(m) '});
lgd2 = legend(ax2, 'i = 1', 'i = 2', 'i = 3', 'i = 4', 'i = 5', 'i = 6','boundary');
set(lgd2, 'Position', [legend_position, 0.1, legend_width, legend_height], 'NumColumns', 1);

ax3 = subplot('Position', pos3);
plot(t_sample, EQI(3,:,1), '--', t_sample, EQI(3,:,2), ':',...
    t_sample, EQI(3,:,3), '-.', t_sample, EQI(3,:,4), '-',...
    t_sample, EQI(3,:,5), '--', t_sample, EQI(3,:,6), '-.',...
    t_sample, Rho_t_z,'k',t_sample, -Rho_t_z,'k',...
    'LineWidth', 2.5);
ylabel(ax3, {'r_{e}(rap) '});
xlabel(ax3, 'Times (s)');
lgd3 = legend(ax3,'i = 1', 'i = 2', 'i = 3', 'i = 4', 'i = 5', 'i = 6', 'boundary');
set(lgd3, 'Position', [legend_position, 0.1, legend_width, legend_height], 'NumColumns', 1);

figure (3)
plot(X0(1,:),X0(2,:),'LineWidth',1);
plot(Q(1,:,1),Q(2,:,1),Q(1,:,2),Q(2,:,2),Q(1,:,3),Q(2,:,3),Q(1,:,4),Q(2,:,4),Q(1,:,5),Q(2,:,5),Q(1,:,6),Q(2,:,6),'LineWidth',1);

hold on
for i=1:(length(formation)-1)/7 

    formationx = formation(1,:); 
    formationy = formation(2,:);
    formationtheta = formation(3,:);

    for ii=1:length(formationx)
    AirplanePlot([formationx(ii+1),formationy(ii+1)],formationtheta(ii+1)-pi/2,0.5); 
    end
    hold on;
end
