%% 模型初始化
clear all;
clc;

%% 初始化参数
% 初始位置、速度
R0 = [36 748 0];
V0 = [-8.6 -69.6 0];

% 目标位置、速度
Rf = [0 0 0];
Vf = [0 0 0];

% 初始欧拉角、欧拉角速度
Euler0 = [0 0 -3]/57.3;% 57.3=180/Π为度转弧度
Omega0 = [0 0 0]/57.3;

% 目标欧拉角、欧拉角速度
Eulerf = [0 0 90]/57.3;
Omegaf = [0 0 0]/57.3;

% 初始火箭质量
m0 = 155975;

% 初始状态向量
xx0 = [R0 V0 m0];

% 求解结束判断条件
options = odeset('Events',@(t,x)odeEventFun(t,x));
% 四次多项式求解初始轨迹和终端时间
[t_sj,xx] = ode45(@Polynomial_Guidance,0:1:100,xx0,options);

%% 仿真参数
cvx_status='Infeasible';

mwet = m0;% 火箭总质量
mfuel = 35352;% 剩余燃料质量
mdry = mwet-mfuel;% 火箭干重量

Pe = 2450000;% 单台发动机额定推力
dG = 580;% 单台发动机秒耗量
g0 = 9.8;
Isp = Pe/(dG*g0);% 单台发动机比冲

Tmax = Pe*1.05;% 单台发动机最大推力
Tmin = Pe*0.5;% 单台发动机最小推力
Smax = 8/57.3;% 单台发动机最大摆角

a = -1/(Isp*g0); % 变量替换系数α
Sref = 63.617;% 气动参考面积
rho = 1.2;% 大气密度
% beita=1.3785e-4;
count=0;

tf = ceil(t_sj(end))-1;% 终端时间取整
dt = 1;
N = tf/dt+1;% 离散点个数
t_zd(1,1) = 0;
for i = 2:N
    t_zd(i,1) = 0+dt*(i-1);
end

% 落点系与箭体系转换矩阵
cphi = cos(Euler0(3)); sphi = sin(Euler0(3));
cpsi = cos(Euler0(2)); spsi = sin(Euler0(2));
cgam = cos(Euler0(1)); sgam = sin(Euler0(1));

Cbody_land = [ cphi*cpsi, sphi*cpsi, -sphi;
               cphi*spsi*sgam-sphi*cgam, sphi*sgam*spsi+cphi*cgam, cpsi*sgam;
               cphi*spsi*cgam+sphi*sgam, sphi*cgam*spsi-cphi*sgam, cpsi*cgam ];

Vbody = Cbody_land*V0;% 箭体系下速度
alpha = tan(Vbody(3)/Vbody(1));% 攻角

% 气动力参数
[Cx, Cn, Cm, Xcp] = fcn(0,0,alpha);

% 动力学方程输入矩阵
B = [0 0 0 0 0 0 0
     0 0 0 0 0 0 0
     0 0 0 0 0 0 0
     1 0 0 1 0 0 0
     0 1 0 0 1 0 0
     0 0 1 0 0 1 0
     0 0 0 0 0 0 -a];

%% 凸优化迭代更新参考轨迹
while 1
    % 提取初始轨迹（上一时刻轨迹）
    rx = xx(:,1);
    ry = xx(:,2);
    rz = xx(:,3);
    vx = xx(:,4);
    vy = xx(:,5);
    vz = xx(:,6);
    z = log(xx(:,7));

    % 凸优化问题求解
    cvx_solver SeDuMi
    cvx_begin
    variables PKX(N) PKY(N) PKZ(N)   % 无损凸化T(t)/m(t)
    variable QK(N)                   % 无损凸化Γ(t)/m(t)
    variables RX(N) RY(N) RZ(N)      % 位置
    variables VX(N) VY(N) VZ(N)      % 速度
    variable Z(N)                    % 对数质量Z(k)=ln(m(t))
    minimize sum(QK*dt)              % 目标函数ΣσΔt
    subject to

    % 初始条件约束
    RX(1) == R0(1); 
    RY(1) == R0(2);
    RZ(1) == R0(3);
    VX(1) == V0(1); 
    VY(1) == V0(2);
    VZ(1) == V0(3);
    Z(1) == log(m0);

    % 构造采用最大推力导致的最小质量的对数Z0
    for i = 1:N
       Z0(i) = log(mwet-a*Tmax*i*dt);
    end

    % 离散动力学约束
    for i = 1:N-1
        % 利用初始轨迹（上一时刻轨迹）计算气动阻力参数
        f(i) = -0.5*Cx*Sref*rho*norm([vx(i) vy(i) vz(i)])/m(i);

        % 离散方程状态转移矩阵
        G = [1 0 0 (exp(f(i)*dt)-1)/f(i) 0 0 0
             0 1 0 0 (exp(f(i)*dt)-1)/f(i) 0 0
             0 0 1 0 0 (exp(f(i)*dt)-1)/f(i) 0
             0 0 0 exp(f(i)*dt) 0 0 0
             0 0 0 0 exp(f(i)*dt) 0 0
             0 0 0 0 0 exp(f(i)*dt) 0
             0 0 0 0 0 0 1]; 

        % 状态转移矩阵积分
        G_int = [dt 0 0 (exp(f(i)*dt)-1)/(f(i)^2)-dt/f(i) 0 0 0
                 0 dt 0 0 (exp(f(i)*dt)-1)/(f(i)^2)-dt/f(i) 0 0
                 0 0 dt 0 0 (exp(f(i)*dt)-1)/(f(i)^2)-dt/f(i) 0
                 0 0 0 (exp(f(i)*dt)-1)/f(i) 0 0 0
                 0 0 0 0 (exp(f(i)*dt)-1)/f(i) 0 0
                 0 0 0 0 0 (exp(f(i)*dt)-1)/f(i) 0
                 0 0 0 0 0 0 dt];
        H = G_int*B;% 离散方程输入矩阵

        % 离散方程
        [RX(i+1);RY(i+1);RZ(i+1);VX(i+1);VY(i+1);VZ(i+1);Z(i+1)] == ...
        G*[RX(i);RY(i);RZ(i);VX(i);VY(i);VZ(i);Z(i)]+H*[PKX(i);PKY(i);PKZ(i);0;-g0;0;QK(i)];

        % 映入松弛变量
        norm([PKX(i) PKY(i) PKZ(i)]) <= QK(i);

        % 推力下界
        QK(i) >= Tmin*exp(-Z0(i))*(1-(Z(i)-Z0(i))+((Z(i)-Z0(i))^2)/2);
        % 推力上界
        QK(i) <= Tmax*exp(-Z0(i))*(1-(Z(i)-Z0(i)));

        % z的范围约束
        Z0(i) <= Z(i) <= log(mwet-a*Tmin*i*dt);
    end

    % 终端时刻约束
    norm([PKX(N) PKY(N) PKZ(N)]) <= QK(N);
    QK(N) >= Tmin*exp(-Z0(N))*(1-(Z(N)-Z0(N))+((Z(N)-Z0(N))^2)/2);
    QK(N) <= Tmax*exp(-Z0(N))*(1-(Z(N)-Z0(N)));
    PKX(N) == 0;PKZ(N) == 0;
    Z0(N) <= Z(N) <= log(mwet-a*Tmin*N*dt);
    RX(N) == 0;RY(N) == 0;RZ(N) == 0;
    VX(N) == 0;VY(N) == 0;VZ(N) == 0;
    Z(N) >= log(mdry);
    
    cvx_end

    % 推力结构转换
    for i = 1:N
        M(i,1) = exp(Z(i));% 质量结构转换
        TX(i,1) = PKX(i)*M(i);
        TY(i,1) = PKY(i)*M(i);
        TZ(i,1) = PKZ(i)*M(i);
    end

    % 优化得到的第1个时刻作为初始状态x0
    xx = [RX(1) RY(1) RZ(1) VX(1) VY(1) VZ(1) M(1)];
    x0 = [RX(1) RY(1) RZ(1) VX(1) VY(1) VZ(1) M(1)];

    % 将实际动力学按每个时间步积分
    for i = 1:N-1
        tspan = [t_zd(i) t_zd(i+1)];
        Tx = TX(i);Ty = TY(i);Tz = TZ(i);
        PKx = PKX(i);PKy = PKY(i);PKz = PKZ(i);
        f1 = f(i);
        [tt,x] = ode45(@zd,tspan,x0,[],Tx,Ty,Tz,PKx,PKy,PKz,f1);
        xx = [xx;x(end,:)];
        x0 = x(end,:);
    end

    count=count+1;
    if count>=8
        break
    end
end

%% 四次多项式求解初始轨迹和终端时间
function dx = Polynomial_Guidance(t,x)
    dx = zeros(7,1);
    rx = x(1);
    ry = x(2);
    rz = x(3);
    vx = x(4);
    vy = x(5);
    vz = x(6);
    m = x(7);

    % 仿真参数
    Pe = 2450000;% 单台发动机额定推力
    dG = 580;% 单台发动机秒耗量
    g0 = 9.8;
    Isp = Pe/(dG*g0);% 单台发动机比冲
    
    rfx = 0; rfy = 0; rfz = 0;   % 期望终端位置矢量
    vfx = 0; vfy = 0; vfz = 0;   % 期望终端速度矢量
    afx = 0; afy = 9.8; afz = 0; % 期望终端加速度矢量

    vf = [vfx vfy vfz];
    v0 = [vx vy vz];
    rf = [rfx rfy rfz];
    r0 = [rx ry rz];

    % 估算终端时间tg
    p = [g^2 0 -4*(norm(v0+vf)^2+(v0.')*vf) 24*(((rf-r0).')*(v0+vf)) -36*((rf-r0).')*(rf-r0)];
    r = roots(p);
    r = r(r == real(r)); % 取实根
    tg = r(r > 0);       % 取正根

    % 求解四次多项式得实时指令加速度矢量a0
    a0x = afx + 12*(rfx - rx)/(tg^2) - 6*(vfx + vx)/tg;
    a0y = afy + 12*(rfy - ry)/(tg^2) - 6*(vfy + vy)/tg;
    a0z = afz + 12*(rfz - rz)/(tg^2) - 6*(vfz + vz)/tg;

    % 推力
    Tx = a0x*m;
    Ty = a0y*m;
    Tz = a0z*m;

    % 返回微分方程
    dx(1) = x(4);
    dx(2) = x(5);
    dx(3) = x(6);
    dx(4) = ax;
    dx(5) = ay;
    dx(6) = az;
    dx(7) = -sqrt(Tx^2+Ty^2+Tz^2)/(Isp*g0);
end

%% 求解停止条件
function [value,isterminal,direction] = odeEventFun(t,x)
    value = x(2)-0.0001; % 高度阈值 H1
    isterminal = 1;      % 触发后终止积分
    direction = -1;      % 仅在 value 从正变负时触发
end

%% 
function dy=zd(t,y,Tx,Ty,Tz,PKX,PKY,PKZ,f1)
g0=9.8;
Isp=282;
CD=2.2;
Sref=10;
rho0=1.225;
beita=1.3785e-4;
dy=zeros(7,1);
rx=y(1);
ry=y(2);
rz=y(3);
vx=y(4);
vy=y(5);
vz=y(6);
m=y(7);
rho=rho0*exp(-beita*ry);
% adx=-CD*Sref*rho*norm([vx vy vz])*vx/(2*m);
% ady=-CD*Sref*rho*norm([vx vy vz])*vy/(2*m);
% adz=-CD*Sref*rho*norm([vx vy vz])*vz/(2*m);
% adx=0;
% ady=0;
% adz=0;
adx=f1*vx;
ady=f1*vy;
adz=f1*vz;
dy(1)=y(4);
dy(2)=y(5);
dy(3)=y(6);
dy(4)=0+Tx/m+adx;
dy(5)=-g0+Ty/m+ady;
dy(6)=0+Tz/m+adz;
% dy(4)=0+PKX+adx;
% dy(5)=-g0+PKY+ady;
% dy(6)=0+PKZ+adz;
dy(7)=-norm([Tx Ty Tz])/(Isp*g0);
end

%% 气动力参数计算
function [Cx, Cn, Cm, Xcp] = fcn(Forward, Rudder, alpha)
% 输入参数：
%   Forward：前舵舵角
%   Rudder：后舵舵角
%   alpha：攻角
% 输出参数：
%   总的气动参数 4*1 


% 数据拟合
%================= 基础数据 =================
x = [0,4,8,30,50,60,70,90,150]/57.3;   % 攻角
y = [0,20,40,60,80]/57.3;              % 舵角

% ----------- 前舵数据 -----------------
z_f = [0.2983	0.3105	0.3109	0.34	0.3102	0.36	0.3968	0.3115	-1.178;
0.0333	0.6536	1.442	5.3229	7.9765	8.5087	9.1668	9.894	3.4019;
0.0004	0.3464	0.8089	3.1887	4.5519	4.9497	5.3774	6.0923	2.8019;
0.0126	0.53	0.5609	0.599	0.5707	0.5817	0.5866	0.6158	0.8236;
0.2965	0.3117	0.3068	0.3731	0.3353	0.3577	0.4056	0.3175	-1.1556;
0.0433	0.6561	1.4587	5.4398	7.9116	8.4858	9.1227	9.8364	3.4092;
-0.0002	0.3129	0.7542	3.1904	4.5143	4.9434	5.3699	6.0813	2.8162;
-0.0044	0.4769	0.517	0.5865	0.5706	0.5826	0.5886	0.6182	0.8261;
0.2922	0.3164	0.3073	0.3742	0.3142	0.3435	0.3922	0.3478	-1.1666;
0.0595	0.6854	1.5233	5.6129	8.1435	8.7701	9.4561	10.1558	3.6688;
0.0169	0.3286	0.7791	3.2016	4.5616	4.985	5.4273	6.1567	2.9435;
0.2844	0.4794	0.5114	0.5704	0.5602	0.5684	0.5739	0.6062	0.8023;
0.2862	0.3039	0.2954	0.3019	0.3771	0.3267	0.4269	0.3565	-1.1929;
0.0006	0.6651	1.5116	6.0214	8.2621	9.098	9.5941	10.405	3.7881;
-0.0185	0.3289	0.7783	3.4086	4.5714	5.0651	5.4516	6.2152	2.9431;
0	0.4945	0.5149	0.5661	0.5533	0.5567	0.5682	0.5973	0.7769;
0.242	0.2554	0.2475	0.2483	0.3807	0.3942	0.4246	0.3628	-1.1474;
-0.0001	0.6199	1.4786	6.4426	8.5527	9.1594	9.8809	10.7777	4.0663;
-0.0002	0.29	0.7364	3.7257	4.6952	5.0971	5.57	6.3597	3.0003;
1.9663	0.4678	0.4981	0.5783	0.549	0.5565	0.5637	0.5901	0.7378;
];

% ----------- 后舵数据 -----------------
z_r = [0.242	0.2554	0.2475	0.2483	0.3807	0.3942	0.4246	0.3628	-1.1474;
-0.0001	0.6199	1.4786	6.4426	8.5527	9.1594	9.8809	10.7777	4.0663;
-0.0002	0.29	0.7364	3.7257	4.6952	5.0971	5.57	6.3597	3.0003;
1.9663	0.4678	0.4981	0.5783	0.549	0.5565	0.5637	0.5901	0.7378;
0.2966	0.3011	0.295	0.2658	0.3711	0.3814	0.4245	0.3793	-1.2536;
-0.0292	0.5531	1.3313	5.5525	7.9164	8.7014	9.3596	10.3076	3.832;
-0.023	0.242	0.6264	3.0118	4.2753	4.7378	5.2039	6.0424	2.7939;
0.787	0.4375	0.4705	0.5424	0.5401	0.5445	0.556	0.5862	0.7291;
0.2845	0.3062	0.2959	0.213	0.3556	0.3813	0.4158	0.4182	-1.2058;
0.0027	0.5215	1.2005	5.2407	7.2245	7.8964	8.6652	9.3737	3.6004;
0.0037	0.2172	0.5194	2.7472	3.675	4.096	4.6013	5.3534	2.5415;
1.3613	0.4164	0.4327	0.5242	0.5087	0.5187	0.531	0.5711	0.7059;
0.2868	0.2984	0.3035	0.1714	0.3074	0.2816	0.4045	0.4315	-1.1809;
0.0032	0.4744	1.0469	4.7343	6.759	7.2923	7.7609	8.5337	3.2481;
0.0054	0.179	0.3915	2.2943	3.2371	3.5137	3.8759	4.6227	2.1779;
1.7111	0.3774	0.3739	0.4846	0.4789	0.4818	0.4994	0.5417	0.6705;
0.2865	0.2988	0.3125	0.2126	0.3141	0.3044	0.3788	0.4087	-1.2033;
0.0065	0.4291	0.9333	4.0827	6.8615	7.3098	7.8867	8.5997	3.1052;
0.0094	0.1405	0.2985	1.7559	3.3187	3.5934	3.978	4.6746	2.0867;
1.435	0.3275	0.3198	0.4301	0.4837	0.4916	0.5044	0.5436	0.672
];

%================= 数据拆分 =================
Cx_fr  = flipud(z_f(1:4:end, :));
Cn_fr  = flipud(z_f(2:4:end, :));
Cm_fr  = flipud(z_f(3:4:end, :));
Xcp_fr = flipud(z_f(4:4:end, :));

Cx_rr  = z_r(1:4:end, :);
Cn_rr  = z_r(2:4:end, :);
Cm_rr  = z_r(3:4:end, :);
Xcp_rr = z_r(4:4:end, :);

% 数据插值结果
%================= 二维插值 =================
Cx_fr_  = interp2(x, y, Cx_fr,  alpha, Forward, 'spline');
Cn_fr_  = interp2(x, y, Cn_fr,  alpha, Forward, 'spline');
Cm_fr_  = interp2(x, y, Cm_fr,  alpha, Forward, 'spline');
Xcp_fr_ = interp2(x, y, Xcp_fr, alpha, Forward, 'spline');

Cx_rr_  = interp2(x, y, Cx_rr,  alpha, Rudder, 'spline');
Cn_rr_  = interp2(x, y, Cn_rr,  alpha, Rudder, 'spline');
Cm_rr_  = interp2(x, y, Cm_rr,  alpha, Rudder, 'spline');
Xcp_rr_ = interp2(x, y, Xcp_rr, alpha, Rudder, 'spline');

%================= 合成 =================
Cx  = Cx_fr_ + Cx_rr_;
Cn  = Cn_fr_ + Cn_rr_;
Cm  = Cm_fr_ + Cm_rr_;

if (Cn_fr_ + Cn_rr_) ~= 0
    Xcp = (Cn_fr_*Xcp_fr_ + Cn_rr_*Xcp_rr_) / (Cn_fr_ + Cn_rr_);
else
    Xcp = 0;
end
end