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
Euler0 = [-3 0 0]/57.3;% 57.3=180/Π为度转弧度
Omega0 = [0 0 0]/57.3;

% 目标欧拉角、欧拉角速度
Eulerf = [90 0 0]/57.3;
Omegaf = [0 0 0]/57.3;

% 初始火箭质量
m0 = 155975;

% 初始状态向量
xx0 = [R0 V0 Euler0 Omega0 m0];

% 求解结束判断条件
options = odeset('Events',@(t,x)odeEventFun(t,x));
% 四次多项式求解初始轨迹和终端时间
[t_sj,xx] = ode45(@Polynomial_Guidance,0:1:100,xx0,options);



%% 仿真参数
cvx_status='Infeasible';

mwet = m0;% 火箭总质量
mfuel = 35352;% 剩余燃料质量

Pe = 2450000;% 单台发动机额定推力
dG = 580;% 单台发动机秒耗量
g0 = 9.8;
Isp = Pe/(dG*g0);% 单台发动机比冲

Tmax = Pe*1.05;% 单台发动机最大推力
Tmin = Pe*0.5;% 单台发动机最小推力
Smax = 8/57.3;% 单台发动机最大摆角

% alpha=1/(Isp*g0); % 变量替换系数 alpha
% thetaTmax=60*pi/180; % 推力偏离 y 轴的最大角度（rad）
% CD = 2.2;
% Sref=10;
% rho0=1.225;
% beita=1.3785e-4;
% count=0;

tf = ceil(t_sj(end))-1;% 终端时间取整
dt = 1;
N = tf/dt+1;% 离散点个数
t_zd(1,1) = 0;
for i = 2:N
    t_zd(i,1) = 0+dt*(i-1);
end

% BC 矩阵：用于状态输入耦合（论文里称为 BC）
% BC=[0 0 0 0 0 0 0
%     0 0 0 0 0 0 0
%     0 0 0 0 0 0 0
%     1 0 0 1 0 0 0
%     0 1 0 0 1 0 0
%     0 0 1 0 0 1 0
%     0 0 0 0 0 0 -alpha];

% 对数质量的参考（用于某些近似）
Z=log(xx(:,7));

% tic

%% 凸优化迭代更新参考轨迹
while 1
    % 提取初始轨迹（上一时刻轨迹）
    rx = xx(:,1);
    ry = xx(:,2);
    rz = xx(:,3);
    vx = xx(:,4);
    vy = xx(:,5);
    vz = xx(:,6);
    m = xx(:,7);
    mp = log(m); % 取自然对数（向量）

    % 利用初始轨迹（上一时刻轨迹）计算阻力相关参数
    for i=1:N-1
       rho(i)=rho0*exp(-beita*ry(i));   % 随高度变化的大气密度
       % f(i) = - (1/2) * CD * Sref * rho * ||v|| / m
       % 这是阻力对速度的一种“线性化因子”（见论文 A 矩阵项）
       f(i)=-1/2*CD*Sref*rho(i)*norm([vx(i) vy(i) vz(i)])/m(i); 
       fn2(i)=exp(mp(i+1)-mp(i));      % exp(log(m_{i+1}) - log(m_i)) = m_{i+1}/m_i
       fn3(i)=exp(-mp(i));             % exp(-log(m_i)) = 1/m_i
       % 注：注释中有另一种 f(i) 定义被注释掉（使用 Z）
    end

    % 凸优化问题求解
    cvx_solver SeDuMi
    cvx_begin
    variables PKX(N) PKY(N) PKZ(N)   % 归一化推力方向分量（u = T/m 的分量）
    variable QK(N)                   % 归一化推力幅值 σ (论文记号)
    variables RX(N) RY(N) RZ(N)      % 离散化的状态位置
    variables VX(N) VY(N) VZ(N)      % 离散化的速度
    variable Z(N)                    % 对数质量 z(k) = ln(m(k))
    minimize sum(QK*dt)              % 目标：最小化 Σ σ Δt 即燃料消耗近似
    subject to
        % 初始条件约束（与 yy0 对应）
        RX(1)==3000; RY(1)==4500;
        VX(1)==-150; VY(1)==-320;
        Z(1)==log(22200+411000*0.03);
        RZ(1)==2600; VZ(1)==-260;
    % PKX(1)>=0; % 被注释的可能冗余约束

    % 构造 Z0（数值参考质量对数轨迹）——使用最大推力导致的最小质量近似
    for i=1:N
       Z0(i)=log(mwet-alpha*Tmax*i*dt);
    end

    % 逐步写离散动力学约束（用前面计算的 f(i) 代入）
    for i=1:N-1
        % AD 是状态转移矩阵（基于 f(i) 的指数形式），这里已按论文推导填写
        AD=[1 0 0 (exp(f(i)*dt)-1)/f(i) 0 0 0
           0 1 0 0 (exp(f(i)*dt)-1)/f(i) 0 0
           0 0 1 0 0 (exp(f(i)*dt)-1)/f(i) 0
           0 0 0 exp(f(i)*dt) 0 0 0
           0 0 0 0 exp(f(i)*dt) 0 0
           0 0 0 0 0 exp(f(i)*dt) 0
           0 0 0 0 0 0 1]; 
       % AINT 用于计算输入积分项 BD = AINT * BC
       AINT=[dt 0 0 (exp(f(i)*dt)-1)/(f(i)^2)-dt/f(i) 0 0 0
             0 dt 0 0 (exp(f(i)*dt)-1)/(f(i)^2)-dt/f(i) 0 0
             0 0 dt 0 0 (exp(f(i)*dt)-1)/(f(i)^2)-dt/f(i) 0
             0 0 0 (exp(f(i)*dt)-1)/f(i) 0 0 0
             0 0 0 0 (exp(f(i)*dt)-1)/f(i) 0 0
             0 0 0 0 0 (exp(f(i)*dt)-1)/f(i) 0
             0 0 0 0 0 0 dt];
       BD=AINT*BC;

       % 状态转移（等式约束）：
       [RX(i+1);RY(i+1);RZ(i+1);VX(i+1);VY(i+1);VZ(i+1);Z(i+1)] == ...
           AD*[RX(i);RY(i);RZ(i);VX(i);VY(i);VZ(i);Z(i)] + BD*[PKX(i);PKY(i);PKZ(i);0;-g0;0;QK(i)];

       % 下面是一些被注释的等价差分形式，可作为调试参考
%        [RX(i+1);RY(i+1);RZ(i+1)]==[RX(i);RY(i);RZ(i)]+dt*([VX(i);VY(i);VZ(i)]);
%        [VX(i+1);VY(i+1);VZ(i+1)]==[VX(i);VY(i);VZ(i)]+dt*(f(i)*[VX(i);VY(i);VZ(i)]+[PKX(i);PKY(i);PKZ(i)]+[0;-g0;0]);
%        Z(i+1)==Z(i)+dt*(-alpha*QK(i));

       % 力分量与 σ 的关系与角度约束
%        PKX(i)+PKY(i)<=QK(i);
%        (PKX(i))<=QK(i);
%       0<=PKX(i)<=10;

       PKY(i)>=0;
       % 二阶锥约束：||[PKX PKY PKZ]|| <= QK
       norm([PKX(i) PKY(i) PKZ(i)]) <= QK(i);

       % 推力下界（使用泰勒展开的近似，参考论文公式）
       QK(i) >= Tmin*exp(-Z0(i))*(1-(Z(i)-Z0(i))+((Z(i)-Z0(i))^2)/2);
       % 推力上界（线性化近似）
       QK(i) <= Tmax*exp(-Z0(i))*(1-(Z(i)-Z0(i)));

       RY(i) >= 0;
       % z 的范围约束（参考质量下界与上界）
       Z0(i) <= Z(i) <= log(mwet-alpha*Tmin*i*dt);  

       % 推力方向角约束： PKY >= QK * cos(thetaTmax)
       PKY(i) >= QK(i)*cos(thetaTmax);
    end

    % 速率和平滑性约束（相邻时间步的推力变化率约束等）
    for i=1:N-1
       -3 <= PKX(i+1)-PKX(i) <= 3;
       -5 <= PKY(i+1)-PKY(i) <= 5;
       -3 <= PKZ(i+1)-PKZ(i) <= 3;
       % 注意：fn2, fn3 之前已计算，是数字因子 (m_{k+1}/m_k) 和 1/m_k
       abs(QK(i+1)*fn2(i)-QK(i)) <= deltaTmax*fn3(i);
    end

    % 终端时刻的一些约束
    norm([PKX(N) PKY(N) PKZ(N)]) <= QK(N);
    QK(N) >= Tmin*exp(-Z0(N))*(1-(Z(N)-Z0(N))+((Z(N)-Z0(N))^2)/2);
    QK(N) <= Tmax*exp(-Z0(N))*(1-(Z(N)-Z0(N)));
    PKY(N) >= QK(N)*cos(thetaTmax);
    abs(QK(N)*fn2(N-1)-QK(N-1)) <= deltaTmax*fn3(N-1)
    PKX(N) == 0;
    PKZ(N) == 0;
    Z0(N) <= Z(N) <= log(mwet-alpha*Tmin*N*dt);
    RX(N) == 0; RY(N) == 0;
    VX(N) == 0; VY(N) == 0;
    Z(N) >= log(22200);
    VZ(N) == 0; RZ(N) == 0;
    % cvx_end 下文
    cvx_end

    % 将优化结果中的对数质量恢复为实际质量
    for i=1:N
        M(i,1)=exp(Z(i));
    end

    % 将归一化推力 u 转换为实际推力 T = u * m
    for i=1:N
        TX(i,1)=PKX(i)*M(i);
        TY(i,1)=PKY(i)*M(i);
        TZ(i,1)=PKZ(i)*M(i);
    end

    % 用优化得到的第 1 个时刻作为初始状态 y0
    xx=[RX(1) RY(1) RZ(1) VX(1) VY(1) VZ(1) M(1)];
    t1=0;
    y0=[RX(1) RY(1) RZ(1) VX(1) VY(1) VZ(1) M(1)];

    % 用 ode45 将实际动力学按每个时间步积分（或仿真真实动力学）
    for i=1:N-1
        tspan=[t_zd(i) t_zd(i+1)];
        Tx=TX(i);
        Ty=TY(i);
        Tz=TZ(i);
        PKx=PKX(i);
        PKy=PKY(i);
        PKz=PKZ(i);
        f1=f(i); % 阻力线性化因子
        [tt,y]=ode45(@zd,tspan,y0,[],Tx,Ty,Tz,PKx,PKy,PKz,f1);
        t1=[t1;tt(end,:)];
        xx=[xx;y(end,:)];
        y0=y(end,:);
    end

    count=count+1;
    if count>=8
        break
    end
end

%% 四次多项式求解初始轨迹和终端时间
function dx = Polynomial_Guidance(t,x)
    dx = zeros(13,1);
    rx = x(1);ry = x(2);rz = x(3);
    vx = x(4);vy = x(5);vz = x(6);
    m = x(13);

    % 仿真参数
    Pm = 10.2;
    Isp = 282;% Isp=
    g = -9.8; % 注意：这里把 g 定为负值（向下），与主程序中 g0=9.8 的定义需注意符号一致性
    urx=0; ury=0; urz=0;   % 航向修正项
    rfx=0; rfy=0; rfz=0;   % 目标位置（四次多项式文中通常是 rf）
    vfx=0; vfy=0; vfz=0;   % 目标速度
    afx=0; afy=9.8; afz=0; % 终端加速度指令（注意这里 afy = 9.8 用意？）
    tao=0.01;

    % 求解四次多项式剩余时间：解一个四次多项式（按照论文推导）
    vf2=[vfx vfy vfz];
    vf1=vf2*vf2.'; % ||v_f||^2
    v2=[vx vy vz];
    v1=v2*v2.';    % ||v||^2
    rf2=[rfx rfy rfz];
    r2=[rx ry rz];

    % p 为多项式系数，见论文中 H*(tg) 的展开
    p=[g^2 0 -4*(v1+vf1+vf2*(v2.')) 24*((vf2+v2)*((rf2-r2).')) -36*(rf2-r2)*((rf2-r2).')];
    r=roots(p);
    r=r(r==real(r)); % 取实根
    tg=r(r>0);       % 取正根（可能有多个，这里默认取所有正实根）

    % 计算四次多项式制导得到的加速度指令（见论文公式）
    upx=afx+12*(rfx-rx)/(tg^2)-6*(vfx+vx)/tg-0;
    upy=afy+12*(rfy-ry)/(tg^2)-6*(vfy+vy)/tg-g;
    upz=afz+12*(rfz-rz)/(tg^2)-6*(vfz+vz)/tg-0;
    % 若 tg 非常小（接近零）应采用替代策略（注释掉的部分）
    % if tg<=tao
    %    upx=afx-12*rx/(tao^2)-6*(vfx+vx)/tao-0;
    %    upy=afy-12*ry/(tao^2)-6*(vfy+vy)/tao-g;
    %    upz=afz-12*rz/(tao^2)-6*(vfz+vz)/tao-0;
    % end

    ux=upx+urx;
    uy=upy+ury;
    uz=upz+urz;

    ax=0+ux;
    ay=g+uy;
    az=0+uz;

    Tx=upx*m;
    Ty=upy*m;
    Tz=upz*m;

    % 制导角和油门开度（输出/监测变量）
    phic=atand(uy/ux);
    psic=-asind(uz/sqrt(ux^2+uy^2+uz^2));
    Kc=m*sqrt(ux^2+uy^2+uz^2)/Pm;

    % 动力学模型（返回微分方程）
    dx(1)=x(4);
    dx(2)=x(5);
    dx(3)=x(6);
    dx(4)=ax;
    dx(5)=ay;
    dx(6)=az;
    % 下面这一句有问题：按动力学 m_dot = -||T||/(Isp*g0) ，应为负号
    dx(13)=sqrt(Tx^2+Ty^2+Tz^2)/(Isp*g);
end

%% 求解停止条件
function [value,isterminal,direction] = odeEventFun(t,x)
    value = x(2)-0.0001; % 高度阈值 H1
    isterminal = 1;      % 触发后终止积分
    direction = -1;      % 仅在 value 从正变负时触发
end