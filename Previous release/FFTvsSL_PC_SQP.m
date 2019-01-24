% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFTvsSL子函数
% 利用帕塞瓦尔定理约束最小二乘优化结果
% SQP法（牛顿-拉格朗日法改进版）
% 全局二次收敛
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = phi_mat;
y = h_ls;
AtA = A'*A;
Aty = A'*y;
norm_y = norm(y)*sqrt(2/length(y)) * param_y_div_x;

k = 1;
x = LS_temp / norm(LS_temp)*norm_y;
miu = 0.5;
x_len = length(x);

round_max = 100;
z = zeros(x_len+length(miu), round_max);
p = z;
rho = 0.5;
gamma = 0.2;
tao = 0.1;
mk = 0;
tic

while k<round_max
    z(:,k) = [x; miu];
    Hesse = [[2*AtA-2*miu*eye(x_len)+2/tao*(x*x'), -2*x];[2*x', 0]];
    dL = [2*AtA*x-2*Aty-2*miu*x;norm(x)-norm_y];
    dL_x = [2*AtA*x-2*Aty;norm(x)-norm_y];
    p(:,k) = Hesse\(-dL);
    
    x_add = p(1:x_len,k);
    miu_add = p(x_len+1:end,k) - miu -1/2/tao*(2*x')*x_add;
    
    
    % Amijio check
    m = floor(mk/100);
    mk = m;
    while m<100000
        if norm([2*AtA*(x+x_add*rho.^m)-2*Aty-2*(miu+miu_add*rho.^m)*(x+x_add*rho.^m);...
                -norm(x+x_add*rho.^m)+norm_y]) < sqrt(1-gamma*rho^m)*norm(dL)
            mk = m;
            break
        end
        if m<50
            m = m+1;
        elseif m<400
            m = m+10;
        else
            m = m+100;
        end
        mk = m;
    end
    if mk >= 100000
        mk
        break
    end
    x = x+rho^mk * x_add;
    miu = miu+rho^mk *miu_add;
    k = k+1;
end
toc


err2_min_SQP = norm(A*x-y)
norm_LS = norm(x)
LS_temp_new = x;





