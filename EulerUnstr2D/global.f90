module global
    implicit none
    real(8),parameter::pi=3.1415926535
    real(8),parameter::res_global=1e-6      !最小收敛残值，达到之后就终止迭代
    real(8),parameter:: gamma=1.4
    
    integer::iteration_max,iterationreal   !iteration_max代表实际迭代步数,iterationreal代表实际迭代步数
    integer::ncell,nnode,nedge        
    integer,allocatable::near(:,:),iedge(:,:),icell(:,:)!网格信息
    real(8),allocatable::xy(:,:),vol(:) !网格信息
    real(8),allocatable::vect(:,:)             !vect存储法矢信息
    
    real(8)::dens_init,u_init,v_init,p_init,e_init,velo_init,s_init !初始化的变量
    real(8)::cfl,k2,k4
    real(8)::alpha,ma   !来流角和来流马赫数
    
    real(8),allocatable::w_orig(:,:),w_new(:,:)   ! 用于r-k 5步迭代法，更新守恒变量
    real(8),allocatable::inv_flux(:,:)
    real(8),allocatable::vis_art(:,:),res(:,:)
    real(8),allocatable::pre(:),loc_time(:) !压力和当地时间步长
    real(8),allocatable::alpha_factor(:)    !人工粘性项的比例因子
    real(8),allocatable::dens_res(:)        !每一步迭代之后的密度残值
end module
