program main
    use global
    implicit none
    real(8)::start_time,finish_time,running_time
    call read_grid()         !读入网格
    call initialize()       !流场初始化
    call cpu_time(start_time)!记录求解开始前的时间
    call solver()           !求解
    call cpu_time(finish_time)!记录求解完成后的时间
    call output()           !输出结果
    open(60,file="60cl_cd.dat")
    write(60,*)'ma    alpha(deg)  cl  cd'
    call plot()             !将读入的网格信息转化为可识别的网格信息，用以显示网格，绘制相关曲线
    close(60)
    running_time=finish_time-start_time
    write(*,*)'*******************************'
    write(*,*) "Running time:"
    write(*,*) running_time
    write(*,*)'*******************************'
    pause
end

