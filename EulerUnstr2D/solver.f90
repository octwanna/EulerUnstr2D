subroutine solver()
    use global
    implicit none
    integer::i,j,iter,r_k
    real(8)::step(5)
    real(8)::total_drho,average_drho
    step(1)=0.25
    step(2)=1.0/6.0
    step(3)=3.0/8.0
    step(4)=1./2.0
    step(5)=1.0
    do iter=1,iteration_max
      !r_k五步迭代
       do r_k=1,5
	       call flux
           if(r_k==1) then
                call timestep
                call dissipation
           end if
           call residual
           do j=1,ncell
		       do i=1,4
                   w_new(i,j)=w_orig(i,j)+step(r_k)*loc_time(j)*res(i,j)
               end do
           end do
        end do   
       !计算密度残值大小
        total_drho=0.0
	    do i=1,ncell
	       total_drho=total_drho+abs(w_new(1,i)-w_orig(1,i))
	    end do
        average_drho=1./ncell
	    average_drho=average_drho*total_drho
	    dens_res(iter)=average_drho
	    if(mod(iter,100)==0) write(*,*) iter,dens_res(iter)
        !r-k迭代作完更新变量
        w_orig=w_new
        !设定最小收敛残值，当达到此值时直接跳出迭代循环
	    if(dens_res(iter)<res_global) then
	        write(*,*) iter,dens_res(iter)
	        iterationreal=iter-1
	        exit
	    end if
    end do
end subroutine
subroutine initialize()
    use global
    implicit none
    integer::i
    !对global里面的变量进行赋值，用于整个程序过程
    open(12,file='12input.dat')
        read(12,*)
        read(12,*) ma,alpha,iteration_max
    close(12)
    write(*,*)'*********************************************'
    write(*,*)'The input mach number and incidence angle is:'
    write(*,*)'ma=',ma,'alpha=',alpha
    write(*,*)'*********************************************'
    
    cfl=3.0 
    k2=0.8   !两个根据经验选择的常数,根据文献给出的范围得到
    k4=1./70 !
    !********************************************************************
    !对整个流场赋初值
    dens_init=1.0
    velo_init=1.0
    u_init=velo_init*cosd(alpha)!注意是角度制
    v_init=velo_init*sind(alpha)
    p_init=dens_init/(gamma*ma*ma)
    e_init=p_init/((gamma-1.)*dens_init)+0.5*velo_init**2
    s_init=p_init/(dens_init**gamma)
    !*********************************************************************
    allocate(w_orig(4,ncell),w_new(4,ncell))
    allocate(inv_flux(4,ncell),vis_art(4,ncell),res(4,ncell))
    allocate(pre(ncell),loc_time(ncell))
    allocate(alpha_factor(nedge))
    allocate(dens_res(iteration_max))

    do i=1,ncell
       w_new(1,i)=dens_init
       w_new(2,i)=dens_init*u_init
       w_new(3,i)=dens_init*v_init
       w_new(4,i)=dens_init*e_init
    end do
    
    w_orig=w_new
    write(*,*) "The initialization of flowfield has been finished!"
end subroutine
subroutine dissipation()!计算人工耗散项
    use global
    implicit none
    integer::i,j,k,left,right
    real(8)::w(4),rho,uu,vv,pp
    real(8)::e2,e4,pl,pr,d2,d4
    real(8)::vk!振动传感器
    real(8)::ww(4,ncell)

    do j=1,ncell
       do i=1,4
          vis_art(i,j)=0.0
	      ww(i,j)=0.0
       end do
    end do

    !更新守恒变量，求解每个单元的压强
    do i=1,ncell
       w(1)=w_new(1,i)
       w(2)=w_new(2,i)
       w(3)=w_new(3,i)
       w(4)=w_new(4,i)
       call pressure(w,rho,uu,vv,pp)
       pre(i)=pp
    end do
    !求解ww
    do i=1,nedge
       left=near(1,i)
       right=near(2,i)
       if(right>0) then
           do k=1,4
	          ww(k,left)=ww(k,left)+0.5*(w_new(k,right)-w_new(k,left))
		      ww(k,right)=ww(k,right)+0.5*(w_new(k,left)-w_new(k,right))
	       end do
       end if
    end do

    !按边循环求解每个单元的人工粘性
    do i=1,nedge
       left=near(1,i)
       right=near(2,i)
       if(right>0) then
          pr=pre(right)
	      pl=pre(left)
	      vk=abs(pr-pl)/abs(pr+pl)
	      e2=k2*vk
	      e4=max(0.0,k4-e2) 
      
	      do k=1,4
	         d2=alpha_factor(i)*e2*(w_new(k,right)-w_new(k,left))
		     d4=-1.0*alpha_factor(i)*e4*(ww(k,right)-ww(k,left)) 
		     vis_art(k,left)=vis_art(k,left)+d2+d4
		     vis_art(k,right)=vis_art(k,right)-d2-d4
	      end do
       end if
    end do
end subroutine
subroutine timestep !求解当地时间步长
    use global
    implicit none
    integer::i,left,right !左右单元编号
    !赋初值
    do i=1,ncell
       loc_time(i)=0.0
    end do
    !按边循环
    do i=1,nedge
       left=near(1,i)
       right=near(2,i)
       loc_time(left)=loc_time(left)+alpha_factor(i)
       if(right>0) then !针对右单元为内部单元的时候
          loc_time(right)=loc_time(right)+alpha_factor(i)
       end if 
    end do
    do i=1,ncell
       loc_time(i)=vol(i)*cfl/loc_time(i)
    end do
end subroutine
subroutine flux()!按边循环求解通量
    use global
    implicit none
    integer::left,right !左右单元
    integer::i,jalpha_factor
    real(8)::deltay,deltax,delta,vn
    real(8)::w(4) !用来存放边的四个通量
    real(8)::rho,uu,vv,pp
    real(8)::aa,aa2  !代表音速与音速平方
    real(8)::f1,f2,f3,f4 !四个无粘通量

    inv_flux=0.0
    alpha_factor=0.0

    !按边循环
    do i=1,nedge
       left=near(1,i)
       right=near(2,i)
       deltay=vect(1,i)
       deltax=vect(2,i)
       delta=deltax**2+deltay**2
  
       if(right>0) then
            w(1)=0.5*(w_new(1,left)+w_new(1,right))
            w(2)=0.5*(w_new(2,left)+w_new(2,right))
            w(3)=0.5*(w_new(3,left)+w_new(3,right))
            w(4)=0.5*(w_new(4,left)+w_new(4,right))
            call pressure(w,rho,uu,vv,pp)
		    vn=uu*deltay+vv*deltax
            aa2=gamma*pp/rho
            aa=sqrt(aa2)
            alpha_factor(i)=abs(vn)+aa*sqrt(delta)
        
        else if(right==-1) then
            w(1)=w_new(1,left)
            w(2)=w_new(2,left)
            w(3)=w_new(3,left)
            w(4)=w_new(4,left)
        
		    call pressure(w,rho,uu,vv,pp) 
            vn=0!壁面边界条件，法向速度为0
            aa2=gamma*pp/rho
            aa=sqrt(aa2)
            alpha_factor(i)=abs(vn)+aa*sqrt(delta)
        else if(right==-2) then
		    w(1)=w_new(1,left)
            w(2)=w_new(2,left)
            w(3)=w_new(3,left)
            w(4)=w_new(4,left)
            call far_field(deltay,deltax,w,rho,uu,vv,pp)
		    vn=uu*deltay+vv*deltax   !远场边界条件
            aa2=gamma*pp/rho
            aa=sqrt(aa2)
            alpha_factor(i)=abs(vn)+aa*sqrt(delta)
	    end if

            f1=vn*rho
            f2=vn*rho*uu+pp*deltay
            f3=vn*rho*vv+pp*deltax
            f4=vn*(w(4)+pp)
   
            inv_flux(1,left)=inv_flux(1,left)+f1
            inv_flux(2,left)=inv_flux(2,left)+f2
            inv_flux(3,left)=inv_flux(3,left)+f3
            inv_flux(4,left)=inv_flux(4,left)+f4
   
       if(right>0) then
            inv_flux(1,right)=inv_flux(1,right)-f1
            inv_flux(2,right)=inv_flux(2,right)-f2
            inv_flux(3,right)=inv_flux(3,right)-f3
            inv_flux(4,right)=inv_flux(4,right)-f4
       end if
    end do
end subroutine
subroutine far_field(deltay,deltax,w,rho,uu,vv,pp)
    use global
    implicit none
    integer::i
    real(8)::deltay,deltax,rho,uu,vv,pp,w(4)  !w用来存放外边界左单元的信息
    real(8)::ww(4)
    real(8)::dy,dx !dy,dx单位化法矢量
    real(8)::vn,vn_init,vne,ue,ve !代表法向速度,e代表单元速度
    real(8)::aa,aa2     !音速和音速的平方
    real(8)::r_plus,r_minus,s,af !s熵,af代表边界的音速
    dy=deltay/sqrt(deltay**2+deltax**2)
    dx=deltax/sqrt(deltay**2+deltax**2)
    call pressure(w,rho,uu,vv,pp)  !用远场左单元的信息，来计算
       ue=uu
       ve=vv
       vne=ue*dy+ve*dx          !左单元的法向速度
       vn_init=u_init*dy+v_init*dx !无穷远处法向速度
       aa2=gamma*pp/rho     
       aa=sqrt(aa2)       !左单元音速
       r_plus=vne+2.*aa/(gamma-1.)                  !r_plus
       r_minus=vn_init-2./ma/(gamma-1.)               !rm initinite
       vn=0.5*(r_plus+r_minus)                         !远场的法向速度和音速(难点)
       af=(gamma-1.)*0.25*(r_plus-r_minus)

    if(abs(vn_init/(1./ma))<1.) then

         if(vn<0.0) then          !亚音速入流
	        s=s_init
		    uu=u_init+(vn-vn_init)*dy
		    vv=v_init+(vn-vn_init)*dx
	     else                    !亚音速出流
	        s=pp/rho**gamma
		    uu=ue+(vn-vne)*dy
		    vv=ve+(vn-vne)*dx
	     end if

    !亚音速
	     rho=(af**2/(gamma*s))**(1./(gamma-1.))
	     pp=rho*af*af/gamma
	     ww(1)=rho
	     ww(2)=rho*uu
	     ww(3)=rho*vv
	     ww(4)=pp/(gamma-1.)+0.5*rho*(uu**2+vv**2)      
    else 
    !超音速
        if(vn<0.) then  !超音速入流，取来流值
	       ww(1)=dens_init
	       ww(2)=dens_init*u_init
	       ww(3)=dens_init*v_init
	       ww(4)=dens_init*e_init
        else  !超音速出流，取内场值
            ww=w
        end if
    end if
    call pressure(ww,rho,uu,vv,pp)
    return
end subroutine
subroutine pressure(w,rho,u,v,p)
    use global
    implicit none
    real(8)::w(4),rho,u,v,p
    rho=w(1)
    u=w(2)/w(1)
    v=w(3)/w(1)
    p=(gamma-1.0)*(w(4)-0.5*rho*(u**2+v**2))!根据理想气体总能公式求解单元压强

    return
end subroutine

subroutine residual()
    use global
    implicit none
    integer::i,j
    res=0.0		!无意义
    do j=1,ncell
       do i=1,4
	      res(i,j)=(vis_art(i,j)-inv_flux(i,j))/vol(j)
       end do
    end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output()
    use global
    implicit none
    integer::i
    open(20,file="20output_w.dat")
    write(20,*) "alpha=",alpha,"ma=",ma
    do i=1,ncell
       write(20,"(1X,F16.8,\)") w_new(1,i),w_new(2,i),w_new(3,i),w_new(4,i)
       write(20,*)' '
    end do
    open(30,file="30residual_history.dat")
    do i=1,iterationreal,2
       write(30,*) i,dens_res(i)
    end do
    close(20)
    close(30)
end subroutine
