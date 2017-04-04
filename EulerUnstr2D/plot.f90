subroutine plot()
    use global
    implicit none
    integer::i,j
    integer::a,b,left,right,m
    real(8)::c,d,u,v,e,x,cp,p,dx,dy,cl1,cd1,cl,cd
    real(8)::wma(5,nnode),wp(4,ncell),sum_vol(nnode)
    real(8),parameter::b1=1.0
    wma=0
    sum_vol=0
    do i=1,ncell
       d=w_new(1,i)
       u=w_new(2,i)/d
       v=w_new(3,i)/d
       e=w_new(4,i)
       wp(1,i)=d
       wp(2,i)=u
       wp(3,i)=v
       wp(4,i)=(gamma-1)*(e-0.5*d*(u**2+v**2))
    end do

    do i=1,nedge
       a=iedge(1,i)
       b=iedge(2,i)
       left=iedge(3,i)
       right=iedge(4,i)
       if(right<0) then
         wma(1,a)=wma(1,a)+vol(left)*wp(1,left)
	     wma(2,a)=wma(2,a)+vol(left)*wp(2,left)
	     wma(3,a)=wma(3,a)+vol(left)*wp(3,left)
	     wma(4,a)=wma(4,a)+vol(left)*wp(4,left)
	     sum_vol(a)=sum_vol(a)+vol(left)

	     wma(1,b)=wma(1,b)+vol(left)*wp(1,left)
	     wma(2,b)=wma(2,b)+vol(left)*wp(2,left)
	     wma(3,b)=wma(3,b)+vol(left)*wp(3,left)
	     wma(4,b)=wma(4,b)+vol(left)*wp(4,left)
	     sum_vol(b)=sum_vol(b)+vol(left)
       else
         wma(1,a)=wma(1,a)+vol(left)*wp(1,left)+vol(right)*wp(1,right)
	     wma(2,a)=wma(2,a)+vol(left)*wp(2,left)+vol(right)*wp(2,right)
	     wma(3,a)=wma(3,a)+vol(left)*wp(3,left)+vol(right)*wp(3,right)
	     wma(4,a)=wma(4,a)+vol(left)*wp(4,left)+vol(right)*wp(4,right)
	     sum_vol(a)=sum_vol(a)+vol(left)+vol(right)

	     wma(1,b)=wma(1,b)+vol(left)*wp(1,left)+vol(right)*wp(1,right)
	     wma(2,b)=wma(2,b)+vol(left)*wp(2,left)+vol(right)*wp(2,right)
	     wma(3,b)=wma(3,b)+vol(left)*wp(3,left)+vol(right)*wp(3,right)
	     wma(4,b)=wma(4,b)+vol(left)*wp(4,left)+vol(right)*wp(4,right)
	     sum_vol(b)=sum_vol(b)+vol(left)+vol(right)
       end if
    end do
    !针对tecplot输出网格和变量信息文件，然后再在其中显示出图像
    open(31,file="31plot0012_pressure_ma.dat")
    open(21,file="21plot0012_rho_velocity.dat")
    !
    write(21,*)'TITLE = "plot0012_rho_velocity"'
    write(21,*)'VARIABLES="X","Y","rho","U","V"'
    write(21,*)'ZONE N=5306 ,E=10382 , F=FEPOINT, ET=TRIANGLE'
    !
    write(31,*)'TITLE = "plot0012_p_ma"'
    write(31,*)'VARIABLES="X","Y","P","M"'
    write(31,*)'ZONE N=5306 ,E=10382 , F=FEPOINT, ET=TRIANGLE'

    do i=1,nnode
       wma(1,i)=wma(1,i)/sum_vol(i)
       wma(2,i)=wma(2,i)/sum_vol(i)
       wma(3,i)=wma(3,i)/sum_vol(i)
       wma(4,i)=wma(4,i)/sum_vol(i)
       c=sqrt(gamma*wma(4,i)/wma(1,i))
       wma(5,i)=sqrt((wma(2,i))**2+(wma(3,i))**2)/c
       write(21,*)xy(1,i),xy(2,i),wma(1,i),wma(2,i),wma(3,i)
       write(31,*)xy(1,i),xy(2,i),wma(4,i),wma(5,i)
    end do
    do i=1,ncell
       write(21,*)icell(1,i),icell(2,i),icell(3,i)
       write(31,*)icell(1,i),icell(2,i),icell(3,i)
    end do
    close(21)
    close(31)

    !求翼型上下表面压力系数
    open(40,file="40cp_low.dat")
    open(50,file="50cp_up.dat")
    do i=1,nedge
      a=iedge(1,i)
      b=iedge(2,i)
      left=iedge(3,i)
      right=iedge(4,i)
      if(right==-1) then    !判断壁面边界条件（即翼型边界）
        x=0.5*(xy(1,a)+xy(1,b))
	    p=wp(4,left)
	    cp=(p-p_init)/(0.5*dens_init)
        if(xy(2,a)<0) then      !判断翼型下边界
	        write(40,"(2F10.6)") x,cp
	    else
	        write(50,"(2F10.6)") x,cp
	    end if
      end if
    end do
    close(40)
    close(50)

    !求翼型总的升力和阻力
    !open(60,file="60cl_cd.dat")
    cl1=0
    cd1=0
    do i=1,nedge
        a=iedge(1,i)
        b=iedge(2,i)
        left=iedge(3,i)
        right=iedge(4,i)
        p=wp(4,left)
        if(right==-1) then
            dx=vect(1,i)
            dy=vect(2,i)
            cl1=cl1+p*dy
            cd1=cd1+p*dx
      end if
    end do
    cl=cl1*cosd(alpha)/(0.5*dens_init*b1)-cd1*sind(alpha)/(0.5*dens_init*b1)
    cd=cl1*sind(alpha)/(0.5*dens_init*b1)+cd1*cosd(alpha)/(0.5*dens_init*b1)
    !write(60,*)'ma    alpha(deg)  cl  cd'
    write(60,*) ma,alpha,cl,cd
    write(*,*)'ma    alpha(deg)  cl  cd'
    write(*,"(1X,F12.8)") ma,alpha,cl,cd
    close(60)
    return
end subroutine plot