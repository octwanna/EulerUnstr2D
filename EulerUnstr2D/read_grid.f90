 !读入网格
subroutine read_grid()
    use global
    implicit none
    integer:: i,j,a,b
		!nnode     :   node numbers
		!ncell     :   cell numbers
		!nedge     :   edge numbers
		!iedge      :   matrix for edge
		                !iedge(1,i) : edge's start point 'a'
		                !iedge(2,i) : edge's end point 'b'
		                !iedge(3,i) : edge's left cell  (all positive)
		                !iedge(4,i) : positive for right cell
		                            !iedge(4,i) -1 for wall boundary
		                            !iedge(4,i) -2 for farfiled boundary
		!xy         :   cartesian coordinates  xy(1,i) for x  xy(2,i) for y
		!icell      :   triangle cell's three nodes index
		!vol        :   cell's volume(area in 2d)
    open(10,file='10naca0012.grd')
        read(10,*)  nnode,nedge,ncell
        !allocate memory
        allocate(iedge(4,nedge))
        allocate(xy(2,nnode))
        allocate(icell(3,ncell))
        allocate(vol(ncell))
        allocate (near(2,nedge))
        allocate(vect(2,nedge))
        do i=1,nnode
            read(10,*)   xy(1,i),xy(2,i)
        end do
        do i=1,nedge
            read(10,*) iedge(1,i),iedge(2,i),iedge(3,i),iedge(4,i)
        end do
        do i=1,ncell
            read(10,*) icell(1,i),icell(2,i),icell(3,i)
        end do
        do i=1,ncell
            read(10,*) vol(i)
        end do
        do i=1,nedge
            do j=1,2
                near(j,i)=iedge(j+2,i)
            end do
        end do
       do i=1,nedge
           a=iedge(1,i)
           b=iedge(2,i)
           vect(1,i)=xy(2,b)-xy(2,a)!外法矢量dy=By-Ay
           vect(2,i)=xy(1,a)-xy(1,b)!-dx=Ax-Bx
      end do
    close(10)	
  return
end subroutine