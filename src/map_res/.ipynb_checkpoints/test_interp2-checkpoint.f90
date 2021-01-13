  !____________________________________________________________________
  !subroutine mesh_interp(cx,cy,cz,pfnx,pfny,pfnz,pfxcb,pfycb,pfzcb,e4d_mfp)
  subroutine mesh_interp(pfxcb,pfycb,pfzcb,e4d_mfp,pfnx,pfny,pfnz)
    implicit none
    !integer :: cx,cy,cz,i
    integer, intent(in) :: pfnx ,pfny ,pfnz
    real, intent(in), dimension(pfnx+1) :: pfxcb
    real, intent(in), dimension(pfny+1) :: pfycb
    real, intent(in), dimension(pfnz+1) :: pfzcb
    character*80, intent(in) :: e4d_mfp

    !f2py directives
!f2py integer intent(hide), depend(pfxcb) :: pfnx=shape(pfxcb,0)-1
!f2py integer intent(hide), depend(pfycb) :: pfny=shape(pfycb,0)-1
!f2py integer intent(hide), depend(pfzcb) :: pfnz=shape(pfzcb,0)-1
	!f2py intent(in) pfnx,pfny,pfnz
!f2py intent(in) pfxcb,pfycb,pfzcb
!f2py intent(in) e4d_mfp
    
    real :: j1,j2
    integer :: nnodes,nelem
    integer, dimension(:,:), allocatable :: elements
    real, dimension(:,:), allocatable :: nodes
    real :: xorig,yorig,zorig

    logical, dimension(:), allocatable :: flags
    integer :: my_nelem,cct,indx,indy,indz,k,j,i
    integer*8 :: cnt,n_tets,cpos
    real :: xmn,ymn,zmn,xmx,ymx,zmx,imn,imx
    real, dimension(:), allocatable :: midx,midy,midz,w
    integer, dimension(:), allocatable :: rw,v
    real, dimension(9,3) :: pts
    real, dimension(4,3) :: mfaces
    real, dimension(3) :: vc
    real, dimension(8) :: wi,C
    integer, dimension(8) :: vi
    real, dimension(72) :: twi
    integer, dimension(72) :: tvi
    write(*,*) pfnx,pfny,pfnz

    write(*,*)
    write(*,*) "Interpolation matrix builder reporting for duty ..."
    write(*,*)
    write(*,*) "Number of PFLOTRAN cells in x: ",pfnx
    do i=1,pfnx+1
       write(*,*) 'x node ',i,pfxcb(i)
    end do
    
    write(*,*)
    write(*,*) "Number of PFLOTRAN cells in y: ",pfny
    do i=1,pfny+1
       write(*,*) 'y node ',i,pfycb(i)
    end do
    
    write(*,*)
    write(*,*) "Number of PFLOTRAN cells in z: ",pfnz
    do i=1,pfnz+1
       write(*,*) 'z node ',i,pfzcb(i)
    end do
    
    write(*,*)
    write(*,*) "E4D mesh file prefix: ",trim(e4d_mfp)

    !read and translate the nodes
    open(11,file=trim(e4d_mfp)//".1.node",status='old',action='read')
    read(11,*) nnodes
    allocate(nodes(nnodes,3))
    do i=1,nnodes
       read(11,*) j1,nodes(i,:)
    end do
    close(11)
    open(11,file=trim(e4d_mfp)//".trn",status='old',action='read')
    read(11,*) xorig,yorig,zorig
    close(11)
    nodes(:,1)=nodes(:,1)+xorig
    nodes(:,2)=nodes(:,2)+yorig
    nodes(:,3)=nodes(:,3)+zorig
    
    !read the elements
    open(10,file=trim(e4d_mfp)//".1.ele",status="old",action="read")
    read(10,*) nelem,j1,j2
    allocate(elements(nelem,4))
    do i=1,nelem
       read(10,*) j1,elements(i,1:4)
    end do
    close(10)
   
    !!get the midpoints of the grid (not the e4d_mesh)
    xmn=1e15
    ymn=xmn
    zmn=xmn
    xmx=-1e15
    ymx=xmx
    zmx=xmx
    allocate(midx(pfnx),midy(pfny),midz(pfnz))
    do i=1,pfnx
       midx(i) = 0.5*(pfxcb(i+1)+pfxcb(i))
       if(midx(i)<xmn) xmn=midx(i)
       if(midx(i)>xmx) xmx=midx(i)
      
    end do
    do i=1,pfny
       midy(i) = 0.5*(pfycb(i+1)+pfycb(i))
       if(midy(i)<ymn) ymn=midy(i)
       if(midy(i)>ymx) ymx=midy(i)
    
    end do
    do i=1,pfnz
       midz(i) = 0.5*(pfzcb(i+1)+pfzcb(i))
       if(midz(i)<zmn) zmn=midz(i)
       if(midz(i)>zmx) zmx=midz(i)
    end do
    my_nelem = nelem
    allocate(flags(my_nelem))
    flags = .true.
    cnt = 0
    n_tets = 0

    !!find which elements we need to map
    do i=1,nelem
       cnt=cnt+1
       imx=maxval(nodes(elements(i,1:4),1))
       imn=minval(nodes(elements(i,1:4),1))
       if((imn<xmn .or. imx>xmx) .and. (pfnx > 1)) then
          flags(cnt)=.false.
          goto 100
       end if
       imx=maxval(nodes(elements(i,1:4),2))
       imn=minval(nodes(elements(i,1:4),2))
       if((imn<ymn .or. imx>ymx) .and. (pfny > 1)) then
          flags(cnt)=.false.
          goto 100
       end if
       imx=maxval(nodes(elements(i,1:4),3))
       imn=minval(nodes(elements(i,1:4),3)) 
       if((imn<zmn .or. imx>zmx) .and. (pfnz>1)) then
          flags(cnt)=.false.
          goto 100
       end if
       n_tets=n_tets+1
100    continue      
    end do
    
    write(*,*) "Building interpolation matrix"
	
    !!allocate and start the interpolation
    allocate(rw(72*n_tets),v(72*n_tets),w(72*n_tets))
    cnt=0
    cct=0
    cpos=0
    do i=1,nelem
       cct=cct+1
       if(flags(cct)) then
          !get the nine points
          cnt=cnt+1
          do j=1,3
             pts(1,j)=.25*sum(nodes(elements(i,1:4),j)) !tet midpoint
          end do
        
          do j=1,3
             mfaces(1,j) = sum(nodes(elements(i,[1,2,3]),j))/3
             mfaces(2,j) = sum(nodes(elements(i,[1,2,4]),j))/3
             mfaces(3,j) = sum(nodes(elements(i,[1,3,4]),j))/3
             mfaces(4,j) = sum(nodes(elements(i,[2,3,4]),j))/3       
          end do
          do j=1,4
             vc=mfaces(j,:)-pts(1,:)
             pts(j+1,:) = pts(1,:) + 0.5*vc
             vc=nodes(elements(i,j),:)-pts(1,:)
             pts(j+5,:) = pts(1,:) + 0.5*vc
          end do
          

          !get the weighting function for each point
          do j=1,9
             if(pfnx .eq. 1) then
                indx = 1
             else
                do k=1,pfnx
                   if(midx(k)>pts(j,1)) then
                      indx=k-1
                      exit
                   end if
                end do
             end if

             if(pfny.eq.1) then
                indy = 1
             else
                do k=1,pfny
                   if(midy(k)>pts(j,2)) then
                      indy=k-1
                      exit
                   end if
                end do
             end if
             
             if(pfnz .eq.1) then
                indz = 1
             else
                do k=1,pfnz
                   if(midz(k)>pts(j,3)) then
                      indz=k-1
                      exit
                   end if
                end do
             end if

             
             if(pfny .eq.1) then
                !for 2D problem in x,z
                vi(5) = indx + (indy-1)*pfnx + (indz-1)*pfnx*pfny
                vi(6) = vi(5)
                vi(8) = vi(5)+1
                vi(7) = vi(8);
                vi(1) = vi(5)+pfnx
                vi(2) = vi(1);
                vi(3) = vi(1)+1;
                vi(4) = vi(3);
                
             else
                vi(5) = indx + (indy-1)*pfnx + (indz-1)*pfnx*pfny
                vi(8) = vi(5)+1
                vi(6) = vi(5)+pfnx
                vi(7) = vi(6)+1
                vi(1:4) = vi(5:8) + pfnx*pfny
             end if

              if(pfnx .eq.1) then
                 C(1)=0
              else
                 C(1)=(pts(j,1)-midx(indx))/(midx(indx+1)-midx(indx));
              end if
              C(2)=C(1);
              C(3)=C(1);
              C(4)=C(1);

              if(pfny .eq.1) then
                 C(5)=0
              else
                 C(5)=(pts(j,2)-midy(indy))/(midy(indy+1)-midy(indy));
              end if
              C(6)=C(5);
              
              if(pfnz .eq.1) then
                 C(7)=0
              else
                 C(7)=(pts(j,3)-midz(indz))/(midz(indz+1)-midz(indz));
              end if

              wi(1)=(1-C(1))*(1-C(6))*C(7);
              wi(2)=(1-C(2))*C(6)*C(7);
              wi(3)=C(2)*C(6)*C(7);
              wi(4)=C(1)*C(7)*(1-C(6));
              wi(5)=(1-C(3))*(1-C(5))*(1-C(7));
              wi(6)=C(5)*(1-C(4))*(1-C(7));
              wi(7)=C(4)*C(5)*(1-C(7));
              wi(8)=C(3)*(1-C(5))*(1-C(7));
           
              !write(*,*) vi
              !write(*,*) wi

              tvi(8*(j-1)+1:8*j)=vi
              twi(8*(j-1)+1:8*j)=wi
              
           end do
           
           !consolidate the weights
           do j=1,71
              do k=j+1,72
                 if(tvi(k)==tvi(j)) then
                    twi(j)=twi(j)+twi(k)
                    twi(k)=0
                    tvi(k)=0;
                 end if
              end do
           end do
           twi=twi/9
          
           do j=1,72
              if(tvi(j) .ne. 0) then
                 cpos=cpos+1;
                 rw(cpos)=i
                 v(cpos)=tvi(j)
                 w(cpos)=twi(j)
              end if
           end do
        end if
     
    end do
    
    !do i=1,cpos
    !   write(*,*) i,rw(i),v(i),w(i)
    !end do
    open(10,file=trim(e4d_mfp)//"_map.bin",form='unformatted')
    write(10) pfnx
    write(10) pfny
    write(10) pfnz
    write(10) nelem
    write(10) cpos
    write(10) rw(1:cpos)
    write(10) v(1:cpos)
    write(10) w(1:cpos)
    close(10)
    write(*,*) "Writing interpolation matrix"
    write(*,*) nelem,maxval(rw),maxval(v)
  end subroutine mesh_interp
  !__________________________________________________________________

 
