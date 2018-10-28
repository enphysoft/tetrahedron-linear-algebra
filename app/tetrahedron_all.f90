use parmkind
use tetrahedron_geom
      implicit none
      type(tetrahedron) :: TTR
      integer (ikind) :: ii, jj, ierr
      real (rkind)  :: Jmat(4,4),Jinv(4,4),det  


      ! setting positions of four nodes
      open(9,file='postn.dat',status='old')
      do i = 1,4; read(9,*) (TTR%postn(i,j),j=1,3) ;enddo
      write(*,*) "TTR%postn"
      do i=1,4; write(*,"(4(2X,F12.8))") TTR%postn(i,:); enddo ; write(*,*)

      ! setting node values 
      open(10,file='ndval.dat',status='old')
      do i = 1,4; read(10,*) TTR%ndval(i) ;enddo
      write(*,*) "TTR%ndval"
      write(*,"(4(2X,F12.8))") (TTR%ndval(j),j=1,4) ; write(*,*)

      ! setting up Jacobian mattrix 
      call setup_ttrhdrn (TTR) 

      write(*,*) "TTR%Jamt0"
      do i=1,4; write(*,"(4(2x,F12.8))") (TTR%Jmat0(i,j), j=1,4)
      enddo   ; write(*,*)

      write(*,*) "det         = ", TTR%detJ0
      write(*,*) "TTR%ttvol   = ", TTR%ttvol
      write(*,*) "Analytic vol= ", (2.0d0)**3 / (6.0d0*sqrt(2.0d0))
      write(*,*)  

      ! calcuating inverse matrix of J
        write(*,*) "TTR%Jinv0"
        do ii=1,4
          write(*,"(4(2x,F12.8))") (TTR%Jinv0(ii,jj), jj=1,4)
        enddo
      write(*,*)  


      ! setting interpolation coordinates
      open(11,file='locat.dat',status='old')
      do j = 1,3; read(11,*) TTR%locat(j) ;enddo
      write(*,*) "TTR%locat"
      write(*,"(4(2X,F12.8))") (TTR%locat(j),j=1,3) ; write(*,*)

      !call pvolume_ttrhdrn (TTR) 
      !call check_inside_ttrhdrn (TTR)
      call calc_ttrhdrn (TTR)
      write(*,*)     "TTR%detJ0  ="
      write(*,"(4(2X,F12.8))")     TTR%detJ0
      do i=1,4
          write(*,*) "TTR%detJ(i)=", TTR%detJ(i) 
      enddo
      write(*,*) "sum_TTR%detJ(i)=", sum(TTR%detJ) 

      

      
      
      write(*,*) TTR%inside
      write(*,"(' TTR%ttcrd =  ',4(2X,F12.8))")  TTR%ttcrd
      write(*,"(' TTR%ndval =  ',4(2X,F12.8))")  TTR%ndval
      write(*,*) 
      write(*,*) "TTR%apprx = ", TTR%apprx
      write(*,"(' TTR%pderv =  ',4(2X,F12.8))")  TTR%pderv


      stop
      end

