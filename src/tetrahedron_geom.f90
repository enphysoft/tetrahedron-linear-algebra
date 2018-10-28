module tetrahedron_geom
      use parmkind
      implicit none
      type tetrahedron 
              ! values to be inputted
              real(rkind)  :: ndval(4)   = 0.0D0 
              real(rkind)  :: postn(4,3) = 0.0D0 

              real(rkind)  :: ttvol      = 0.0D0 
              real(rkind)  :: detJ0      = 0.0D0 
              real(rkind)  :: Jmat0(4,4) = 0.0D0
              real(rkind)  :: Jinv0(4,4) = 0.0D0 

              ! values to be set
              real(rkind)  :: locat(3)   = 0.0D0 
              ! values to be calculated
              real(rkind)  :: ttloc(4)   = 0.0D0 ! ttloc(1) = 1, ttloc(i+1)=locat(i) 
              real(rkind)  :: ttcrd(4)   = 0.0D0 
              real(rkind)  :: apprx      = 0.0D0 
              real(rkind)  :: pderv(3)   = 0.0D0 
              real(rkind)  :: ptvol(4)   = 0.0D0 
              real(rkind)  :: Jmat(4,4,4)= 0.0D0
              real(rkind)  :: detJ(4)    = 0.0D0 
              logical      :: inside
      end type tetrahedron
      public :: make_Jmatdet, interpol_ttrhdrn 

   contains
      
      subroutine calc_ttrhdrn (TTRcalc)
              implicit none 
              type(tetrahedron) :: TTRcalc
              call check_inside_ttrhdrn (TTRcalc)
              if(TTRcalc%inside.eqv..true.) then
                  call interpol_ttrhdrn     (TTRcalc)
                  call gradient_ttrhdrn     (TTRcalc)
              else
              endif
      end subroutine calc_ttrhdrn 

      subroutine setup_ttrhdrn (TTRsetup)
              implicit none
              type(tetrahedron) :: TTRsetup
              integer (ikind) :: ierr 
              TTRsetup%Jmat0(1,:) = 1.0d0
              do i = 2,4; do j = 1,4
                   TTRsetup%Jmat0(i,j) = TTRsetup%postn(j,i-1) 
              enddo ; enddo
              call make_Jmatdet (TTRsetup%Jmat0, TTRsetup%detJ0)
              TTRsetup%ttvol = abs( TTRsetup%detJ0 ) / 6.0d0 
              call make_Jmatinv (TTRsetup%Jmat0, TTRsetup%Jinv0, ierr)
              if (ierr < 0 ) write(*,*) "Zero determiant in make_Jmatinv "
      end subroutine setup_ttrhdrn

      subroutine make_Jmat (Jmat,postn)
      implicit none
      real(rkind) , intent(inout) :: Jmat(4,4)
      real(rkind) , intent(in)    :: postn(4,3)
      Jmat(1,1:4) = 1.0d0
      Jmat(2,1)   = postn(1,1)
      Jmat(3,1)   = postn(1,2)
      Jmat(4,1)   = postn(1,3)
      do i = 1,3
        do j = 1,4
          Jmat(i+1,j) = postn(j,i)
        enddo
      enddo
      end subroutine make_Jmat

      subroutine make_Jmatdet (Jmat, det)
              implicit none
              real(rkind)  , intent(in)  :: Jmat(4,4)
              real(rkind)   :: det
              real(rkind)   :: A(4,4)
              A = Jmat
              det=&
            A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))   &
                   +A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))   &
                   +A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)) ) &
           -A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))   &
                   +A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))   &
                   +A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)) ) &
           +A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))   &
                   +A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))   &
                   +A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)) ) &
           -A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))   &
                   +A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))   &
                   +A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)) )
      end subroutine make_Jmatdet

      subroutine make_Jmatinv (Jmat, Jinv, ierr)
              implicit none
              real(rkind)  , intent(in)  :: Jmat(4,4)
              real(rkind)  , intent(out) :: Jinv(4,4)
              integer (ikind), intent(out) :: ierr
              real(rkind) :: det, A(4,4), cofactor(4,4)
      
              A = Jmat
              call make_Jmatdet (A, det)

              if(abs(det) > 0.0d0 ) then
              COFACTOR(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
                      +A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))&
                      +A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
              COFACTOR(1,2) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
                      +A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))&
                      +A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
              COFACTOR(1,3) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))&
                      +A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))&
                      +A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
              COFACTOR(1,4) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))&
                      +A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))&
                      +A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
              COFACTOR(2,1) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
                      +A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))&
                      +A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
              COFACTOR(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
                      +A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))&
                      +A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
              COFACTOR(2,3) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))&
                      +A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))&
                      +A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
              COFACTOR(2,4) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))&
                      +A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))&
                      +A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
              COFACTOR(3,1) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))&
                      +A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))&
                      +A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
              COFACTOR(3,2) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))&
                      +A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))&
                      +A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
              COFACTOR(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))&
                      +A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))&
                      +A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
              COFACTOR(3,4) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))&
                      +A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))&
                      +A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
              COFACTOR(4,1) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))&
                      +A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))&
                      +A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
              COFACTOR(4,2) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))&
                      +A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))&
                      +A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
              COFACTOR(4,3) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))&
                      +A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))&
                      +A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
              COFACTOR(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
                      +A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))&
                      +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
              Jinv = transpose(cofactor) / det
              ierr = 0
              else
              Jinv= 0.0d0
              ierr = -1
              endif        
      end subroutine make_Jmatinv

      subroutine interpol_ttrhdrn (TTRitpl)
              implicit none
              type (tetrahedron), intent(inout) :: TTRitpl
                    TTRitpl%ttloc(1) = 1.0d0
                    TTRitpl%ttloc(2) = TTRitpl%locat(1)
                    TTRitpl%ttloc(3) = TTRitpl%locat(2)
                    TTRitpl%ttloc(4) = TTRitpl%locat(3)
                    TTRitpl%ttcrd    = matmul(TTRitpl%Jinv0,(TTRitpl%ttloc))
                    TTRitpl%apprx    = dot_product (TTRitpl%ndval,TTRitpl%ttcrd)
      end subroutine interpol_ttrhdrn

      subroutine gradient_ttrhdrn (TTRgrad) 
              implicit none
              type (tetrahedron), intent(inout) :: TTRgrad
                    TTRgrad%pderv(1) = dot_product (TTRgrad%ndval, TTRgrad%Jinv0(:,2))
                    TTRgrad%pderv(2) = dot_product (TTRgrad%ndval, TTRgrad%Jinv0(:,3))
                    TTRgrad%pderv(3) = dot_product (TTRgrad%ndval, TTRgrad%Jinv0(:,4))
      end subroutine gradient_ttrhdrn 

      subroutine pvolume_ttrhdrn (TTRpvol) 
              implicit none
              type (tetrahedron), intent(inout) :: TTRpvol
              real(rkind) :: Jmat(4,4,4), detval(4)

              TTRpvol%Jmat(1,:,:) =  TTRpvol%Jmat0
              TTRpvol%Jmat(2,:,:) =  TTRpvol%Jmat0
              TTRpvol%Jmat(3,:,:) =  TTRpvol%Jmat0
              TTRpvol%Jmat(4,:,:) =  TTRpvol%Jmat0

              TTRpvol%Jmat(1,2:4,1) = TTRpvol%locat(1:3)
              TTRpvol%Jmat(2,2:4,2) = TTRpvol%locat(1:3)
              TTRpvol%Jmat(3,2:4,3) = TTRpvol%locat(1:3)
              TTRpvol%Jmat(4,2:4,4) = TTRpvol%locat(1:3)

              Jmat = TTRpvol%Jmat
              do i=1,4
                  call make_Jmatdet (Jmat(i,:,:), detval(i))
                  TTRpvol%detJ(i) = detval(i)
              enddo
      end subroutine pvolume_ttrhdrn

      subroutine check_inside_ttrhdrn (TTRisdchk)
              implicit none
              type (tetrahedron), intent(inout) :: TTRisdchk
              real (rkind) :: ttvol, ptvol(4)
              call pvolume_ttrhdrn (TTRisdchk) 
              TTRisdchk%ttvol = TTRisdchk%detJ0/6.0d0
              TTRisdchk%ptvol = TTRisdchk%detJ/6.0d0
              ttvol = TTRisdchk%ttvol
              ptvol = TTRisdchk%ptvol

              !write(*,*) ttvol
              !write(*,*) ptvol 
              if (ttvol*ptvol(1) > 0.0d0 .and. &
                  ttvol*ptvol(2) > 0.0d0 .and. &
                  ttvol*ptvol(3) > 0.0d0 .and. &
                  ttvol*ptvol(4) > 0.0d0 ) then
                 TTRisdchk%inside = .true.
!                 call  interpol_ttrhdrn (TTRisdchk)
!                 call  gradient_ttrhdrn (TTRisdchk)
              else 
                 TTRisdchk%inside = .false.
              endif
      end subroutine check_inside_ttrhdrn


end module tetrahedron_geom      
