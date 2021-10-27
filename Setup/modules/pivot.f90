module pivot
  use variables
  implicit none
  real(dp) :: omega,beta
  integer  :: resrot,dihedral,dy,dz,w
contains
  !     ##################### READ LIST ######################
  Subroutine readlist(m,config)
    implicit none
    integer             :: ierror=999,i
    integer,intent(out) :: m,config
    real(dp)            :: annel

    nlines_pivot=0
!    write(10,*) 'Reading pivotables file'
    
    !   open(21,file='Energies.dat', status='unknown')
    open(13,file='pivot_list.dat', status='old')
    open(23,file='clean_list.dat', status='unknown')
    
    read(13,*) turn, m, config, seed
    read(13,*) temp_min, temp_max, minus
    read(13,*) eta_lj, eta_chi, eta_rg, eta_bond

    annel=(temp_max - temp_min)/minus
    
   
!    print*, n_passos,turn,bin,annel
    
    allocate(list(1))
    
    do while (.true.)
       read(13, "(A1,I8)", iostat=ierror) list
       if (ierror .eq. -1) exit
! '       print*, list(1)%ss_type
       
       select case (list(1)%ss_type)
       case ('C')
         write(23, "(A1,1X,I6)") list(1)
          nlines_pivot=nlines_pivot+1
       case ('S','H','T')
          cycle
       case default
          write(10,"(A)") "Something is wrong! Check the list. I've got ", list(1)%ss_type
          stop "Error detected. Check the pivot list."
       end select
       
    end do
    bin=bin*nlines_pivot
    n_passos=((2*turn) + 1)*(bin*annel)
    close(13)
    
!    write(10,*) "There are ",nlines_pivot," pivotable residues in this structure."
    
    deallocate(list)
    
    allocate(list(1:nlines_pivot))
    
   rewind 23
    
    do i=1,nlines_pivot
       read(23,"(A1,1X,I6)") list(i)
    end do
    close(23)
    
  end Subroutine readlist
  
!     ####################### ROTATE #######################
  Subroutine Rotate(angrot)
    implicit none
    real(dp), intent(in) :: angrot
    integer              :: i
    
    beta=angrot*angle1
!!!!!!!!!!!!!Rotação psi
    select case (dihedral)
    case (1)
!       write(10,*) "Rotating psi angle of ", resrot," in", angrot,"degrees."
!       print*, "Rotating psi angle of ", resrot," in", angrot,"degrees."
       
       !     TRANSLAÇÃO CENTRALIZANDO O CA DO ÁTOMO ESCOLHIDO

!       print pdb_format, CA(resrot)
!       print pdb_format, atom(flag_CA(resrot))
       do i=1,nlines_pdb
!          print *, resrot
!          print "(3(F10.3))", atom(i)%x,atom(i)%y,atom(i)%z
          atom(i)%x=atom(i)%x-CA(resrot)%x
          atom(i)%y=atom(i)%y-CA(resrot)%y
          atom(i)%z=atom(i)%z-CA(resrot)%z
!         print "(3(F10.3))", atom(i)%x,atom(i)%y,atom(i)%z
       end do
!       print pdb_format, atom(flag_CA(resrot))
       
       !     PROJEÇÃO NO PLANO XZ
       
       w = 1
       
       call vetor()
       
       !     ROTAÇÃO EM TORNO DE Y
       
       do i=1,nlines_pdb
          call euler(i)
       end do
       
!       print pdb_format, atom(flag_CA(resrot))
!            write(*,*) 'euler'
       
       !     PROJEÇÃO NO PLANO XY
       
       w=2
       call vetor()          
!            write(*,*) 'vetor'
       
       !     ROTAÇÃO EM TORNO DE Z
       
       do i=1,nlines_pdb
          call euler(i)
       end do
       
!            write(*,*) 'euler'         
       
       !     ROTAÇÃO EM TORNO DE X
       
!       print pdb_format, atom(flag_CA(resrot))
       w=3
       !       print*, N(resrot)%residue
       
!       print pdb_format, atom(flag_CA(resrot))
       if ((N(resrot)%residue) .eq. 'PRO') return
!       print *, N(resrot)%residue,flag_CA(resrot)
       do i=1,flag_CA(resrot)
          call euler(i)
       end do
       !            write(*,*) 'euler'

       if ((N(resrot)%residue) .eq.'GLY') return
!       print*, N(resrot)%residue,flag_CB(resrot)
       do i=flag_CB(resrot), flag_last(resrot)
          call euler(i)
       end do
     
!            write(*,*) 'euler'      
       
!!!!!!!!!!!!!Rotação phi
       
    case (2)
!       write(10,*) "Rotating phi angle of ", resrot," in",angrot,"degrees."
!       print*, "Rotating phi angle of ", resrot," in",angrot,"degrees."
       
       !     TRANSLAÇÃO CENTRALIZANDO O N DO ÁTOMO ESCOLHIDO

       
       !print *, resrot
       !print pdb_format, N(resrot)
       do i=1,nlines_pdb
        !  print *, resrot
          !print "(3(F10.3))", atom(i)%x,atom(i)%y,atom(i)%z
          atom(i)%x=atom(i)%x-N(resrot)%x
          atom(i)%y=atom(i)%y-N(resrot)%y
          atom(i)%z=atom(i)%z-N(resrot)%z
       end do
       
       !     PROJEÇÃO NO PLANO XZ
       
       w = 1
       call vetor() 
       
       !     ROTAÇÃO EM TORNO DE Y
       
       do i=1,nlines_pdb
          call euler(i)
       end do
       
       !     PROJEÇÃO NO PLANO XY
       
       w=2
       call vetor() 
       
       !     ROTAÇÃO EM TORNO DE Z
       
       do i=1,nlines_pdb
          call euler(i)
       end do
       
       !     ROTAÇÃO EM TORNO DE X
       w = 3

!       print *, N(resrot)%residue,flag_CA(resrot)
       do i=flag_CA(resrot),nlines_pdb
          call euler(i)
       end do
       
    end select
  end Subroutine Rotate

!     ################## VECTOR ALIGNMENT ##################
  
  Subroutine vetor
    real(dp) :: eixo(3),proj(3) !Vetores posição
    real(dp) :: norm_eixo,norm_proj!,norm_vec !Norma dos vetores
    
!    print*, "dihedral", dihedral
    !     DEFINIÇÃO DO EIXO DE REFERÊNCIA (NO CASO O EIXO X)
    
    eixo(1) = 1
    eixo(2) = 0
    eixo(3) = 0
    norm_eixo=1

    select case (dihedral)
    case (2)
       !     DEFINIÇÃO DO VETOR N-CA
        proj(1)=atom(flag_CA(resrot))%x-atom(flag_N(resrot))%x
        proj(2)=atom(flag_CA(resrot))%y-atom(flag_N(resrot))%y
        proj(3)=atom(flag_CA(resrot))%z-atom(flag_N(resrot))%z
                    
    case (1)
       !     DEFINIÇÃO DO VETOR CA-C
       
       proj(1)=atom(flag_C(resrot))%x-atom(flag_CA(resrot))%x
       proj(2)=atom(flag_C(resrot))%y-atom(flag_CA(resrot))%y
       proj(3)=atom(flag_C(resrot))%z-atom(flag_CA(resrot))%z
       
    end select
    
!    norm_vec=sqrt((proj(1)*proj(1))+(proj(2)*proj(2))+(proj(3)*proj(3)))

!    print*, "norm_vec", norm_vec
!    print*, 'w', w
    select case (w)
    case (1)       
       !     PROJEÇÃO NO PLANO XZ Y=0
       
       proj(2)=0
       norm_proj=sqrt((proj(1)*proj(1))+(proj(2)*proj(2))+ &
            &(proj(3)*proj(3)))
       
       !     CÁLCULO DE DELTA Z PARA DETERMINAR ORIENTAÇÃO DE ROTAÇÃO

       dz=-1
       if (proj(3) .ge. 0) dz=1
       
    case (2)
       
       proj(3)=0          
       norm_proj=sqrt((proj(1)*proj(1))+(proj(2)*proj(2))+ &
            &(proj(3)*proj(3)))
       
       !     CÁLCULO DE DELTA Y PARA DETERMINAR ORIENTAÇÃO DE ROTAÇÃO

       dy=1
       if (proj(2) .ge. 0) dy=-1
       
    end select
    
    !     DETEMINAÇÃO DO ÂNGULO DE ROTAÇÃO
    omega =dacos(((proj(1)*eixo(1))+(proj(2)*eixo(2))+ &
         &(proj(3)*eixo(3)))/(norm_proj*norm_eixo))
    
!    print*, "omega", omega
    
  end Subroutine vetor
  
!     #################### EULER MATRIX ####################
  Subroutine euler(i)
    implicit none
    integer, intent(in) :: i
    real(dp)            :: b(3,3),a(3),coord(3)   !Matrizes
    
    b(:,:)=0.d0
    coord(:)=0.d0
    a(1)=atom(i)%x
    a(2)=atom(i)%y
    a(3)=atom(i)%z
    
    select case (w)
       !     Rotação em torno do eixo y
    case (1)
       
       b(1,1)=cos(dz*omega)
       b(1,3)=sin(dz*omega)
       b(2,2)=1
       b(3,1)=-sin(dz*omega)
       b(3,3)=cos(dz*omega)
       
       !     Rotação em torno do eixo z
    case (2)
       
       b(1,1)=cos(dy*omega)
       b(1,2)=-sin(dy*omega)
       b(2,1)=sin(dy*omega)
       b(2,2)=cos(dy*omega)
       b(3,3)=1
       
       !     Rotação em torno do eixo X
    case(3)
       
       b(1,1)=1
       b(2,2)=cos(beta)
       b(2,3)=sin(beta)
       b(3,2)=-sin(beta)
       b(3,3)=cos(beta)
       
    end select
    
!    print "(3F8.3)", a(:)
    !    print "(3F8.3)", b(:,:)

    coord=matmul(b,a)

    atom(i)%x=coord(1)
    atom(i)%y=coord(2)
    atom(i)%z=coord(3)
    
  end subroutine euler
     
end module pivot
  
