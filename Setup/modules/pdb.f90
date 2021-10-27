Module pdb_file
  use variables
  implicit none
  real(dp), parameter :: water_box_size=119.7, thickness=3.d0, close_dist=3.5
  integer, parameter  :: n_waters=60565
  
  character(40)        :: cleanpdb,PDB_output,EXP_output
  integer              :: CA_index=0,j,nhbond
  real(dp),allocatable :: dist(:,:), radii(:),rvw(:,:)
  real(dp)             :: max_coord(3), min_coord(3),d_ete
  
contains 
  !     #################### READ PDB #######################
  Subroutine readpdb()
    implicit none
    integer :: i=1
    
    i=1
    nlines_pdb=0
!    write(10,*) "Reading PDB file"
    open(11, file=PDB_input, status='old') !Input
    
    write(cleanpdb, '( "clean_", A34 )' )  PDB_input
    
    open(12, file=cleanpdb, status='Unknown') !Input
    
    allocate(atom(1))
    
    atom(1)%dummy="start"
    i=1
    do while (atom(1)%dummy .ne. 'END')
       read(11,"(A6)") atom%dummy
       if (atom(1)%dummy .ne. "ATOM") then
          cycle
       else          
          backspace 11
          read(11,pdb_format) atom
          if (atom(1)%alternate .eq. 'B') cycle
          atom(1)%serialnumber=i
          write(12,pdb_format) atom
          i=i+1
       end if
    end do
    close(11)
    
    write(12,'(A3)') "END"
    nlines_pdb=i-1 !Defining the number of lines and allow to use i
    !again
    
    deallocate (atom)
    allocate(atom(0:nlines_pdb))
    
    rewind 12
    
    do i=1,nlines_pdb
       read(12,pdb_format) atom(i)
      
       if (atom(i)%name .eq. " CA ") then
          CA_index=CA_index+1
          res_count=res_count+1
       end if
    end do
    close (12)

   
   
    CA_index=CA_index+atom(1)%ressqnumber-1

  end Subroutine readpdb
  
!     ##################### READ CA ########################
  Subroutine readCA(CA_count,PDB_CA)
    implicit none
    integer, intent(in)       :: CA_count
    character(40), intent(in) :: PDB_CA
    integer                   :: ierr,index=0

!    write(10,*) "Reading CA file"
    allocate(CA(0:(atom(1)%ressqnumber+CA_count)))
    
    index=0
        
    open(13, file=PDB_CA, status='old')
      
    do while (CA(index)%dummy .ne. "END")
       index=index+1
       read(13,pdb_format,iostat=ierr) CA(index)
!       print *, CA(index)%dummy
       if (ierr .ne. 0)  stop "Something is wrong! Check CA file"
      
    end do
    
    close(13)
  end Subroutine readCA
  !     ####################### LABEL ########################
  Subroutine label(PDB_CA,index,CA_count)
    implicit none
    integer                   :: j=0,start=0,i
    integer,intent(out)       :: index,CA_count
    character(40),intent(out) :: PDB_CA
    logical                   :: ok
    
!    write(10,*) "Labeling the atoms"
    write(PDB_CA, '( "CA_", A37 )' )  PDB_input

    
!    open(12, file="pdb.tmp", action='write', position='append') !Input
    open(13, file=PDB_CA, status='unknown')
    
    start=atom(1)%ressqnumber
    atom(0)%serialnumber=0
    atom(0)%ressqnumber=0
    CA_count=0
    index=atom(1)%ressqnumber
    max_coord(:)=0
    min_coord(:)=100000
    
!    print *, CA_index, atom(1)%ressqnumber
!    print*, first
    select case (first)
    case (0)
       allocate(CA(start:CA_index),N(CA_index),flag_CB(CA_index))
       allocate(flag_C(CA_index),flag_last(0:(CA_index+1)),flag_CA(CA_index))
       allocate(flag_N(CA_index),radii(1:nlines_pdb),rvw(1:nlines_pdb,1:nlines_pdb))
       
!       write(12,pdb_format)atom(1:)
!       write(12,"(A3)") "TER"
       
       do i=1,nlines_pdb

          select case (atom(i)%element)
          case ("H")
             radii(i)=1.0
          case ("C")
             radii(i)=1.7
          case ("N")
             radii(i)=1.37
          case ("O")
             radii(i)=1.16
          case ("S")
             radii(i)=1.8
          end select
          
          index=atom(1)%ressqnumber+CA_count
          select case (atom(i)%name)
          case (' N  ')
             j=i-1
             flag_N(index)=atom(i)%serialnumber
             flag_last(index-1)=atom(j)%serialnumber
             N(index)=atom(i)
!             print*, atom(i)%residue, flag_last(index-1),atom(i)%serialnumber
          case (' CA ')
             CA(index)=atom(i)
             flag_CA(index)=atom(i)%serialnumber
!             print*, atom(i)%residue, flag_CA(index),atom(i)%serialnumber,index
             CA_count=CA_count+1
             CA(index)%serialnumber=CA_count
             CA(index)%ressqnumber=CA_count
             
             if (CA(index)%x > max_coord(1)) max_coord(1) = CA(index)%x
             if (CA(index)%y > max_coord(2)) max_coord(2) = CA(index)%y
             if (CA(index)%z > max_coord(3)) max_coord(3) = CA(index)%z
             
             if (CA(index)%x < min_coord(1)) min_coord(1) = CA(index)%x
             if (CA(index)%y < min_coord(2)) min_coord(2) = CA(index)%y
             if (CA(index)%z < min_coord(3)) min_coord(3) = CA(index)%z
             
             write(13,pdb_format) CA(index)
          case (' C  ') 
             flag_C(index-1)=atom(i)%serialnumber
             if (atom(i)%residue .eq. 'GLY')  flag_CB(index-1)=atom(i+1)%serialnumber
!             print*, atom(i)%residue, flag_C(index-1),atom(i)%serialnumber
          case (' CB ')
             
             flag_CB(index-1)=atom(i)%serialnumber
!             print*, atom(i)%residue, flag_CB(index-1),atom(i)%serialnumber
          case default
             cycle
          end select
          
       end do
     
!     print*, flag_CB(:)
    
       first=1
       
    case (1)
        
       deallocate(CA)
       allocate(CA(start:CA_index))
              
!       write(12,pdb_format)atom(:)
!       write(12,"(A3)") "TER"
       
       do i=1,nlines_pdb
          
          index=atom(1)%ressqnumber+CA_count
          
          select case (atom(i)%name)
          case (' N  ')
             N(index)=atom(i)
             
          case (' CA ')
             CA(index)=atom(i)
             CA_count=CA_count+1
             CA(index)%serialnumber=CA_count
             CA(index)%ressqnumber=CA_count
             
             if (CA(index)%x > max_coord(1)) max_coord(1) = CA(index)%x
             if (CA(index)%y > max_coord(2)) max_coord(2) = CA(index)%y
             if (CA(index)%z > max_coord(3)) max_coord(3) = CA(index)%z
             
             if (CA(index)%x < min_coord(1)) min_coord(1) = CA(index)%x
             if (CA(index)%y < min_coord(2)) min_coord(2) = CA(index)%y
             if (CA(index)%z < min_coord(3)) min_coord(3) = CA(index)%z
             
             write(13,pdb_format) CA(index)
          case default
             cycle
          end select
          
       end do
    end select

    close(13)
    
  end Subroutine label
  
!     ##################### DISTANCE #######################
  Subroutine distance()
    integer  :: j=0,i
    real(dp) :: r,s,t

    
!    write(10,*) "Calculating the interatomic distances"
    allocate(dist(nlines_pdb,nlines_pdb))

    do i=1,nlines_pdb
       do j=1,nlines_pdb
          r=0
          s=0
          t=0
          r=atom(i)%x-atom(j)%x
          s=atom(i)%y-atom(j)%y
          t=atom(i)%z-atom(j)%z
          dist(i,j)=sqrt(r*r+s*s+t*t)
       end do
    end do
        
!    select case (int(d_ete))
!    case (15:40)
!       E_ete=0
!    case(:14)
!       E_ete=0.1*d_ete
!    case (41:)
!       E_ete=10*d_ete
!    end select
         
  end Subroutine distance

!     ######################### LJ #########################
  Subroutine  energy_LJ
    implicit none
    integer  :: i,k,l  !Controle
    real(dp) :: var=0.d0
    
!    write(10,*) "Calculating the Lennard-Jones energy"

    E_lj=0.d0
!    print *, E_lj    
    !     Todos resíduos não adjacentes
    do k= atom(1)%ressqnumber, (atom(nlines_pdb)%ressqnumber-2)
!      print *, k
       do l=(k+2), atom(nlines_pdb)%ressqnumber
!         PRINT *, CA(k)%residue, CA(l)%residue,l
          do i=flag_N(k), flag_last(k)
!            print *, flag_N(k), flag_last(k),i
             do j=flag_N(l), flag_last(l)
!               print *, flag_N(l), flag_last(l),j
                var=0
                
                if (dist(i,j) .lt. rvw(i,j)) then
                   var =  1!(((rvw(i,j)/dist(i,j))**12))! &
                   !                     & - ((rvw(i,j)/dist(i,j))**6))
!                   print "(3(F10.3))", dist(i,j), rvw(i,j), E_lj
                end if
                E_lj=E_lj+var
!                print *, flag_N(l), flag_last(l),j
             end do
          end do
       end do
    end do
    
!    print *, "Energia tds tds", E_lj
    
    !     Resíduos adjacentes
    do k=atom(1)%ressqnumber,(atom(nlines_pdb)%ressqnumber-1)
       l=k+1
!       PRINT *, CA(k)%residue, CA(l)%residue  
       do i=flag_N(k), flag_CB(k)
!          print *, flag_N(k), flag_CB(k)
          do j=flag_N(l), flag_CB(l)
!             print *, flag_N(l), flag_CB(l)
             if (dist(i,j) .lt. rvw(i,j)) then
                var = 1 !(((rvw(i,j)/dist(i,j))**12))!-((rvw(i,j)/dist(i,j))**6))
                
                if ((atom(i)%name .eq. ' CA ') .and. &
                     &(atom(j)%name .eq. ' N  ')) then
                   var=0
                   cycle
                end if
                
                if ((atom(i)%name .eq. ' CA ') .and. &
                     &(atom(j)%name .eq. ' CA ')) then
                   var=0
                   cycle
                end if

                if ((atom(i)%name .eq. ' C  ') .and. &
                     &(atom(j)%name .ne. ' O  ')) then
                   var=0
                   cycle
                end if
                
                if ((atom(i)%name .eq. ' O  ') .and. &
                     &(atom(j)%name .eq. ' N  ')) then
                   var=0
                   cycle
                end if
                
                if (((atom(i)%name .eq. ' O  ') .and. &
                     &(atom(j)%name .eq. ' CA '))) then
                   var=0
                   cycle
                end if
                
                if ((atom(i)%name .eq. ' CB ') .and. &
                     &(atom(j)%name .eq. ' N ')) then
                   var=0
                   cycle
                end if
                E_lj=E_lj+var
             end if
          end do
       end do
    end do
!    print *, "Energia adjacentes", E_lj

!    print *, E_lj
  end Subroutine energy_LJ
!     ####################### HBOND ########################
  Subroutine hbond
    implicit none
    integer,allocatable :: atm1(:),atm2(:),atm3(:),atm5(:),atm4(:)
    integer  :: i,k
    real(dp) :: r1(5),r2(5),r3(5),eta,theta,r,s,Falpha,sub,times
    real(dp) :: v1,v1x,v1y,v1z,v2,v2x,v2y,v2z,v3,v3x,v3y,v3z,pt999=0.9990d0
    real(dp) :: alpha(3),lambda,beta(3),mi(3),max1,max2,min1,min2
    
    E_hbond=0.d0
    allocate(atm1(nhbond),atm2(nhbond),atm3(nhbond),atm4(nhbond),atm5(nhbond))
    
    i=1
    open(13,file="hbonds.dat",status="old")
    
    do while (.true.)
!      print *, i
       read(13,8,end=15)atm2(i),atm3(i),atm4(i),atm1(i),atm5(i)
       i=i+1
       
    end do
15  continue
    
    close(13)
      
!      PRINT *, nhbond
      
    i=1
    do i=1, nhbond
       
       !     Acceptor coordinates
       r1(1)=atom(atm1(i))%x
       r2(1)=atom(atm1(i))%y
       r3(1)=atom(atm1(i))%z
       
       !     Donor coordinates
       r1(2)=atom(atm2(i))%x
       r2(2)=atom(atm2(i))%y
       r3(2)=atom(atm2(i))%z
       
       !     C donor coordinates
       r1(3)=atom(atm3(i))%x
       r2(3)=atom(atm3(i))%y
       r3(3)=atom(atm3(i))%z
       
       !     C-alpha donor coordinates
       r1(4)=atom(atm4(i))%x
       r2(4)=atom(atm4(i))%y
       r3(4)=atom(atm4(i))%z
       
       !     C-alpha acceptor coordinates
       r1(5)=atom(atm5(i))%x
       r2(5)=atom(atm5(i))%y
       r3(5)=atom(atm5(i))%z
              
       !     Distance Check
       
       v1x=r1(2)-r1(1)
       v1y=r2(2)-r2(1)
       v1z=r3(2)-r3(1)
       v1=sqrt((v1x**2)+(v1y**2)+(v1z**2))
       
       alpha(1)=v1
       beta(1)=100
       mi(1)=3.0
       
       !     if (alpha(1) .gt. mi(1)) cycle
       
       !     First Angle check
       v2x=r1(3)+r1(4)-r1(2)-r1(2)
       v2y=r2(3)+r2(4)-r2(2)-r2(2)
       v2z=r3(3)+r3(4)-r3(2)-r3(2)
       v2=sqrt((v2x**2)+(v2y**2)+(v2z**2))
       
       r=(r1(1)*r1(2))+(r2(1)*r2(2))+(r3(1)*r3(2))
       eta=cos(r/(v1*v2))
       
       max1=max(-pt999,eta)
       min1=min(pt999,max1)
       
       alpha(2)=acos(min1)
       beta(2)=100
       mi(2)=0.5
       
       !     if (alpha(2) .gt. mi(2)) cycle
       
       !     Second Angle check
       v3x=r1(1)-r1(5)
       v3y=r2(1)-r2(5)
       v3z=r3(1)-r3(5)
       v3=sqrt((v3x**2)+(v3y**2)+(v3z**2))
       
       s=(r1(1)*r1(3))+(r2(1)*r2(3))+(r3(1)*r3(3))
       theta=cos(s/(v1*v3))
       
       max2=max(-pt999,theta)
       min2=min(pt999,max2)
       
       alpha(3)=acos(min2)
       beta(3)=100
       mi(3)=0.6
       
       !     if (alpha(3) .gt. mi(3)) cycle
       
       !         print *, 'VALORES: ', alpha(1),alpha(2),alpha(3)
       
       if (alpha(1) .lt. mi(1) .and. alpha(2) .lt. mi(2) &
            &.and. alpha(3) .lt. mi(3)) then
          lambda=1.0d0
          cycle
       endif
       
       Falpha=0.d0
       
       do k=1,3
          if (alpha(k) .lt.mi(k)) then
             lambda=1.0d0
          else
             sub=alpha(k)-mi(k)
             times=beta(k)*sub
             Falpha=1/(1+exp(times))
             lambda=nint(lambda*Falpha)
          endif
!          WRITE(100,*) sub,times,Falpha,lambda
       end do
       
       !         PRINT *,'FUNÇÃO CALCULADA: ', lambda, i,j
       
       E_hbond=E_hbond-lambda
       
    end do
    
!      PRINT *, 'Energia: ',E_hbond
    
8   format(5I7)
    
  end Subroutine hbond
  !     ###################### ACC-DON #######################
  Subroutine accdon
    implicit none
    integer :: i,j,k,RES3,p

    nhbond=1
    open(12,file="hbonds.dat",status="unknown")
    do i=1,nlines_pdb-1
       do p=1,nlines_pivot
!            print *, list(p), p
          if (atom(i)%ressqnumber .eq. list(p)%residue) then
             if(atom(i)%name .eq.' N  ' .and. atom(i+1)%name .eq.' CA ' &
                  & .and. atom(i)%ressqnumber .gt.atom(1)%ressqnumber .and. &
                  & atom(i)%residue .ne.'PRO') then
                do j=1,i
                   if(atom(j)%name .eq.' C  ' .and. atom(j)%ressqnumber &
                        & .eq. atom(i)%ressqnumber-1) then
                      do k=2,nlines_pdb
                         RES3=ABS(atom(i)%ressqnumber-atom(k)%ressqnumber)
                         !     RES3 define the distance between the residues 
                         !     I think should be 3, but the default in the original code is 1
                         if(atom(k)%name.eq.' O  '.and.atom(k-1)%name .eq.' C  ' &
                              & .and. RES3.gt.3) then
                            write(12, '(5I7)') atom(i)%serialnumber,atom(i+1)%serialnumber, &
                                 & atom(j)%serialnumber,atom(k)%serialnumber, &
                                 & atom(k-1)%serialnumber
                            nhbond=nhbond+1
                         end if
                      end do
                   end if
                end do
             end if
          end if
       end do
    end do
    
    close(12)
  end Subroutine accdon

!     ##################### WRITE PDB ######################
  Subroutine writepdb(flag_final,m)
    integer, intent(in) :: flag_final,m
    character :: C
    
!    write(10,*) "Writing output file"    
    
    if (flag_final .eq. 1) then
       PDB_output="ressub.pdb"
       open(42,file=PDB_output,status='unknown')
       write(42,pdb_format)atom(1:)
       write(42,'(A3)')'END'
       close (42)
    end if
    
    Select case (m)
    case (:9) 
       write(PDB_output, '( "conformation", I1, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I1, ".dat" )' )  m
    case (10:99) 
       write(PDB_output, '( "conformation", I2, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I2, ".dat" )' )  m
    case (100:999)
       write(PDB_output, '( "conformation", I3, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I3, ".dat" )' )  m
    case (1000:9999)
       write(PDB_output, '( "conformation", I4, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I4, ".dat" )' )  m
    case (10000:99999)
       write(PDB_output, '( "conformation", I5, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I5, ".dat" )' )  m
    case (100000:999999)
       write(PDB_output, '( "conformation", I6, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I6, ".dat" )' )  m
       write(10,'(A)') new_line(C),"Too many outputs for this run. Reset."
    case (1000000)
       write(PDB_output, '( "conformation", I7, ".pdb" )' )  m
       write(EXP_output, '( "conformation", I7, ".dat" )' )  m
       write(10,'(A)') new_line(C), "The run reached output limit"
    case default
       write(10,'(A)') new_line(C), "Something is wrong!! Check."
       stop "Error detected"
    end Select

!    print *, PDB_output, EXP_output
!    
    open(42,file="gene.pdb",status='unknown')
    write(42,pdb_format)atom(1:)
    write(42,'(A3)')'END'
    close (42)

    open(41,file="gene.dat",status='unknown')
    write(41,exp_format)theo_data(1:)
    write(41,'(A3)')'END'
    close (41)
    
    open(42,file=PDB_output,status='unknown')
    write(42,pdb_format)atom(1:)
    write(42,'(A3)')'END'
    close (42)
  End Subroutine writepdb
     
end Module pdb_file
