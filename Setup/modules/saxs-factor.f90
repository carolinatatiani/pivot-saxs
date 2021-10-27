!     ######################################################
!     Protein solvation and transformation into c-alpha pdb
!    
!     #######################################################
!     #              ##########################             #
!     ############### Sintaxe: ./solvate.e pdb ##############
!     #              ##########################             #        
!     #######################################################
!
!     Carolina Tatiani Alves Ferreira - carolina.tatiani@unesp.br
 
program saxs_factor
  use variables
  use pdb_file
  use exp_file
  use water_box
  use saxs_theo
  implicit none
  integer       :: flag_final=0,m=0
  integer       :: index,CA_count,water_count
  character(40) :: arg,PDB_CA,char2num!,logexp
  real(dp)      :: protsize(3)=0,Rg,Rg_theo
  
  open(10, file="log.out",action='write',position='append')
  i = 0
  do
     call getarg(i, arg)
     if (len_trim(arg) == 0) exit    
     select case (arg)
     case ("-pdb")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) then
           write(10,"(A)") "Something is wrong! I've got no pdb file name."
           stop "Error detected"
        end if
        
        !write (*,*) trim(arg)
        PDB_input=trim(arg)
     case ("-exp")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) then
           write(10,"(A)") "Something is wrong! I've got no experimental file name."
           stop "Error detected"
        end if

        !write (*,*) trim(arg)
        EXP_input=trim(arg)
     case ("-rg")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) cycle
        !write (*,*) trim(arg)
        char2num=trim(arg)
        read(char2num,*)Rg
     end select
     i=i+1         
  end do
    
  call readpdb()
  !print*, "PDB read."
  close(10)
  open(10, file="log.out",action='write',position='append')  
  call readexp() 
  !print*, "EXP read."  
  close(10)
  open(10, file="log.out",action='write',position='append')
  
  call label(PDB_CA,index,CA_count)
  !print*, "Protein labeled."
  close(10)
  open(10, file="log.out",action='write',position='append')
  
  protsize(1)=max_coord(1)-min_coord(1)
  protsize(2)=max_coord(2)-min_coord(2)
  protsize(3)=max_coord(3)-min_coord(3)

  allocate(water(1))
  call  solvate(protsize,min_coord,CA_count,water_count)
  !print*, "Protein solvated."
  
  close(10)
  open(10, file="log.out",action='write',position='append')
  deallocate(CA)  
  call readCA(CA_count,PDB_CA)
  allocate(factor(1:CA_count,1:nlines_exp))
  
  close(10)
  open(10, file="log.out",action='write',position='append')

  call distance()

  !print*, "Interatomic distances calculated."
  
  close(10)
  open(10, file="log.out",action='write',position='append')
    
!  d_ete=dist(end_C,end_N)
  !print exp_format, exp_data(:)
!  allocate(log_exp(1:nlines_exp))    

!  log_exp%q=exp_data%q
!  log_exp%intense=log10(exp_data%intense)
!  log_exp%error=log10(exp_data%error)
  
  !  write(logexp, '( "log_", A36 )' ) EXP_input
  !  open(30, file=logexp, status='Unknown') !Input
  !  write(30, exp_format) log_exp
  !  close(30)
  
  !print*, CA_count,water_count

  !  print exp_format, exp_data(:)
  first=0
  call theo(CA_count,Rg,water_count,Rg_theo)
  
  close(10)
  open(10, file="log.out",action='write',position='append')
    
  write(10,"(A,1X,F10.3)") 'Experimental gyration radius = ', Rg
  write(10,"(A,1X,F10.3)") 'Theoretical gyration radius = ', Rg_theo
  write(10,"(A, f10.3, A1)") "Theoretical profile calculated with chi=", Chi,"."
  
  close(10)
  
end program saxs_factor


!rm *.mod *.e clean_* conformation* log.out ressub.* CA_* theo.dat fort.60 #Clean diretory
