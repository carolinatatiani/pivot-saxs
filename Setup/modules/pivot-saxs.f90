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

program pivot_saxs
  use variables
  use pivot
  use pdb_file
  use exp_file
  use water_box
  use saxs_theo
  implicit none
  integer  :: flag_final=0,m=0,config=0
  integer  :: index,CA_count,water_count
  real(dp) :: protsize(3)=0,kT=0.d0,E_old=100000
  real(dp) :: angrot,randm(3),E_total
  integer  :: step,flag=0,config_m
  real(dp) :: chi_m,rnd,Rg=0.0,Rg_theo!,dist_m
  real(dp) :: hbond_m,lj_m,rest,E_rg_m,acc,hb
  character(40) :: PDB_CA,arg,char2num!,logexp

  open(10, file="log.out",action='write',position='append')
  open(40, file="Energies.out",status="unknown")

  bin=10
  
  i = 0
  call readlist(m,config)
  
  do
     call getarg(i, arg)
     if (len_trim(arg) == 0) exit
     select case (arg)
     case ("-pdb")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) then
           write(10,"(A)") "Something is wrong! I've got no pdb file name."
           stop "Error detected. I've got no pdb file name."
        end if
        PDB_input=trim(arg)
     case ("-exp")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) then
           write(10,"(A)") "Something is wrong! I've got no experimental file name."
           stop "Error detected. I've got no experimental file name."
        end if
        EXP_input=trim(arg)
     case ("-rg")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) cycle
        write (*,*) trim(arg)
        char2num=trim(arg)
        read(char2num,*)Rg
     case ("-min")
        i = i+1
        call getarg(i, arg)
        if (len_trim(arg) == 0) stop ("Something is wrong! I have got 0 minimization steps.")
        !write (*,*) trim(arg)
        char2num=trim(arg)
        read(char2num,*)n_open
     case ("-hb")
        i = i+1
        call getarg(i, arg)
        char2num=trim(arg)
        read(char2num,*)eta_bond
     case ("-lj")
        i = i+1
        call getarg(i, arg)
        char2num=trim(arg)
        read(char2num,*)eta_lj
     case ("-chi")
        i = i+1
        call getarg(i, arg)
        char2num=trim(arg)
        read(char2num,*)eta_chi  
     end select
     i=i+1         
  end do
    
!  open(10, file="log.out",action='write',position='append')
  call readpdb()
  !print*, "PDB read."
 ! close(10) 
  !open(10, file="log.out",action='write',position='append')
  
  call label(PDB_CA,index,CA_count)

  flag_last(atom(nlines_pdb)%ressqnumber)=atom(nlines_pdb)%serialnumber
  
!  end_N=N(atom(1)%ressqnumber)%serialnumber
!  end_C=flag_C(index-1)
  
  do i=1,nlines_pdb
     do j=1,nlines_pdb
        rvw(i,j)=radii(i)+radii(j)
     end do
  end do

  call distance()
  call energy_lj
  deallocate(dist)
  
  !print*, "Protein labeled."

  kT=kb*temp_max
  !write(10,"(A)") "Minimizing Lennard-Jones."

  !print*, n_open

  E_old=E_lj

  
  !print pdb_format, atom(1:)
  !print*, n_open
  
  do step=1,n_open
     call label(PDB_CA,index,CA_count)
    
     do i=1,3
        call generator(rnd,seed)
        randm(i)=rnd
     end do

     i=int(nlines_pivot*randm(1))+1 !Escolha aleatória do resíduo
     resrot=list(i)%residue
     
     angrot=(nint(12*randm(2)-6))*5 !Escolha aleatória do angulo de rotação
     
     dihedral=nint((randm(3))+1) !1 para psi 2 para phi

     !print*, resrot, angrot,dihedral

     call rotate(angrot)
     call distance()
     call energy_lj
     deallocate(dist)
     
     call metropolis(E_old,E_lj,kT,flag,seed)
!     if (flag .eq. 1) then
        E_old=E_lj
        open(20, file="open.pdb",status="unknown")  
        write(20, pdb_format) atom(1:)
        write(20, "(A3)") "END"
        close (20)
!     end if
     
  end do
  
  open(20, file="open.pdb",status="unknown")  
  write(20, pdb_format) atom(1:)
  write(20, "(A3)") "END"
  close (20)

  
  deallocate(atom)

    !print*, "Lennard-Jones minimized,", E_old
  
  PDB_input="open.pdb"
!  close(10)
!  open(10, file="log.out",action='write',position='append')
  call readpdb()
  !print*, "PDB read."
  call accdon()
  !print*, "Donors and acceptors calculated."
  !close(10) 
!  open(10, file="log.out",action='write',position='append')
  call readexp()
  !print*, EXP_input
  !print*, "EXP read."
!  close(10)
  
!  open(10, file="log.out",action='write',position='append')
  
  call label(PDB_CA,index,CA_count)
  !print*, "Protein labeled."
!  close(10)
!  open(10, file="log.out",action='write',position='append')
 
  call distance()

  !print*, "Interatomic distances calculated."
    
  protsize(1)=max_coord(1)-min_coord(1)
  protsize(2)=max_coord(2)-min_coord(2)
  protsize(3)=max_coord(3)-min_coord(3)

  allocate(water(1))   
  call solvate(protsize,min_coord,CA_count,water_count)
  !print*, "Protein solvated."
!  close(10)
!  open(10, file="log.out",action='write',position='append')
  deallocate(CA)  
  call readCA(CA_count,PDB_CA)
  !print*, "CA file read."
  allocate(factor(1:CA_count+water_count,1:nlines_exp))

  first=0
  factor(:,:)=0.d0   
  call theo(CA_count,Rg,water_count,Rg_theo)
  !print*, 'chi=', chi,'rg=',Rg
  !print*, "Theoretical profile calculated."
  first=1
  !  close(10)
  !  open(10, file="log.out",action='write',position='append')
  !   
  !  write(10,"(A37, f10.3, A1)") "Theoretical profile written and chi=", Chi,"."
  
!  close(10)
!  open(10, file="log.out",action='write',position='append')
  
  call energy_LJ
!  print*, "Lennard-Jones calculated."
!  call hbond
!  hb=E_hbond
!   print*, "Hidrogen bond calculated."
  E_old= 1000!(eta_lj*E_lj)+(eta_chi*chi)+(eta_rg*E_rg)+(eta_bond*E_hbond)
  write(40,'(5(A10, 1X))') "E_LJ","E_chi","E_rg","E_hbond","E_total"
!  print*, E_old
  
  !n_passos=nint(bin*turn*(temp_max-temp_min))
!  print*, n_passos
  seed=1
  tun=0.d0
  temp=temp_max
  
  do step=1,n_passos
     !print*, step, temp
     !print'(5(F10.3,1X))',E_lj,chi,E_rg,E_hbond,E_total
     rest=mod(step,bin)
     
     if (rest .eq. 0) then
        temp=temp-minus
        kT=kb*temp
        !     write(*,*)temp,minus
        if ((temp .eq. temp_max) .or. (temp .eq. temp_min)) then
           minus=minus*(-1)
           tun=tun+0.5
           eta_chi=eta_chi+0.25
           if (tun .eq. turn) exit
        end if
     end if
 
     do i=1,3
        call generator(rnd,seed)
        randm(i)=rnd
     end do
     
     i=int(nlines_pivot*randm(1))+1 !Escolha aleatória do resíduo
     resrot=list(i)%residue
     
     angrot=nint(60*randm(2))-30 !Escolha aleatória do angulo de rotação
     
     dihedral=nint((randm(3))+1) !1 para psi 2 para phi
     
     !print*, resrot,angrot,dihedral
     
     call rotate(angrot)
!     close(10)
!     open(10, file="log.out",action='write',position='append')
     !print*, "New configuration generated."
     
     config=config+1
     
     call label(PDB_CA,index,CA_count)
     !print*, "Protein labeled."
!     close(10)
!     open(10, file="log.out",action='write',position='append')
     
     deallocate(dist)
     call distance()
     
     !print*, "Interatomic distances calculated."
     
     protsize(1)=max_coord(1)-min_coord(1)
     protsize(2)=max_coord(2)-min_coord(2)
     protsize(3)=max_coord(3)-min_coord(3)
     
     call solvate(protsize,min_coord,CA_count,water_count)
     !print*, "Protein solvated."
     !close(10)
     !open(10, file="log.out",action='write',position='append')
     deallocate(CA)  
     call readCA(CA_count,PDB_CA)
     !print*, "CA file read."
     
     call theo(CA_count,Rg,water_count,Rg_theo)
     !print*,'new chi=', chi,'new rg=', Rg
     !print*, "Theoretical profile calculated."
    ! close(10)
    ! open(10, file="log.out",action='write',position='append')
     
     !    write(10,"(A37, f10.3, A1)") "Theoretical profile written and chi=", Chi,"."
     !    close(10)
     !    open(10, file="log.out",action='write',position='append')
     
     call energy_LJ
     !print*, "Lennard-Jones calculated."
     !call hbond
     !E_hbond=E_hbond-hb
     !print*, "Hidrogen bond calculated."
     E_total=(eta_lj*E_lj)+(eta_chi*chi)+(eta_rg*E_rg)+(eta_bond*E_hbond)
     
     call metropolis(E_old,E_total,kT,flag,seed)
     !print'(5(F10.3,1X),I1)', E_lj,chi,E_rg,E_old,E_total,flag
     if (flag .eq. 1) then
        E_old=E_total
        chi_m=chi
        lj_m=E_lj
        config_m=config
        E_rg_m=E_rg
        hbond_m=E_hbond

        m=m+1
        !       print '(4(F10.3,1X))', E_lj,chi,E_hbond, E_total
        !       write(10,'(4(F10.3,1X))')E_lj,chi,E_hbond, E_total
        write(40,'(5(F10.3,1X))')E_lj,chi_m,E_rg,hbond_m,E_total
        
        call writepdb(flag_final,m)
        
        open (30, file=EXP_output, status='unknown')
        write(30, exp_format) theo_data(1:)
        close(30)
     end if
     if (chi .le. (0.2*chi_m)) exit
     !    print "(E10.3, 1X, I1)", chi, flag
  end do
  
  flag_final=1
  call writepdb(flag_final,m)
  open(60, file="ressub.out",status='unknown')
  write(60,"(f10.3, 3I15)") turn, m, n_passos, seed
  write(60,"(3I10)") temp_min, temp_max, minus
 ! write(60,"(4f10.3)") eta_lj, eta_chi, eta_rg, eta_bond
  
  acc=m/config  
  write(10,'(A16,F6.3)') "Acceptance rate=", acc
  write(10,'(5(A10, 1X))') "E_LJ","E_chi","E_rg","E_hbond","E_total"
!  print '(5(F10.3,1X))', lj_m,chi_m,E_rg_m,hbond_m, E_old
  write(10,'(5(F10.3,1X))') lj_m,chi_m,E_rg_m,hbond_m, E_old
  
  close(10)
  close(40)
  close(60)
  
end program pivot_saxs


!     ####################  SUBROUTINE #####################

      
!     ####################### RANDOM #######################
Subroutine generator(ran2,idum) !NUMERICAL RECIPIES
  integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  real*8 ran2,AM,EPS,RNMX
  parameter (IM1=2147483563,IM2=2147483399,AM=1./real(IM1),IMM1=IM1-1, &
       &IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
       &IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  
  integer idum2,j,k,iv(NTAB),iy
  save iv,iy,idum2
  data idum2/123456789/, iv/NTAB*0/, iy/0/
  
  if (idum.le.0) then       
     idum=max(-idum,1)      
     idum2=idum
     do j=NTAB+8,1,-1    
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     end do
     iy=iv(1)
  endif
  k=idum/IQ1 
  idum=IA1*(idum-k*IQ1)-k*IR1 
  if(idum.lt.0) idum=idum+IM1 
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2 
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV 
  iy=iv(j)-idum2 
  iv(j)=idum 
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX) 

end Subroutine generator

!     ##################### METROPOLIS #####################
Subroutine metropolis(E_old,E_total,kT,flag,seed)
  integer,parameter    :: dp = selected_real_kind(15, 307)
  integer, intent(out) :: flag
  real(dp), intent(in) :: E_total,E_old,kT
  real(dp)             :: random,delta,Pb
  integer              :: seed
    
  delta=E_total-E_old
  flag=0
  
  if (delta .lt. 0.d0) flag=1
  
  if (delta .ge. 0.d0) then

     call generator(random,seed)
     Pb=exp(-delta/kT)
     if (random .lt. Pb) flag=1
     
  end if

end Subroutine metropolis

!    ######################################################
!rm *.mod *.e clean_* conformation* log.out ressub.* CA_* theo.dat fort.80 hbonds.dat Energies.out *# *~
