Module water_box
  use variables
  implicit none
  real(dp), parameter :: water_box_size=119.7, thickness=3.d0, close_dist=3.5
  integer, parameter  :: n_waters=60655
  
contains
  !     ###################### SOLVATE #######################
  Subroutine solvate(protsize,min_coord,CA_count,water_count)
    implicit none
    
    real(dp), intent(in)   :: protsize(3),min_coord(3)
    real(dp)               :: box_size(3),dist,r,s,t
    real(dp)               :: shell,min
    integer                :: i,j,k,m=0,n=0,CA_count,water_count
    character(40)          :: PDB_CA
!    write(10,*) "Protein solvation"
    i=0
    j=0
    k=0
    water_count=0
    
    box_size(1)=(protsize(1)/water_box_size)+1
    box_size(2)=(protsize(2)/water_box_size)+1
    box_size(3)=(protsize(3)/water_box_size)+1

    open(11, file="water_box.pdb", status="old")
    write(PDB_CA, '( "CA_", A37 )' )  PDB_input
    
    open(13, file=PDB_CA,action='write',position='append')
    
    do while (i < box_size(1))
       do while (j < box_size(2))
          do while (k < box_size(3))
             do m=1,n_waters
              
                read(11,pdb_format) water

                water%x=water%x+min_coord(1)+(i*water_box_size)
                water%y=water%y+min_coord(2)+(j*water_box_size)
                water%z=water%z+min_coord(3)+(k*water_box_size)
                
                min=1000000
                do n=1,CA_count
                   
                   r=CA(n)%x-water(1)%x
                   s=CA(n)%y-water(1)%y
                   t=CA(n)%z-water(1)%z
                   dist=sqrt(r*r+s*s+t*t)
                   if (dist < 3.0) then
                      min=0
                      exit
                   end if
                   if (dist < min) min=dist
                end do
                
                if (min < close_dist) cycle
                
                shell=close_dist+thickness
                
                if (min > shell) cycle
                
                water_count=water_count+1
              
                water%serialnumber=CA_count+water_count
                water%ressqnumber=CA_count+water_count
                write(13,pdb_format) water
                
             end do
             
             rewind 11
             k=k+1
          end do
          k=0
          j=j+1
       end do
       j=0
       i=i+1
    end do

    CA_count=CA_count+water_count
    close(11)
    
    write(13, "(A3)") "END"
    close(13)
    
    
  end Subroutine solvate
end Module water_box

