Module exp_file
  use variables
  implicit none
  integer                :: i
  character(40)          :: cleanexp
   
contains
!     ###################### READ EXP ######################
  Subroutine readexp
    implicit none
    integer       :: ierror=999,k,num=1
    character(80) :: junk
    real          :: ncol(3)=0

    i=0
    nlines_exp=0
    write(10,*) "Reading EXP file"    
    open(12, file=EXP_input, status='old') !Input
   
    write(cleanexp, '( "clean_", A34 )' ) EXP_input

    open(22, file=cleanexp, status='Unknown') !Input

    do while (.true.)
       read(12,'(I1)',iostat=ierror) num
       select case (ierror)
       case (0) 
          backspace 12
          read(12,'(A)') junk
          call getvals(junk,ncol,k)
          if ((k .eq. 2) .or. (ncol(3) .eq. 0)) ncol(3)=ncol(1)*ncol(2)
           write(22,exp_format) ncol
          i=i+1
       case (:-1)
         exit
       case default
          backspace 12
          read(12,'(A)') junk
       end select
    end do
    close(12)

    rewind 22

    nlines_exp=i !Defining the number of lines and allow to use i again

!    write(10,*) "The file contains ",k, " columns, and ", nlines_exp," lines."

    allocate(exp_data(0:nlines_exp),theo_data(0:nlines_exp))
    exp_data(0)%q=0.d0
    
    do i=1,nlines_exp
       read(22,exp_format) exp_data(i)
    end do
    close(22)
  end Subroutine readexp

!     ###################### GETVALS ####################### 
  subroutine getvals(line,values,icount,ierr)
    implicit none
!       character(len=*),parameter  :: ident='@(#)getvals:read arbitrary number of values from a &
!            &character variable up to size of values'
    ! JSU 20170831
    
    character(len=*),intent(in)  :: line
    real                         :: values(:)
    integer,intent(out)          :: icount
    integer,intent(out),optional :: ierr
    
    character(len=len(line)+1)   :: buffer
    character(len=len(line))     :: words(size(values))
    integer                      :: ios, i, ierr_local
    
    ierr_local=0

    words=' '                            ! make sure words() is initialized to null+blanks
    buffer=trim(line)//"/"               ! add a slash to the end so how the read behaves with missing values is clearly defined
    read(buffer,*,iostat=ios) words      ! undelimited strings are read into an array
    icount=0
    do i=1,size(values)                  ! loop thru array and convert non-blank words to numbers
       if(words(i)(1:1).eq.' ')exit
       read(words(i),*,iostat=ios)values(icount+1)
       if(ios.eq.0)then
          icount=icount+1
       else
          ierr_local=ios
          write(10,*)'*getvals* WARNING:['//trim(words(i))//'] is not a number'
       endif
    enddo
    
    if(present(ierr))then
       ierr=ierr_local
    elseif(ierr_local.ne.0)then        ! error occurred and not returning error to main program to print message and stop program
       write(10,*)'*getval* error reading line ['//trim(line)//']'
       stop 2
    endif
    
  end subroutine getvals
  
End Module exp_file
