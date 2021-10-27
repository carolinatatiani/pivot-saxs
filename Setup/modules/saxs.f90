Module saxs_theo
  use variables
  implicit none
  
contains
  !    ####################### THEO #######################
  Subroutine theo(CA_count,Rg,water_count,Rg_theo)
    real(dp), parameter   :: pi=3.1415,grau=(1+sqrt(5.0))/2
    real(dp)              :: arg1=0,arg2=0,qx,qy,qz,M,N,qdotr
    real(dp)              :: offset,var1,var2,var
    real(dp)              :: qb,qa,Rg_theo,Rg
    real(dp)              :: tablef(7)
    integer               :: j=0,i=0,k=0,nstop,nstart,l
    integer, intent(in)   :: CA_count,water_count

!    write(10,*) "Calculating the theoretical intensity"
    
    l=201
    nstart=-(l-1)/2
    nstop=(l-1)/2
    
    theo_data%intense=0.d0
    theo_data%error=0.d0
    theo_data%q=exp_data%q
    
    i=(CA_count-water_count+2)

                

!    write(100,pdb_format) CA(:)
    
    do k=1,nlines_exp
!       print*, i,k
       if (first .eq. 0) then
          do i=1,(res_count+1)
             call ff_table(CA(i)%residue,tablef)
             factor(i,k)=ff(tablef(1),tablef(2),tablef(3),tablef(4), &
                  &tablef(5),tablef(6),tablef(7),exp_data(k)%q)
             !print*, ca(i)%residue,factor(i,k)
          end do
       else
          i=(res_count+1)
       end if
       factor(i:,k)=factor(i-1,k)
 !      print*, factor(:,k)
 !      print*, k
    end do

    !print'(F8.3)', factor(1,k-1)
    !print'(F8.3)', factor(i,k-1)
    do k=1,nlines_exp
       !       print *, k
       do j=nstart,nstop
          
          arg1= (2.0*j)/l
          arg2= (2.0*pi*j)/grau
          
          qx = theo_data(k)%q * cos(asin(arg1)) * cos(arg2)
          qy = theo_data(k)%q * cos(asin(arg1)) * sin(arg2)
          qz = theo_data(k)%q * arg1
          
          M=0.d0 
          N=0.d0
          
          do i=1,CA_count
             
             qdotr=qx*CA(i)%x+qy*CA(i)%y+qz*CA(i)%z
             
             M = M + factor(i,k)*cos(qdotr)
             N = N + factor(i,k)*sin(qdotr)
             
          end do
          
          M=M*M
          N=N*N                                    
          
          theo_data(k)%intense = theo_data(k)%intense + M + N
          
       end do
!       print*, k
       !print exp_format , theo_data(k)
    end do

    !write(100,exp_format)theo_data
    theo_data(:)%intense = (theo_data(:)%intense) / l
    !write(200,exp_format)theo_data
    
    qa=0.035
    qb=0.04
    if (Rg .eq.0.d0)  Rg=RgwI(qa,qb,exp_data%q,exp_data%intense)
    
    qb=1.3/Rg
    
    if (qb .le. qa) qa=theo_data(nlines_exp)%q*0.01 
    
    k=0
    var1=0.d0
    var2=0.d0
    exp_data(0)%q=0.d0
    
    do while (exp_data(k)%q .le. qb)
       k=k+1
       if (exp_data(k)%q .le. qa) cycle
       
       var=0.d0
       var=((theo_data(k)%intense)*(exp_data(k)%intense))/(exp_data(k)%error**2)
       var1=var1+var
       
       var=0.d0
       var=((exp_data(k)%intense)**2)/((exp_data(k)%error)**2)
       var2=var2+var
    end do
    
    offset=var2/var1
    !print*, var2, var1
    !print*, 'offset=', offset
    theo_data(:)%intense = theo_data(:)%intense * offset
    
    theo_data(:)%error= 0.1*theo_data(:)%intense
    
    chi=0.d0
    do k=1,nlines_exp
       var=0.d0
       var=((theo_data(k)%intense)-(exp_data(k)%intense))**2
       chi=chi+(var/(exp_data(k)%error**2))
    end do
    
    Rg_theo=RgwI(qa,qb,theo_data%q,theo_data%intense)
    chi=chi/(k-1)
    
    !write(10,"(A,1X,F10.3)") 'Theoretical gyration radius = ', Rg_theo
    open(20, file="theo.dat", status="unknown")
    write(20, exp_format) theo_data(1:)
    close(20)

    E_rg=sqrt((Rg-Rg_theo)**2)

    !write(300,exp_format) exp_data(1:)
    
  end Subroutine theo
  
  !     #################### FF_TABLE ######################
  Subroutine ff_table(resname,tablef)
    character(3),intent(in) :: resname
    real(dp),intent(out)    :: tablef(:)
    
    tablef(:)=0.d0
    
    select case (resname)
    case ("GLY")
       
       tablef(1)=  9.97297e+00
       tablef(2)= -6.14483e-02
       tablef(3)=  4.11083e+01
       tablef(4)= -4.31734e+02
       tablef(5)=  2.07260e+03
       tablef(6)= -4.44561e+03
       tablef(7)=  3.42549e+03
       
    case ("ALA")
       
       tablef(1)=  9.06253e+00
       tablef(2)= -5.79704e+00
       tablef(3)=  2.03678e+02
       tablef(4)= -2.09348e+03
       tablef(5)=  8.17094e+03
       tablef(6)= -1.36443e+04
       tablef(7)=  8.30723e+03
       
    case ("VAL")
       
       tablef(1)=  7.25180e+00
       tablef(2)= -1.22282e+01
       tablef(3)=  2.39638e+02
       tablef(4)= -2.40334e+03
       tablef(5)=  8.70960e+03
       tablef(6)= -1.24978e+04
       tablef(7)=  5.99891e+03
       
    case ("ILE")
       
       tablef(1)=  6.30945e+00
       tablef(2)= -1.22530e+01
       tablef(3)=  2.83650e+02
       tablef(4)= -3.27665e+03
       tablef(5)=  1.30824e+04
       tablef(6)= -2.09082e+04
       tablef(7)=  1.16636e+04
       
    case ("LEU")
       
       tablef(1)=  6.30477e+00
       tablef(2)= -8.63718e+00
       tablef(3)=  1.80797e+02
       tablef(4)= -2.21694e+03
       tablef(5)=  8.92552e+03
       tablef(6)= -1.39515e+04
       tablef(7)=  7.45445e+03
       
    case ("MET")
       
       tablef(1)=  1.66478e+01
       tablef(2)= -1.93055e+01
       tablef(3)=  5.69130e+02
       tablef(4)= -5.99367e+03
       tablef(5)=  2.39864e+04
       tablef(6)= -4.13177e+04
       tablef(7)=  2.59618e+04

    case ("PHE")
       
       tablef(1)=  9.16913e+00
       tablef(2)=  1.69064e-01
       tablef(3)=  5.53998e+01
       tablef(4)= -1.73820e+03
       tablef(5)=  8.89734e+03
       tablef(6)= -1.64236e+04
       tablef(7)=  1.03611e+04
       
    case ("TRP")
       
       tablef(1)=  1.57239e+01
       tablef(2)= -4.86337e+00
       tablef(3)=  2.53387e+02
       tablef(4)= -3.27155e+03
       tablef(5)=  1.36615e+04
       tablef(6)= -2.38064e+04
       tablef(7)=  1.50413e+04
       
    case ("PRO")
       
       tablef(1)=  8.65194e+00
       tablef(2)= -1.94281e+00
       tablef(3)= -8.67886e+00
       tablef(4)=  6.77632e+02
       tablef(5)= -3.75813e+03
       tablef(6)=  7.51809e+03
       tablef(7)= -5.13230e+03
       
    case ("SER")
       
       tablef(1)=  1.40636e+01
       tablef(2)= -7.77211e+00
       tablef(3)=  1.51783e+02
       tablef(4)= -1.08221e+03
       tablef(5)=  3.65908e+03
       tablef(6)= -5.78759e+03
       tablef(7)=  3.42510e+03

    case ("THR")
       
       tablef(1)=  1.31250e+01
       tablef(2)= -5.68283e+00
       tablef(3)=  1.83111e+02
       tablef(4)= -1.78792e+03
       tablef(5)=  7.19676e+03
       tablef(6)= -1.26408e+04
       tablef(7)=  8.08658e+03

    case ("CYS")
       
       tablef(1)=  1.88068e+01
       tablef(2)= -7.68263e+00
       tablef(3)=  2.61934e+02
       tablef(4)= -2.71683e+03
       tablef(5)=  1.05382e+04
       tablef(6)= -1.77110e+04
       tablef(7)=  1.09546e+04
       
    case ("TYR")
       
       tablef(1)=  1.42182e+01
       tablef(2)= -8.50753e+00
       tablef(3)=  1.41348e+02
       tablef(4)= -1.63382e+03
       tablef(5)=  6.22461e+03
       tablef(6)= -9.83758e+03
       tablef(7)=  5.59816e+03
       
    case ("ASN")
       
       tablef(1)=  1.99148e+01
       tablef(2)=  7.28892e+00
       tablef(3)= -1.56404e+02
       tablef(4)=  1.21412e+03
       tablef(5)= -4.11987e+03
       tablef(6)=  6.00434e+03 
       tablef(7)= -3.13654e+03
       
    case ("GLN")
       
       tablef(1)=  1.90477e+01
       tablef(2)= -1.47112e+00
       tablef(3)= -2.26318e+00
       tablef(4)=  2.30599e+01
       tablef(5)= -4.43784e+02
       tablef(6)=  1.10422e+03
       tablef(7)= -7.93490e+02
       
    case ("ASP")
       
       tablef(1)=  2.02363e+01
       tablef(2)= -7.72888e+00
       tablef(3)=  1.51911e+02
       tablef(4)= -8.10751e+02
       tablef(5)=  2.25145e+03
       tablef(6)= -3.87827e+03
       tablef(7)=  2.86620e+03

    case ("GLU")
       
       tablef(1)=  1.92537e+01
       tablef(2)= -6.19618e+00
       tablef(3)=  3.58322e+01
       tablef(4)=  6.03373e+01
       tablef(5)= -9.50129e+02
       tablef(6)=  1.74192e+03
       tablef(7)= -8.43541e+02

    case ("LYS")
       
       tablef(1)=  1.09498e+01
       tablef(2)=  4.90064e+00
       tablef(3)= -1.03591e+02
       tablef(4)=  9.05869e+02
       tablef(5)= -3.01518e+03
       tablef(6)=  4.13753e+03
       tablef(7)= -1.95128e+03
       
    case ("ARG")
       
       tablef(1)=  2.33617e+01
       tablef(2)= -4.66830e+00
       tablef(3)=  6.09538e+01
       tablef(4)= -7.94534e+02
       tablef(5)=  2.32886e+03
       tablef(6)= -2.85862e+03
       tablef(7)=  1.23861e+03
       
    case ("HIS")
       
       tablef(1)=  2.13928e+01
       tablef(2)= -4.11645e+00
       tablef(3)=  4.89328e+01
       tablef(4)= -3.14819e+02
       tablef(5)=  9.54856e+02
       tablef(6)= -1.99606e+03
       tablef(7)=  1.72680e+03
       
    case ("SOL")
       
       tablef(1)= .04*(9.99891e+00)
       tablef(2)= .04*(7.18475e-03)
       tablef(3)= .04*(-1.1436e+00)
       tablef(4)= .04*(1.31469e-01)

    case default
       
       write(10,"(A)") "Something is wrong! Check the residue names. I've got ", resname
       stop "Something is wrong! Check the residue names."
    end Select
 
  end Subroutine ff_table
  
!     ################## FORM FACTOR #####################
  real(dp) function ff(a,b,c,d,f,g,h,q)
    implicit none
    real(dp), intent(in):: a,b,c,d,f,g,h,q
    
    ff=0.d0
    ff=a+(b*q)+(c*q*q)+(d*q*q*q)+(f*q*q*q*q)+(g*q*q*q*q*q)+(h*q*q*q*q*q*q)

  end function ff
  
!     ################## GYRATION RADIUS ###################
  real(dp) function RgwI(qa,qb,q,I)
    implicit none
    real(dp), intent(in) :: qa,qb,q(1:nlines_exp),I(1:nlines_exp)
    real(dp) :: q_guinier(150),I_guinier(150),slope,I_0
    real(dp) :: sumI,sumq,num,den,avgQ,avgI,arg1,arg2
    integer  :: k,cont
    
    q_guinier(:)=0.d0
    I_guinier(:)=0.d0
    sumq=0
    sumI=0
    num=0
    den=0
    
    k=1
    cont=0
    do while ((cont .le. 100) .and. (q(k) .le. qb))
       k=k+1
       if (q(k) .le. qa) cycle
       cont=cont+1
       q_guinier(cont)=q(k)*q(k)
       I_guinier(cont)=log(I(k))
       sumI=sumI+I_guinier(cont)
       sumq=sumq+q_guinier(cont)
    end do
    
    avgQ=sumq/cont
    avgI=sumI/cont
    do k=1,cont
       arg1= ((q_guinier(k)-avgQ)*(I_guinier(k)-avgI))
       arg2= ((q_guinier(k)-avgQ)*(q_guinier(k)-avgQ))
       num = num + arg1
       den = den + arg2
    end do
    
    slope=num/den !ln(I)/q
    I_0=exp(avgI-(slope*avgQ))
    RgwI = sqrt(-3.0*slope)
    
  end function RgwI
  
end Module saxs_theo
