module variables
  implicit none
  integer,parameter       :: dp = selected_real_kind(15, 307)
  real(dp),parameter      :: kb = 0.00138, pi = 4.d0*datan(1.d0), angle1=pi/180
  character(60),parameter :: pdb_format = '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,11X,A1,A2)'
  character(40),parameter :: exp_format = '(3(E11.4,1X))'
  
   type pdb
     character(6) :: dummy
     integer      :: serialnumber
     character(4) :: name
     character    :: alternate
     character(3) :: residue
     character    :: chainID
     integer      :: ressqnumber     
     character    :: codeinsert
     real(dp)     :: x,y,z,occupance,tempfactor
     character(2) :: element,charge
  end type pdb
  type esp
     real(dp)     :: q, intense,error
  end type esp
  type pivot
     character(1) :: ss_type
     integer      :: residue
  end type pivot
  
  type(pdb),allocatable   :: atom(:), CA(:),N(:), water(:)
  type(esp),allocatable   :: exp_data(:), theo_data(:),log_exp(:)
  type(pivot),allocatable :: list(:)
  integer,allocatable     :: flag_CA(:),flag_N(:)
  integer,allocatable     :: flag_C(:),flag_CB(:),flag_last(:)
  real(dp),allocatable    :: factor(:,:)

  integer       :: nlines_pdb,nlines_exp,nlines_pivot,k,seed,n_open=1000
  integer       :: n_passos,temp_min=273,temp_max=373,minus=10,temp
  integer       :: bin,first=0,res_count
  real(dp)      :: eta_chi,eta_rg,eta_bond,eta_lj,turn=1,tun
  real(dp)      :: chi=0,E_rg,E_hbond,E_lj
  character(40) :: PDB_input,EXP_input

  
end module variables
