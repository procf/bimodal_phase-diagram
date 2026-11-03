! Fortran module for analyzing a BD colloid sim GSD file
! NOTE: requires compilation with compile-module
! NOTE: is run from a matching sim-analysis-bi-BD.py Python script
! NOTE: ASSUMES 2 colloid types 
! (Deepak Mangal and Rob Campbell)

! Calculates:
! * coordination number distribution (Z counts) for each frame
! * mean squared displacement (MSD) of colloids for each frame
! * radial distribution function g(r) for colloid-colloid distribution in one frame
! * particle counts to approximate the pair correlation function h(r1,r2) 
!     for the colloid-colloid distribution in one frame
! * static structure factor S(q) (and associated g(r)) for colloids in the final frame
!     or a selection of frames
! * edgelist creation (for network analysis) for all colloid-colloid connections in each frame
! * corrected kinetic temperature and pressure (removing solvent flow velocity) in one frame
!   NOTE: only to be run on DPD data with solvents
! * velocity profile (for shear analysis) for selected particles in one frame
! * fabric tensor for colloids in one frame
! * pore size calculation for the final frame or a selection of frames, using two methods: 
!     Torquato's Pore Size Distribution and Gubbibn's Pore Size Distribution  
!     ! NOTE: This is the most complicated code, it has MULTIPLE subroutines
!     !       uses a linked-list implementation, and REQUIRES an additional module
!     !       generated with solvopt.f90 (which is also compiled in compile module)

  ! NOTE: f2py only understands limited number of KIND parameters
  ! use real(kind=8) not real64

  ! NOTE: fortran indexes from 1 by default, switch to base 0 to match C/Python
  ! for indexing from 0, 0:ncolloid-1 is of length ncolloid, 0:2 is length 3

  ! NOTE: must tell f2py about any variable dependencies or you will get Value Errors


!###########################################
subroutine coordination_number(nframes,Lbox,ncolloids,typeid,&
R_C,m_xys,allpos_allframe,cut_off,npairs,n_allcombos,pairs,Zs_array)
! counts the number of contacts (Z), AKA cooordination number, of each colloid in each frame

  implicit none

  !!! INPUTS
  ! the number of frames 
  integer, intent(in) :: nframes
  ! the simulation box size (L_X, L_Y, L_Z)
  real(8), intent(in) :: Lbox(0:2)

  ! total number of colloids 
  integer, intent(in) :: ncolloids
  ! the typeid of all particles
  integer, intent(in) :: typeid(0:ncolloids-1)
  ! the radius of each colloid particle
  real(8), intent(in) :: R_C(0:ncolloids-1)
  ! xyz position of each colloid particle
  real(8), intent(in) :: allpos_allframe(0:nframes-1,0:ncolloids-1,0:2)
  ! all the pair combinations (ex: 1-1, 1-2, ...) : npairs, pairs
  ! and with the any-combos (ex: any-any, 1-any, 2-any, ...) : n_allcombos
  integer, intent(in) :: npairs, n_allcombos
  integer, intent(in) :: pairs(0:npairs-1,0:1)
  ! the cut-off distance for defining a bond
  real(8), intent(in) :: cut_off
  ! the xy tilt factor per frame (box deformation from shear in x-direction)
  real(8), intent(in) :: m_xys(0:nframes-1)

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(ncolloids) typeid,R_C
  !f2py depend(nframes,ncolloids) allpos_allframe
  !f2py depend(npairs) pairs
  !f2py depend(nframes) m_xys


  !!! OUTPUTS
  ! contact number counts per particle per frame
  integer, intent(out) :: Zs_array(0:nframes-1,0:ncolloids-1,0:n_allcombos-1)


  !!! INTERNAL VALUES
  ! the center-center interparticle distance (r_ij), 
  ! the surface-surface interparticle distance (h_ij), 
  ! 1/box-size, and the modification for particles
  ! interacting across a sheared/tilted box boundary (img)   
  real(8) :: r_ij(0:2), h_ij, inv_Lbox(0:2), img
  ! frame tag (f), particle tags (i, j)
  integer :: f, i, j
  ! m_xy and allpos for each frame
  real(8) :: m_xy
  real(8) :: allpos(0:ncolloids-1,0:2)
  ! reverse pair-combination counters (given pair i-j, find pair j-i)
  integer :: c(0:1)
  ! index of i-j pair (p) and j-i pair (d)
  integer :: p, d


  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  Zs_array = 0

  !!! calculate Z for each colloid particle in every frame
  do f=0,(nframes-1)

    ! reset frame-specific variables and set new values
    allpos = 0
    allpos = allpos_allframe(f,:,:)
    m_xy = 0
    m_xy = m_xys(f)

    ! for all colloid particles
    do i=0,(ncolloids-2)
      ! check the interaction with each other colloid particle
      do j=i+1,(ncolloids-1)
        ! calculate the center-center interparticle distance
        r_ij = allpos(i,:)-allpos(j,:)
        ! adjust for periodic boundaties in all directions (x=0, y=1, z=2)
        ! including possible box tilt from shear in the x-direction (m_xy)
        r_ij(2) = r_ij(2) - Lbox(2) * dnint(r_ij(2)*inv_Lbox(2))
        img = Lbox(1) * dnint(r_ij(1)*inv_Lbox(1))
        r_ij(1) = r_ij(1) - img
        r_ij(0) = r_ij(0) - img * m_xy
        r_ij(0) = r_ij(0) - Lbox(0) * dnint(r_ij(0)*inv_Lbox(0))
        ! convert to surface-surface distance
        h_ij = dsqrt(sum(r_ij*r_ij)) - R_C(i) - R_C(j)

        ! check for bonds and count once for each colloid
        if(h_ij <= cut_off) then

          ! any colloid with any other colloid:
          Zs_array(f,i,0) = Zs_array(f,i,0) + 1
          Zs_array(f,j,0) = Zs_array(f,j,0) + 1

          ! for all other combinations: 
          do p=0,(npairs-1)
            ! find the a-b type pair that describes this interaction
            ! (pairs start at 1, for BD subtract 1 before passing pairs, to match typeid index from 0)
            if((typeid(i) .eq. pairs(p,0)) .and. (typeid(j) .eq. pairs(p,1))) then
              ! find the matching b-a pair
              c(0) = pairs(p,1) ! typeid of j 
              c(1) = pairs(p,0) ! typeid of i
              do d=0,(npairs-1)
                if ((c(0) .eq. pairs(d,0)) .and. (c(1) .eq. pairs(d,1))) then
                  ! increase both counts by 1
                  Zs_array(f,i,p+1) = Zs_array(f,i,p+1) + 1
                  Zs_array(f,j,d+1) = Zs_array(f,j,d+1) + 1
                  ! and Z_a-any and Z_b-any
                  Zs_array(f,i,npairs+c(1)) = Zs_array(f,i,npairs+c(1)) + 1
                  Zs_array(f,j,npairs+c(0)) = Zs_array(f,j,npairs+c(0)) + 1
                endif
              enddo
            endif
          enddo

        endif

      enddo
    enddo

  enddo

end subroutine coordination_number



!###########################################
subroutine msd_calculation(data_outpath,nframes,ncolloids,ncolloid1,ncolloid2,&
allcollpos,allcollpos1,allcollpos2)
! for all frames, calculate the mean-squared displacement of all colloid particles

  implicit none

  !!! INPUTS
  ! the filepath to the output folder 
  character*256, intent(in) :: data_outpath
  ! number of frames, total number of each type of particle
  integer, intent(in) :: nframes, ncolloids, ncolloid1, ncolloid2
  ! xyz positions of all colloids in all frames
  real(8), intent(in) :: allcollpos(0:nframes-1,0:ncolloids-1,0:2) 
  ! xyz positions of all type 1 colloids in all frames
  real(8), intent(in) :: allcollpos1(0:nframes-1,0:ncolloid1-1,0:2)
  ! xyz positions of all type 2 colloids in all frames
  real(8), intent(in) :: allcollpos2(0:nframes-1,0:ncolloid2-1,0:2)

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nframes,ncolloids) allcollpos
  !f2py depend(nframes,ncolloid1) allcollpos1
  !f2py depend(nframes,ncolloid2) allcollpos2

  !!! OUTPUTS
  ! N/A (output file created)

  !!! INTERNAL VALUES
  ! the center-center change in position (dcollpos) of each colloid
  ! and it's value squared (dcollpos_2)
  real(8) :: dcollpos(0:ncolloids-1,0:2), dcollpos_2(0:ncolloids-1)
  ! the center-center change in position of each type 1 colloid
  ! (dpos1), and it's value squared (dpos1_2)
  real(8) :: dcollpos1(0:ncolloid1-1,0:2), dcollpos1_2(0:ncolloid1-1)
  ! the center-center change in position of each type 2 colloid
  ! (dpos2), and it's value squared (dpos2_2)
  real(8) :: dcollpos2(0:ncolloid2-1,0:2), dcollpos2_2(0:ncolloid2-1)
  ! frame tag (i)
  integer :: i
  ! colloid mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: coll_msd_value, coll_msd_corrsmplvar, coll_msd_smplstddev
  ! type 1 mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: coll_msd1_value, coll_msd1_corrsmplvar, coll_msd1_smplstddev
  ! type 2 mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: coll_msd2_value, coll_msd2_corrsmplvar, coll_msd2_smplstddev
  ! output filename/path
  character(len=256) :: filename_msd 

  !!! SET INITIAL VALUES
  ! N/A

  !!! create a file to record the MSD and standard deviation
  filename_msd = trim(adjustl(data_outpath)) // '/msd.txt'
  open(unit=14,file=trim(filename_msd),action='write')
  ! label the output data values
  write(14,'(A,x,A,x,A,x,A,x,A,x,A,x,A)') 'frame','msd_colloid','msd_coll_smplstddev',&
'msd_colloid_1','msdc1_smplstddev','msd_colloid_2','msdc2_smplstddev'

  !!! calculate the msd for each frame
  do i=1,nframes-1
    !! for COLLOIDS
    ! get the change in position relative to frame 0 for all colloids
    dcollpos = allcollpos(i,:,:) - allcollpos(0,:,:)
    ! get the square of the change in position
    !   dcollpos = (particles, positions (x,y,z))
    !   dcollpos*dcollpos = (particles, squared positions (x*x, y*y, z*z)
    !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
    !   -> creates a 1D array of dcollpos_2 for all particles
    dcollpos_2 = sum(dcollpos*dcollpos, dim=2) 
    ! sum all the squares of relative change in position, and take the mean
    coll_msd_value = sum(dcollpos_2)/ncolloids
    ! calculate the sample standard deviation (1/N-1) of the MSD
    coll_msd_corrsmplvar = (sum(dcollpos_2*dcollpos_2)/dble(ncolloids) - coll_msd_value**2) / dble(ncolloids-1) 
    coll_msd_smplstddev = dsqrt(coll_msd_corrsmplvar)

    ! for type 1
    ! get the change in position relative to frame 0 for all type 1 colloids
    dcollpos1 = allcollpos1(i,:,:) - allcollpos1(0,:,:)
    ! get the square of the change in position
    !   dcollpos1 = (particles, positions (x,y,z))
    !   dcollpos1*dcollpos1 = (particles, squared positions (x*x, y*y, z*z)
    !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
    !   -> creates a 1D array of dcollpos1_2 for all particles
    dcollpos1_2 = sum(dcollpos1*dcollpos1, dim=2)
    ! sum all the squares of relative change in position, and take the mean
    coll_msd1_value = sum(dcollpos1_2)/ncolloid1
    ! calculate the sample standard deviation (1/N-1) of the MSD
    coll_msd1_corrsmplvar = (sum(dcollpos1_2*dcollpos1_2)/dble(ncolloid1) - &
coll_msd1_value**2) / dble(ncolloid1-1)
    coll_msd1_smplstddev = dsqrt(coll_msd1_corrsmplvar)

    ! for type 2
    ! get the change in position relative to frame 0 for all type 2 colloids
    dcollpos2 = allcollpos2(i,:,:) - allcollpos2(0,:,:)
    ! get the square of the change in position
    !   dcollpos2 = (particles, positions (x,y,z))
    !   dcollpos2*dcollpos2 = (particles, squared positions (x*x, y*y, z*z)
    !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
    !   -> creates a 1D array of dcollpos2_2 for all particles
    dcollpos2_2 = sum(dcollpos2*dcollpos2, dim=2)
    ! sum all the squares of relative change in position, and take the mean
    coll_msd2_value = sum(dcollpos2_2)/ncolloid2
    ! calculate the sample standard deviation (1/N-1) of the MSD
    coll_msd2_corrsmplvar = (sum(dcollpos2_2*dcollpos2_2)/dble(ncolloid2) - &
coll_msd2_value**2) / dble(ncolloid2-1)
    coll_msd2_smplstddev = dsqrt(coll_msd2_corrsmplvar)
      
    ! write data to the output file
    write(14,'(I5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5)') i, &
coll_msd_value, coll_msd_smplstddev, &
coll_msd1_value, coll_msd1_smplstddev, coll_msd2_value, coll_msd2_smplstddev
  enddo

  ! close output file
  close(14)

end subroutine msd_calculation


!###########################################
subroutine gofr_calc(data_outpath,Lbox,m_xy,nparticle,ncolloids,ncolloidA,ncolloidB,typeid,&
colloidA_typeid,colloidB_typeid,pos,bin_width,nlayers)
! calculate the radial distribution function g(r), the probability density per unit volume,
! for colloid-colloid interactions for a single frame

  implicit none

  !!! INPUTS
  ! the filepath to the output folder 
  character*256, intent(in) :: data_outpath
  ! the simulation box size (L_X, L_Y, L_Z), and the
  ! xy tilt factor (box deformation from shear in x-direction)
  real(8), intent(in) :: Lbox(0:2), m_xy
  ! total number of particles 
  integer, intent(in) :: nparticle
  ! total number of colloids
  integer, intent(in) :: ncolloids, ncolloidA, ncolloidB
  ! the typeid of all particles
  integer, intent(in) :: typeid(0:nparticle-1)
  ! colloid typeids
  integer, intent(in) :: colloidA_typeid, colloidB_typeid
  ! the xyz position of all particles
  real(8), intent(in) :: pos(0:nparticle-1,0:2)
  ! the size of each bin/layer/slab in the search box
  real(8), intent(in) :: bin_width
  ! the number of bins/layers/slabs in the search box
  integer, intent(in) :: nlayers

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nparticle) typeid,pos 

  !!! OUTPUTS
  ! N/A (output files created)

  !!! INTERNAL VALUES
  !! For counting particles in each bin !!
  ! 1/bin-size
  real(8) :: inv_bin_width 
  ! the center-center interparticle distance (dr) and 
  ! corrected center-center interparticle distance (r1),
  ! 1/box-size, and the modification for particles interacting
  ! across a sheared/tilted box boundary (img)
  real(8) :: dr(0:2), r1, inv_Lbox(0:2), img
  ! particle tags (i, j) and bin tag (k)
  integer :: i, j, k
  ! lists for counting the number of other particles at each bin
  ! hist_bin_all = any colloid-colloid
  integer :: hist_bin_all(0:nlayers-1)
  ! colloid type dependent: (a = colloid 1, b = colloid 2)
  integer :: hist_bin_aa(0:nlayers-1), hist_bin_bb(0:nlayers-1)
  integer :: hist_bin_ab(0:nlayers-1)

  !! For the final calculation !!
  ! the number pi, the current and previous volumes
  real(8) :: pi, curr_volume, prev_volume
  ! temporary shell_volume for each bin
  real(8) :: shell_volume
  ! total box volume
  real(8) :: Lbox_vol
  ! number density (rho)
  real(8) :: rho_all, rho_aa, rho_bb, rho_ab
  ! ideal gas AKA Gaussian distribution for each combo
  real(8) :: ideal_aa_in_bin
  real(8) :: ideal_bb_in_bin
  real(8) :: ideal_all_in_bin
  real(8) :: ideal_ab_in_bin

  ! g(r), the particle density calculated for a bin/layer
  real(8) :: all_density  
  real(8) :: aa_density, bb_density, ab_density
  ! output filename/path
  character(len=256) :: filename_gofr 

  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  ! calculate 1/bin-size
  inv_bin_width = 1.d0/bin_width
  ! calculate pi (pi = arctan(1) * 4)
  pi = 4.d0*datan(1.d0)
  ! calculate the total volume
  Lbox_vol = Lbox(0)*Lbox(1)*Lbox(2)
  ! calculate the number density of colloids, rho
  rho_aa  = (ncolloidA / Lbox_vol)
  rho_bb =  (ncolloidB / Lbox_vol)
  rho_ab = (sqrt(real(ncolloidA*ncolloidB)) / Lbox_vol) ! exact
  !rho_ab = ((ncolloidA+ncolloidB) / (2*Lbox_vol)) ! approx
  rho_all = (ncolloids / Lbox_vol)

  ! initialize history bins at 0
  hist_bin_all = 0
  hist_bin_aa = 0
  hist_bin_bb = 0
  hist_bin_ab = 0

  !!! count the number of particles in each bin/layer
  ! for N-1 particles
  do i=0,nparticle-2
    ! interacting with the remaining unbinned particles
    do j=i+1,nparticle-1

      ! for all interactions (COLLOID-COLLOID)
      ! calculate the center-center interparticle distance
      dr = pos(i,:) - pos(j,:)
      ! adjust dr for particles interacting across a boundary:
      ! z-boundary
      dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
      ! y-boundary
      img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
      dr(1) = dr(1) - img
      ! x-boundary
      ! (including accounting for a sheared/tilted frame)
      dr(0) = dr(0) - img * m_xy
      dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
      ! set the corrected center-center interparticle distance
      r1 = dsqrt(sum(dr*dr))
      ! find which bin/layer the interparticle distance is
      ! located in (bins range 0:nlayers-1)
      k = floor(r1*inv_bin_width)
      ! if the interaction is inside the search box
      if(k < nlayers) then
        ! increase the number of particles interactions
        ! in that bin by 2 (once for each particle)
        hist_bin_all(k) = hist_bin_all(k) + 2
        ! if TYPE1-TYPE1
        if((typeid(i) .eq. colloidA_typeid) .AND. (typeid(j) .eq. colloidA_typeid)) then
          hist_bin_aa(k) = hist_bin_aa(k) + 2
        ! if TYPE2-TYPE2
        else if((typeid(i) .eq. colloidB_typeid) .AND. (typeid(j) .eq. colloidB_typeid)) then
          hist_bin_bb(k) = hist_bin_bb(k) + 2
        ! if TYPE1-TYPE2
        else if(typeid(i) .ne. typeid(j)) then
          hist_bin_ab(k) = hist_bin_ab(k) + 2
        endif
      endif

    enddo
  enddo
                       
  ! create file to record gofr data
  filename_gofr = trim(adjustl(data_outpath)) // '/gofr.txt'
  open(unit=15,file=trim(filename_gofr),action='write')
  write(15,'(A,x,A,x,A,x,A,x,A)') 'position', 'gofr_cc', &
'gofr_c1c1', 'gofrc1c2', 'gofr_c2c2'

  ! for each bin/layer
  do i=0,nlayers-1
    ! calculate the volume of the current layer
    curr_volume = (4.d0/3.d0) * pi * ((i+1)*bin_width)**3
    prev_volume = (4.d0/3.d0) * pi * (i*bin_width)**3
    shell_volume = (curr_volume - prev_volume)

    ! type 1-1
    ideal_aa_in_bin = ncolloidA * rho_aa * shell_volume
    aa_density = hist_bin_aa(i) / ideal_aa_in_bin

    ! type 2-2
    ideal_bb_in_bin = ncolloidB * rho_bb * shell_volume
    bb_density = hist_bin_bb(i) / ideal_bb_in_bin

    ! type 1-2 / 2-1 (times 2 for counting all AB/BA, i.e +2)
    ideal_ab_in_bin = 2 * ncolloidA * (ncolloidB/Lbox_vol) * shell_volume
    ab_density = hist_bin_ab(i) / ideal_ab_in_bin

    ! type any-any
    ideal_all_in_bin = ncolloids * rho_all * shell_volume
    all_density = hist_bin_all(i) / ideal_all_in_bin


    ! write the data
    write(15,'(f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5)') (i+0.5d0)*bin_width, &
all_density, aa_density, bb_density, ab_density

  enddo

  ! close files
  close(15)

end subroutine gofr_calc


!###########################################
subroutine pcf_calc(Lbox,m_xy,ncolloids,pos,bin_size,dims_max,array_dim,&
slice_width,typeid,colloid1_typeid,colloid2_typeid,&
counts_XY_CC,counts_XZ_CC,counts_YZ_CC,counts_XY_C1C1,counts_XZ_C1C1,counts_YZ_C1C1,&
counts_XY_C1C2,counts_XZ_C1C2,counts_YZ_C1C2,counts_XY_C2C2,counts_XZ_C2C2,counts_YZ_C2C2)
! collect data about the particle location in 2D slices of the system to approximate 
! the pair correlation function h(r1,r2), the probability density per unit volume
! of finding colloids at two locations, r1 and r2, for one frame

  implicit none

  !!! INPUTS
  ! the simulation box size (L_X, L_Y, L_Z), and the
  ! xy tilt factor (box deformation from shear in x-direction)
  real(8), intent(in) :: Lbox(0:2), m_xy
  ! total number of colloids
  integer, intent(in) :: ncolloids
  ! the xyz position of all particles
  real(8), intent(in) :: pos(0:ncolloids-1,0:2)
  ! the side-length of each bin in the search box
  real(8), intent(in) :: bin_size
  ! the maximum size of the the search box in one dimension
  real(8), intent(in) :: dims_max 
  ! the dimensions of the output arrays
  integer, intent(in) :: array_dim
  ! the width of each 2D slice
  real(8), intent(in) :: slice_width
  ! the typeid of all colloids
  integer, intent(in) :: typeid(0:ncolloids-1)
  ! colloid typeids
  integer, intent(in) :: colloid1_typeid, colloid2_typeid

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(ncolloids) typeid,pos

  !!! OUTPUTS 
  ! lists for counting the number of colloids at each bin
  ! in the xy-plane, xz-plane, and yz-plane
  integer, intent(out) :: counts_XY_CC(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XZ_CC(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_YZ_CC(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XY_C1C1(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XZ_C1C1(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_YZ_C1C1(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XY_C1C2(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XZ_C1C2(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_YZ_C1C2(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XY_C2C2(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_XZ_C2C2(0:array_dim-1,0:array_dim-1)
  integer, intent(out) :: counts_YZ_C2C2(0:array_dim-1,0:array_dim-1)

  !!! INTERNAL VALUES
  ! the center-center interparticle distance (dr),
  ! 1/box-size, and the modification for particles interacting
  ! across a sheared/tilted box boundary (img)
  real(8) :: dr(0:2), inv_Lbox(0:2), img
  ! particle tags (i, j) and bin tags (x,y,z)
  integer :: i, j, x, y, z


  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox

  ! initialize colloid-colloid output arrays at 0
  counts_XY_CC = 0
  counts_XZ_CC = 0
  counts_YZ_CC = 0
  counts_XY_C1C1 = 0
  counts_XZ_C1C1 = 0
  counts_YZ_C1C1 = 0
  counts_XY_C1C2 = 0
  counts_XZ_C1C2 = 0
  counts_YZ_C1C2 = 0
  counts_XY_C2C2 = 0
  counts_XZ_C2C2 = 0
  counts_YZ_C2C2 = 0

  !!! count the number of particles in each bin
  ! for N-1 particles
  do i=0,ncolloids-1
    ! interacting with the remaining unbinned particles
    do j=0,ncolloids-1

      ! if the particles are not the same
      if (i .ne. j) then

        ! calculate the center-center interparticle distance
        dr = pos(i,:) - pos(j,:)
        ! adjust dr for particles interacting across a boundary:
        ! z-boundary
        dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
        ! y-boundary
        img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
        dr(1) = dr(1) - img
        ! x-boundary
        ! (including accounting for a sheared/tilted frame)
        dr(0) = dr(0) - img * m_xy
        dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
  
        ! find the bins in the XY Plane
        ! if the z-distance is within the 2D slice
        if (abs(dr(2)) <= slice_width) then
          ! if the interaction is inside the search box
          if((abs(dr(0)) < dims_max) .AND. (abs(dr(1)) < dims_max)) then
            ! find which XY bin this interparticle distance is in
            x = floor((dr(0) + dims_max)/bin_size)
            y = floor((dr(1) + dims_max)/bin_size)
            ! increase the number of particles in that bin by 1 
            counts_XY_CC(x,y) = counts_XY_CC(x,y) + 1
            ! if TYPE1-TYPE1
            if((typeid(i) .eq. colloid1_typeid) .AND. (typeid(j) .eq. colloid1_typeid)) then
              counts_XY_C1C1(x,y) = counts_XY_C1C1(x,y) + 1
            ! if TYPE2-TYPE2
            else if((typeid(i) .eq. colloid2_typeid) .AND. (typeid(j) .eq. colloid2_typeid)) then
              counts_XY_C2C2(x,y) = counts_XY_C2C2(x,y) + 1
            ! if TYPE1-TYPE2
            else if(typeid(i) .ne. typeid(j)) then
              counts_XY_C1C2(x,y) = counts_XY_C1C2(x,y) + 1
            endif 
          endif
        endif

        ! find the bins in the XZ Plane
        ! if the y-distance is within the 2D slice
        if (abs(dr(1)) <= slice_width) then
          ! if the interaction is inside the search box
          if((abs(dr(0)) < dims_max) .AND. (abs(dr(2)) < dims_max)) then
            ! find which XZ bin this interparticle distance is in
            x = floor((dr(0) + dims_max)/bin_size)
            z = floor((dr(2) + dims_max)/bin_size)
            ! increase the number of particles in that bin by 1 
            counts_XZ_CC(x,z) = counts_XZ_CC(x,z) + 1
            ! if TYPE1-TYPE1
            if((typeid(i) .eq. colloid1_typeid) .AND. (typeid(j) .eq. colloid1_typeid)) then
              counts_XZ_C1C1(x,z) = counts_XZ_C1C1(x,z) + 1
            ! if TYPE2-TYPE2
            else if((typeid(i) .eq. colloid2_typeid) .AND. (typeid(j) .eq. colloid2_typeid)) then
              counts_XZ_C2C2(x,z) = counts_XZ_C2C2(x,z) + 1
            ! if TYPE1-TYPE2
            else if(typeid(i) .ne. typeid(j)) then
              counts_XZ_C1C2(x,z) = counts_XZ_C1C2(x,z) + 1
            endif 
          endif 
        endif 

        ! find the bins in the YZ Plane
        ! if the x-distance is within the 2D slice
        if (abs(dr(0)) <= slice_width) then
          ! if the interaction is inside the search box
          if((abs(dr(1)) < dims_max) .AND. (abs(dr(2)) < dims_max)) then
            ! find which YZ bin this interparticle distance is in
            y = floor((dr(1) + dims_max)/bin_size)
            z = floor((dr(2) + dims_max)/bin_size)
            ! increase the number of particles in that bin by 1 
            counts_YZ_CC(y,z) = counts_YZ_CC(y,z) + 1
            ! if TYPE1-TYPE1
            if((typeid(i) .eq. colloid1_typeid) .AND. (typeid(j) .eq. colloid1_typeid)) then
              counts_YZ_C1C1(y,z) = counts_YZ_C1C1(y,z) + 1
            ! if TYPE2-TYPE2
            else if((typeid(i) .eq. colloid2_typeid) .AND. (typeid(j) .eq. colloid2_typeid)) then
              counts_YZ_C2C2(y,z) = counts_YZ_C2C2(y,z) + 1
            ! if TYPE1-TYPE2
            else if(typeid(i) .ne. typeid(j)) then
              counts_YZ_C1C2(y,z) = counts_YZ_C1C2(y,z) + 1
            endif 
          endif 
        endif 
 
      endif

    enddo
  enddo

end subroutine pcf_calc


!###########################################
subroutine structure_factor(data_outpath,pop_type,nframes,framechoice,Lbox,m_xys,&
ncolloids,ncolloidA,ncolloidB,typeid,colloidA_typeid,colloidB_typeid,pos,rmax,rmin,nq,q)
! calculate the structure factor S(q) -- the fluctuation of the 
! **Fourier component** of the fluctuations of pair density in space
! using the radial correlation function h(r) = g(r)-1
! (note: this also recalculates and saves g(r))

  implicit none

  !!! INPUTS
  ! the filepath to the output folder 
  character*256, intent(in) :: data_outpath
  ! tag for population type
  character*256, intent(in) :: pop_type
  ! the number of frames and the selected frames
  integer, intent(in) :: nframes
  integer, intent(in) :: framechoice(0:nframes-1)
  ! the simulation box size (L_X, L_Y, L_Z)
  real(8), intent(in) :: Lbox(0:2)
  ! the xy tilt factor per frame (box deformation from shear in x-direction)
  real(8), intent(in) :: m_xys(0:nframes-1) 
  ! total number of colloids
  integer, intent(in) :: ncolloids, ncolloidA, ncolloidB
  ! the typeid of all particles
  integer, intent(in) :: typeid(0:ncolloids-1)
  ! colloid typeids
  integer, intent(in) :: colloidA_typeid, colloidB_typeid
  ! the xyz position of all colloids
  real(8), intent(in) :: pos(0:nframes-1,0:ncolloids-1,0:2)

  !! for g(r)
  ! the max distance used in g(r)
  real(8), intent(in) :: rmax
  ! the min distance used in g(r) (usually particle size)
  real(8), intent(in) :: rmin

  !! for S(q)
  ! the number q values used to calculate S(q)
  integer, intent(in) :: nq
  ! the wavenumbers being used to calculate S(q)
  real(8), intent(in) :: q(0:nq-1)
 
  !f2py depend(ncolloids) typeid
  !f2py depend(nframes,ncolloids) pos
  !f2py depend(nframes) framechoice, m_xys
  !f2py depend(nq) q


  !!! INTERNAL VALUES
  ! 1/Lbox
  real(8) :: inv_Lbox(0:2) 
  ! image for adjusting position for shear flow box deformation
  real(8) :: img
  ! frame tag (f), particle tags (i, j) and bin tag (k)
  integer :: f, i, j, k
  
  !! g(r)
  ! the number of bins/layers/slabs in the search box for g(r)
  integer :: nlayers_g
  ! size of each bin/layer/slab (and 1/bin_width) for g(r)
  real(8) :: bin_width_g, inv_bin_width_g
  ! interparticle distance 
  real(8) :: dr(0:2) 
  ! corrected center-center interparticle distance
  real(8) :: r1  
  ! list for counting the number of other colloids at each bin
  integer,allocatable :: hist_bin_ab(:,:)
  integer,allocatable :: hist_bin_aa(:,:)
  integer,allocatable :: hist_bin_bb(:,:)
  integer,allocatable :: hist_bin_all(:,:) 
  ! g(r), the colloid density calculated for a bin/layer 
  real(8),allocatable :: ab_density(:,:)
  real(8),allocatable :: aa_density(:,:)
  real(8),allocatable :: bb_density(:,:)
  real(8),allocatable :: all_density(:,:)  

  ! the number pi, the current and previous volumes
  real(8) :: pi, curr_volume, prev_volume
  ! temporary shell_volume for each bin
  real(8) :: shell_volume
  ! total box volume
  real(8) :: Lbox_vol
  ! number density (rho)
  real(8) :: rho_all, rho_aa, rho_bb, rho_ab
  ! ideal gas AKA Gaussian distribution for each combo
  real(8) :: ideal_aa_in_bin
  real(8) :: ideal_bb_in_bin
  real(8) :: ideal_all_in_bin
  real(8) :: ideal_ab_in_bin

  ! output filename/path for g(r) data
  character(len=256) :: filename_gofr_sq 

  !! S(q) calculated from h(r) = g(r)-1 for all wavenumbers
  real(8) :: Sq_ab(0:nframes-1,0:nq-1)
  real(8) :: Sq_aa(0:nframes-1,0:nq-1)
  real(8) :: Sq_bb(0:nframes-1,0:nq-1)
  real(8) :: Sq_all(0:nframes-1,0:nq-1) 
  ! current position where S(q) is being calculated
  real(8) :: r 
  ! output filename/path for S(q) data
  character(len=256) :: filename_sofq 

  ! Declare ierr for status checks
  integer :: ierr


  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  ! calculate g(r) layers and bin_width
  bin_width_g = 0.001
  nlayers_g = int((rmax-rmin)/bin_width_g)
  bin_width_g = (rmax-rmin)/nlayers_g
  inv_bin_width_g = 1.d0/bin_width_g
  ! calculate pi (pi = arctan(1) * 4)
  pi = 4.d0*datan(1.d0)
  ! calculate the total volume
  Lbox_vol = Lbox(0)*Lbox(1)*Lbox(2)
  ! calculate the number density of colloids, rho
  rho_aa  = (ncolloidA / Lbox_vol)
  rho_bb =  (ncolloidB / Lbox_vol)
  rho_ab = (sqrt(real(ncolloidA*ncolloidB)) / Lbox_vol) ! exact
  !rho_ab = ((ncolloidA+ncolloidB) / (2*Lbox_vol)) ! approx
  rho_all = (ncolloids / Lbox_vol)

  ! dynamically allocate arrays for g(r) counts
  !allocate(hist_bin_cc_g(0:nframes-1,0:nlayers_g-1))
  !allocate(bin_cc_density_g(0:nframes-1,0:nlayers_g-1))
  ! Allocate arrays
  allocate(hist_bin_ab(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate hist_bin_ab"
    stop
  endif
  allocate(hist_bin_aa(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate hist_bin_aa"
    stop
  endif
  allocate(hist_bin_bb(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate hist_bin_bb"
    stop
  endif
  allocate(hist_bin_all(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate hist_bin_all"
    stop
  endif

  allocate(ab_density(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate ab_density"
    stop
  endif
  allocate(aa_density(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate aa_density"
    stop
  endif
  allocate(bb_density(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate bb_density"
    stop
  endif
  allocate(all_density(0:nframes-1,0:nlayers_g-1), stat=ierr)
  if (ierr /= 0) then
    print *, "Error: Unable to allocate all_density"
    stop
  endif

  ! initialize history bins and local bin density to 0 for g(r)
  hist_bin_ab = 0
  hist_bin_aa = 0
  hist_bin_bb = 0
  hist_bin_all = 0

  ab_density = 0.d0
  aa_density = 0.d0
  bb_density = 0.d0
  all_density = 0.d0

  ! initialize Sq values to 0
  Sq_ab = 0.d0
  Sq_aa = 0.d0
  Sq_bb = 0.d0
  Sq_all = 0.d0

  ! create file to record gofr data
  filename_gofr_sq = trim(adjustl(data_outpath)) // '/gofr_sq' // &
trim(adjustl(pop_type)) // '.csv'
  open(unit=14,file=trim(filename_gofr_sq),action='write')
  write(14,fmt="(*(g0:','))") "frame", "position","gofr","gofr_aa","gofr_ab","gofr_bb"

  ! create file to record S(q) data
  filename_sofq = trim(adjustl(data_outpath)) // '/sofq' // &
trim(adjustl(pop_type)) // '.csv'
  open(unit=15,file=trim(filename_sofq),action='write')
  write(15,fmt="(*(g0:','))") "frame", "q","Sq","Sq_aa","Sq_ab","Sq_bb"

  ! for each frame
  do f=0,(nframes-1)

    !!! calculate g(r)
    ! count the number of particles in each bin/layer
    ! for N-1 colloids
    do i=0,ncolloids-2
      ! interacting with the remaining unbinned colloids
      do j=i+1,ncolloids-1
        ! calculate the center-center interparticle distance
        dr = pos(f,i,:) - pos(f,j,:)
        ! adjust dr for particles interacting across a boundary:
        ! z-boundary
        dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
        ! y-boundary
        img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
        dr(1) = dr(1) - img
        ! x-boundary
        ! (including accounting for a sheared/tilted frame)
        dr(0) = dr(0) - img * m_xys(f)
        dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
        ! set the corrected center-center interparticle distance
        r1 = dsqrt(sum(dr*dr))-rmin
        ! find which bin/layer this interparticle distance is
        ! located in (bins range 0:nlayers-1)
        k = int(r1*inv_bin_width_g)
        ! if the interaction is inside the search box
        !if(k < nlayers_g) then
        if (k >= 0 .and. k < nlayers_g) then
          ! increase the number of particles interactions
          ! in that bin by 2 (once for each particle)
          hist_bin_all(f, k) = hist_bin_all(f, k) + 2
          ! if TYPE1-TYPE1
          if((typeid(i) .eq. colloidA_typeid) .AND. (typeid(j) .eq. colloidA_typeid)) then
            hist_bin_aa(f, k) = hist_bin_aa(f, k) + 2
          ! if TYPE2-TYPE2
          else if((typeid(i) .eq. colloidB_typeid) .AND. (typeid(j) .eq. colloidB_typeid)) then
            hist_bin_bb(f, k) = hist_bin_bb(f, k) + 2
          ! if TYPE1-TYPE2
          else if(typeid(i) .ne. typeid(j)) then
            hist_bin_ab(f, k) = hist_bin_ab(f, k) + 2
          endif
        endif
      enddo
    enddo

    ! for each bin/layer
    do i=0,nlayers_g-1
      ! calculate the volume for the current and previous bin/layer
      curr_volume = (4.d0/3.d0) * pi * (rmin+(i+1)*bin_width_g)**3
      prev_volume = (4.d0/3.d0) * pi * (rmin+i*bin_width_g)**3
      shell_volume = (curr_volume - prev_volume)

      ! type 1-1
      ideal_aa_in_bin = ncolloidA * rho_aa * shell_volume
      aa_density(f,i) = hist_bin_aa(f,i) / ideal_aa_in_bin

      ! type 2-2
      ideal_bb_in_bin = ncolloidB * rho_bb * shell_volume
      bb_density(f,i) = hist_bin_bb(f,i) / ideal_bb_in_bin

      ! type 1-2 / 2-1 (times 2 for counting all AB/BA, i.e +2)
      ideal_ab_in_bin = 2 * ncolloidA * (ncolloidB/Lbox_vol) * shell_volume
      ab_density(f,i) = hist_bin_ab(f,i) / ideal_ab_in_bin

      ! type any-any
      ideal_all_in_bin = ncolloids * rho_all * shell_volume
      all_density(f,i) = hist_bin_all(f,i) / ideal_all_in_bin

      ! write the data
      write(14,fmt="(*(g0:','))") framechoice(f), rmin+(i+0.5d0)*bin_width_g, all_density(f,i), &
aa_density(f,i), ab_density(f,i), bb_density(f,i)

    enddo


    !!! calculate S(q)
    ! for each wavenumber q
    do i=0,nq-1
      ! for each bin/layer used for g(r)
      do j=0,nlayers_g-1
        ! find the distance from 0 to the center of the current layer 
        ! (j previous bins*bin-width + 1/2 of the current bin-width)
        r = rmin+bin_width_g*(j+0.5)
        ! calculate the integral for S(q) using h(r)
        Sq_all(f,i) = Sq_all(f,i) + (r * (all_density(f,j)-1) * dsin(q(i)*r))
        Sq_aa(f,i) = Sq_aa(f,i) + (r * (aa_density(f,j)-1) * dsin(q(i)*r))
        Sq_ab(f,i) = Sq_ab(f,i) + (r * (ab_density(f,j)-1) * dsin(q(i)*r))
        Sq_bb(f,i) = Sq_bb(f,i) + (r * (bb_density(f,j)-1) * dsin(q(i)*r))
      enddo

      ! use the integral to finish calculating S(q)
      Sq_all(f,i) = 1.d0 + Sq_all(f,i) * 4.d0 * pi * rho_all * bin_width_g / q(i)
      Sq_aa(f,i) = 1.d0 + Sq_aa(f,i) * 4.d0 * pi * rho_aa * bin_width_g / q(i)
      Sq_ab(f,i) = 1.d0 + Sq_ab(f,i) * 4.d0 * pi * rho_ab * bin_width_g / q(i)
      Sq_bb(f,i) = 1.d0 + Sq_bb(f,i) * 4.d0 * pi * rho_bb * bin_width_g / q(i)

      ! write the data
      write(15,fmt="(*(g0:','))") framechoice(f), q(i), Sq_all(f,i), &
Sq_aa(f,i), Sq_ab(f,i), Sq_bb(f,i)

    enddo

  enddo

  ! close the g(r) file
  close(14)

  ! close the S(q) file
  close(15)

  ! deallocate before ending
  !deallocate(hist_bin_cc_g,bin_cc_density_g)
  ! Deallocate arrays
  if (allocated(hist_bin_ab)) then
    deallocate(hist_bin_ab, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif
  if (allocated(hist_bin_aa)) then
    deallocate(hist_bin_aa, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif
  if (allocated(hist_bin_bb)) then
    deallocate(hist_bin_bb, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif
  if (allocated(hist_bin_all)) then
    deallocate(hist_bin_all, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif

  if (allocated(ab_density)) then
    deallocate(ab_density, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif
  if (allocated(aa_density)) then
    deallocate(aa_density, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif
  if (allocated(bb_density)) then
    deallocate(bb_density, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif
  if (allocated(all_density)) then
    deallocate(all_density, stat=ierr)
    if (ierr /= 0) print *, "Error: Unable to deallocate hist_bin_ab"
  endif

end subroutine structure_factor


!###########################################
subroutine edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)
! calculate the edgelist (every colloid-colloid connection) for all frames in
! a colloid simulation

  implicit none

  !!! INPUTS
  integer, intent(in) :: nframes                                ! number of colloids
  integer, intent(in) :: ncolloids                              ! number of colloids
  real(8), intent(in) :: radii(0:ncolloids-1)                   ! colloid radii
  real(8), intent(in) :: allpos(0:nframes-1,0:ncolloids-1,0:2)  ! colloid positions
  !f2py depend(ncolloids) radii
  !f2py depend(nframes,ncolloids) allpos
  real(8), intent(in) :: lbox(0:2)                              ! simulation box size
  real(8), intent(in) :: rcut                                   ! surface-surface cut-off for a bond
  character*256, intent(in) :: edge_output                       ! output filepath for edgelist

  !!! OUTPUTS
  !N/A

  !!! INTERNAL VARIABLES
  real(8) :: rij(0:2), hij                 ! center-center (rij) and surface-surface (hij) interaction distance
  integer :: f                             ! frame index, f
  integer :: i,j                           ! colloid indices i, j
  character(len=10) :: file_id_frame       ! variable for frame number as string
  character(len=256) :: filename_edge       ! variable for frame-specific filepath/name

  ! for each frame
  do f=0,nframes-1
    ! write the frame number into a string for naming files
    write(file_id_frame, '(i0)') f

    ! construct the edgelist filename
    filename_edge = trim(adjustl(edge_output)) // trim(adjustl(file_id_frame)) // '.csv'
    ! create the file for this frame
    open(unit=14,file=trim(filename_edge),action='write')
    ! write the header
    write(14,fmt="(*(g0:','))") "i","j"

    ! caluclate the edge list, considering for all possible pairs
    do i=0,ncolloids-2
      do j=i+1,ncolloids-1
        ! calculate the center-center interparticle distance
        rij(0) = allpos(f,i,0) - allpos(f,j,0)
        rij(1) = allpos(f,i,1) - allpos(f,j,1)
        rij(2) = allpos(f,i,2) - allpos(f,j,2)
        ! adjust for interactions across the periodic boundaries
        ! NOTE: assuming that there is NO shear deforming the box
        rij = rij - lbox * dnint(rij/lbox)
        ! conver to surface surface distance
        hij = dsqrt(sum(rij*rij)) - (radii(i) + radii(j))
        ! if it is a bond
        if(hij<=rcut) then
          ! record in the edge list
          write(14,fmt="(*(g0:','))") i,j
        endif
      enddo
    enddo

    ! close the file
    close(14)

  enddo

end subroutine edgelist_calc



!###########################################
! Corrected kinetic temperature and pressure
! remove the solvent velocity from the system for one frame
subroutine corrected_temperature(nparticles,nlayerY,pos,vel,mass,typeid,Lbox,pxx,pxy,pxz,kT)
implicit none
integer, intent(in) :: nparticles, nlayerY
real(8), intent(in) :: pos(0:nparticles-1,0:2)
real(8), intent(in) :: vel(0:nparticles-1,0:2)
real(8), intent(in) :: Lbox(0:2)
real(8), intent(in) :: mass(0:nparticles-1)
integer, intent(in) :: typeid(0:nparticles-1)
!f2py depend(nparticles) pos,vel,mass,typeid
real(8), intent(out) :: pxx, pxy, pxz, kT
real(8) :: bin_width, inv_bin_width, vel_x
real(8) :: vel_bin(0:nlayerY-1)
integer :: hist_bin(0:nlayerY-1)
integer :: i, k
pxx=0.d0; pxy=0.d0; pxz=0.d0; kT=0.d0
hist_bin = 0
vel_bin = 0.d0
bin_width = Lbox(1)/nlayerY
inv_bin_width = 1.d0/bin_width

!print *, vel

! get the solvent velocity profile
do i=0,nparticles-1
  if(typeid(i) == 0) then
    k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
    if(k >= nlayerY) then
      k = nlayerY -1
   endif
    hist_bin(k) = hist_bin(k) + 1
    !print *, vel(i,0)
    vel_bin(k) = vel_bin(k) + vel(i,0)
  endif
enddo
do i=0,nlayerY-1
  if(hist_bin(i) == 0) then
    hist_bin(i) = 1
  endif
enddo

!print *, vel_bin
!print *, hist_bin
vel_bin = vel_bin/hist_bin

! remove effect of solvents
do i=0,nparticles-1
  k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
  if(k >= nlayerY) then
    k = nlayerY -1
  endif
  vel_x = vel(i,0) - vel_bin(k)
  !pxx = pxx + mass(i) * (vel_x * vel_x - vel(i,0)*vel(i,0))
  !pxy = pxy + mass(i) * (vel_x * vel(i,1) - vel(i,0)*vel(i,1))
  !pxz = pxz + mass(i) * (vel_x * vel(i,2) - vel(i,0)*vel(i,2))
  pxx = pxx + mass(i) * (vel_x * vel_x)
  pxy = pxy + mass(i) * (vel_x * vel(i,1)) 
  pxz = pxz + mass(i) * (vel_x * vel(i,2))
enddo

! ONLY remove effect of solvents from solvents
!do i=0,nparticles-1
!  if (typeid(i) == 0) then
!    k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
!    if(k >= nlayerY) then
!      k = nlayerY -1
!    endif
!    vel_x = vel(i,0) - vel_bin(k)
!    !pxx = pxx + mass(i) * (vel_x * vel_x - vel(i,0)*vel(i,0))
!    !pxy = pxy + mass(i) * (vel_x * vel(i,1) - vel(i,0)*vel(i,1))
!    !pxz = pxz + mass(i) * (vel_x * vel(i,2) - vel(i,0)*vel(i,2))
!    pxx = pxx + mass(i) * (vel_x * vel_x)
!    pxy = pxy + mass(i) * (vel_x * vel(i,1)) 
!    pxz = pxz + mass(i) * (vel_x * vel(i,2))
!  else:
!     pxx = pxx + mass(i) * (vel(i,0) * vel(i,0))
!     pxy = pxy + mass(i) * (vel(i,0) * vel(i,1)) 
!     pxz = pxz + mass(i) * (vel(i,0) * vel(i,2))
!enddo

! recalculate temperature
kT = pxx/3.d0/nparticles
! normalize pressure by the simulation box volume
pxx = pxx/Lbox(1)**3
pxy = pxy/Lbox(1)**3
pxz = pxz/Lbox(1)**3
end subroutine corrected_temperature


!###########################################
! Calculate the velocity profile of the selected particles for one frame
subroutine vel_profile(nparticles,nlayerY,pos,vel,Lbox,vel_bin,bin_array)
  implicit none

  !!! INPUTS
  integer, intent(in) :: nparticles, nlayerY
  real(8), intent(in) :: pos(0:nparticles-1,0:2)
  real(8), intent(in) :: vel(0:nparticles-1,0:2)
  real(8), intent(in) :: Lbox(0:2)
  !f2py depend(nparticles) pos,vel

  !!! OUTPUTS
  real(8), intent(out) :: vel_bin(0:nlayerY-1)
  real(8), intent(out) :: bin_array(0:nlayerY-1)

  !!! INTERNAL VARIABLES
  real(8) :: bin_width, inv_bin_width
  integer :: hist_bin(0:nlayerY-1)
  integer :: i, k

  bin_width = Lbox(1)/nlayerY
  inv_bin_width = 1.d0/bin_width
  hist_bin = 0
  vel_bin = 0.d0

  do i=0,nparticles-1
    k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
    if(k >= nlayerY) then
      k = nlayerY -1
    endif
    hist_bin(k) = hist_bin(k) + 1
    vel_bin(k) = vel_bin(k) + vel(i,0)
  enddo
  do i=0,nlayerY-1
    bin_array(i) = -0.5*Lbox(1)+(i+0.5d0)
    if(hist_bin(i) == 0) then
      hist_bin(i) = 1
    endif
  enddo
  vel_bin = vel_bin/hist_bin

end subroutine vel_profile


!###########################################
! Calculate the fabric tensor for selected particles in one frame
subroutine fabric_tensor(nparticles,pos,radius,Lbox,cut_off,m_xy,Rxy_avg)
implicit none
integer, intent(in) :: nparticles
real(8), intent(in) :: pos(0:nparticles-1,0:2)
real(8), intent(in) :: radius(0:nparticles-1)
!f2py depend(nparticles) pos,radius
real(8), intent(in) :: Lbox(0:2)
real(8), intent(in) :: cut_off, m_xy
real(8), intent(out) :: Rxy_avg
real(8) :: dr(0:2), r1, r2, inv_Lbox(0:2), img
integer :: i, j!, coord_num
inv_Lbox = 1.d0/Lbox
!coord_num = 0
Rxy_avg = 0.d0
do i=0,nparticles-2
  do j=i+1,nparticles-1
    dr = pos(i,:)-pos(j,:)
    dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
    img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
    dr(1) = dr(1) - img
    dr(0) = dr(0) - img * m_xy
    dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
    r2 = sum(dr*dr)
    r1 = dsqrt(r2) - radius(i) - radius(j)
    if(r1 <= cut_off) then
      Rxy_avg = Rxy_avg + dr(0)*dr(1)/r2
      !coord_num = coord_num + 1
    endif
  enddo
enddo
Rxy_avg = Rxy_avg/nparticles
end subroutine fabric_tensor



!###########################################
! Calculate the pore size for all frames using two methods: 
!   - Torquato’s Pore Size Distribution
!   - Gubbins’s Pore Size Distribution
!! NOTE: This is the most complicated code, it has MULTIPLE
!!       subroutines and uses a linked-list implementation
module pore_size_calculation
implicit none

! Total number of colloids in simulation box and number of frames in trajectory
integer, public :: ncolloid
integer, public :: nframes 

!Framewise particles positions and simulation box size
integer, allocatable, save, public :: framechoice(:)   ! the selected frames 
real(8), allocatable, save, public :: rxi(:,:)         ! framewise x-coordinates
real(8), allocatable, save, public :: ryi(:,:)         ! framewise y-cordinates
real(8), allocatable, save, public :: rzi(:,:)         ! framewise z-coordinates
real(8), allocatable, save, public :: box_length(:,:)  ! framewise box size
real(8), allocatable, save, public :: radii(:)         ! colloid radii, AKA the hard-sphere distance to avoid overlap between probe and colloid

real(8), public :: box_size(3)                         ! simulation box size in a frame
real(8), public :: inv_box_size(3)                     ! inverse of simulation box size in a frame
real(8), public :: rp(3)                               ! pore-size calculation location point
real(8), allocatable, public :: rpos(:,:)              ! colloid positions in a given frame

! Variables for linked-list method
real(8), public :: dcell_init                          ! initial cell-size
integer, save, public :: tot_cell                      ! total number of cells
integer, save, public :: ncell(3)                      ! number of cells in each direction
real(8), save, public:: dcell(3)                       ! cell size in each direction
integer, allocatable, save, public :: head(:,:,:)
integer, allocatable, save, public :: list(:)
integer, allocatable, save, public :: colloid_id(:,:)

contains
!------------------------------------------------------------------------

!subroutine for initial random seed

SUBROUTINE init_random_seed()
implicit none
integer :: i, n, clock
integer, allocatable:: seed(:)
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed
!--------------------------------------------------------------------------

subroutine pore_size_calc(data_outpath, pop_type, ncolloids_py, nframes_py, framechoice_py, &
nprobe, rhs_py, dcell_init_py, rxi_py, ryi_py, rzi_py, box_length_py)
implicit none

! Variables imported from Python file
character*256, intent(in) :: data_outpath
character*256, intent(in) :: pop_type
integer, intent(in) :: ncolloids_py, nframes_py, nprobe
integer, intent(in) :: framechoice_py(0:nframes_py-1)
real(8), intent(in) :: dcell_init_py
real(8), intent(in) :: rhs_py(ncolloids_py)
real(8), intent(in) :: rxi_py(0:nframes_py-1,0:ncolloids_py-1)
real(8), intent(in) :: ryi_py(0:nframes_py-1,0:ncolloids_py-1)
real(8), intent(in) :: rzi_py(0:nframes_py-1,0:ncolloids_py-1)
real(8), intent(in) :: box_length_py(0:nframes_py-1,3)
!f2py depend(ncolloids_py,nframes_py) rxi_py,ryi_py,rzi_py
!f2py depend(nframes_py,3) box_length_py
!f2py depend(nframes_py) framechoice_py
!f2py depend(ncolloids_py) rhs_py

! Internal Variables
integer :: i,j                       ! index for frame (i) and index for particle (j)
real(8) :: f_T, f_G                  ! Fraction of pore volume (as calculated by Torquato's PSD and Gubbin's PSD)
real(8) :: poresize_T, poresize_G    ! pore diameter as (as calculated by Torquato's PSD and Gubbin's PSD)
real(8) :: ri(3)                     ! reference location 
character(len=256) :: filename_pore   ! output filename/path

integer,pointer :: atemp => null()
call init_random_seed()

! assign imported variables to global variables
dcell_init = dcell_init_py
ncolloid = ncolloids_py
nframes = nframes_py
framechoice = framechoice_py

allocate( rxi(0:nframes-1,0:ncolloid-1) )
allocate( ryi(0:nframes-1,0:ncolloid-1) )
allocate( rzi(0:nframes-1,0:ncolloid-1) )
allocate( box_length(0:nframes-1,3) )
allocate( radii(0:ncolloid-1) )

radii(:) = rhs_py(:)

do i=0,nframes-1
  do j=0,ncolloid-1
    rxi(i,j) = rxi_py(i,j)
    ryi(i,j) = ryi_py(i,j)
    rzi(i,j) = rzi_py(i,j)
  enddo
  do j=1,3
    box_length(i,j) = box_length_py(i,j)
  enddo
enddo


! open the output file and write labels
filename_pore = trim(adjustl(data_outpath)) // '/poresize' // trim(adjustl(pop_type)) // '.csv'
open(unit=14,file=trim(filename_pore),action='write')
write(14,fmt="(*(g0:','))") 'frame', 'probe_posx','probe_posy', 'probe_posz', &
'voidcenter_x','voidcenter_y','voidcenter_z','porediameter_T','porediameter_G'

! calculate the pore size
do i=0,nframes-1
  ! read box_size in a given frame
  box_size(:) = box_length(i,:)
  inv_box_size(:) = 1.d0/box_size(:)

  ! read colloid positions for the given frame
  allocate(rpos(3,0:ncolloid-1))
  rpos(:,:)=0.d0
  do j=0,ncolloid-1 
    rpos(1,j)=rxi(i,j)
    rpos(2,j)=ryi(i,j)
    rpos(3,j)=rzi(i,j)
  enddo

  ! Linked list initialization, formation, and check
  call init_list()
  call link_list()
  call check_list()

  ! calculate pore size at every random prob point in the box 
  do j = 1,nprobe
    f_T=0.d0
    ! select random point in simulation box and check overlap with colloids
    do while( f_T == 0.d0)
      call random_number(rp(1))
      call random_number(rp(2))
      call random_number(rp(3))
      ! rescale position from (0,L) to (-L/2,L/2) coordinates to match sim box
      rp(1)=0.5d0*box_size(1)*(2.d0*rp(1)-1.d0)
      rp(2)=0.5d0*box_size(2)*(2.d0*rp(2)-1.d0)
      rp(3)=0.5d0*box_size(3)*(2.d0*rp(3)-1.d0)
      ! use the current position as the center of a pore
      ri(:) = rp(:)
      ! get the fraction of pore volume with Torquato's method 
      call fun(ri,f_T)
    enddo

    ! get the poresize with Torquato's method 
    poresize_T = 2.d0*dsqrt(-f_T)

    ! get the fraction of pore volume with Gubbin's method 
    f_G = 0.d0
    call solvopt(3, ri, f_G, fun, .false., atemp, .true., func, .false., atemp)
    ! get the poresize with Gubbin's method 
    poresize_G = 2.d0*dsqrt(-f_G)

    ! save the data for both poresize calculations 
    write(14,fmt="(*(g0:','))") framechoice(i), rp, ri, poresize_T, poresize_G
    !write(14,'(I8, 8f16.5)')  framechoice(i), rp, ri, poresize_T, poresize_G
  enddo
  deallocate(rpos)
  call finalize_list()
enddo
close(14)
end subroutine pore_size_calc
!--------------------------------------------

!Linked-list initialization - allocation of head and list arrays

subroutine init_list()
implicit none

! calculate the number of cells
ncell(:) = floor(box_size(:) / dcell_init)
! check number of cells
if (any(ncell(:) < 3)) then
  print *, 'system is too small to use cell links'
  stop
endif
! calculate the total number of cells
tot_cell = ncell(1) * ncell(2) * ncell(3)
! calculate the cell dimensions
dcell(:) = box_size(:) / dble(ncell(:))
! initialize variables for head position, list of colloids, and colloid IDs
allocate( head(0 : ncell(1)-1, 0 : ncell(2)-1, 0 : ncell(3)-1) )
allocate( list(0:ncolloid-1) )
allocate( colloid_id(3,0:ncolloid-1) )
head(:,:,:) = 0
list(:) = 0
colloid_id(:,:) = 0
return
end subroutine init_list
!--------------------------------------------------------------------------

!Compute cell numbers of monomers

subroutine list_cell_i(ri,cell_i)
implicit none
real(8), intent(in) :: ri(3)         ! colloid position
integer, intent(out) :: cell_i(3)    ! colloid index, cell at this index

if (any(dabs(ri(:)/box_size(:)) > 0.5d0 )) then
  print *, 'colloid not in the main-box'
  stop
end if

cell_i(:) = 0
cell_i(:) = floor( (ri(:)/box_size(:)+0.5d0) * dble(ncell(:)) )
cell_i(:) = modulo( cell_i(:), ncell(:) )
return
end subroutine list_cell_i
!---------------------------------------------------------------------------

!Formation of head and list arrays

subroutine  link_list()
implicit none
real(8) :: ri(3)
integer :: i, cell_i(3)

do i=0,ncolloid-1
  ri(:) = rpos(:,i)
  cell_i(:) = 0
  call list_cell_i(ri,cell_i)
  list(i) = head(cell_i(1),cell_i(2),cell_i(3))
  head(cell_i(1),cell_i(2),cell_i(3)) = i
  colloid_id(:,i) = cell_i(:)
enddo
return
end subroutine link_list
!----------------------------------------------------------------------

!Check formation of head and list arrays

subroutine check_list()
implicit none
integer :: i, j, k      ! indices for x (i), y (j), and z (k) dimensions of each cell
integer :: c            ! the index of the current cell
integer :: cell_i(3)    ! the current cell
real(8) :: ri(3)        ! individual colloid position 
do i=0,ncolloid-1
  ! read in the colloid position
  ri(:) = rpos(:,i)
  ! check the current cell has the correct contents
  cell_i(:) = 0
  call list_cell_i(ri,cell_i)
  if(any(cell_i(:) .ne. colloid_id(:,i))) then
    print *, 'inconsistency1 found, list_cell_i does not returned assigned colloidIDs:', i, cell_i, colloid_id(:,i)
    stop
  endif
enddo
do i=0, ncell(1)-1
  do j=0, ncell(2)-1
    do k=0, ncell(3)-1
      cell_i(:) = (/i,j,k/)
      c = head(i,j,k)
      do while ( c .ne. 0)
        if(any( cell_i(:) .ne. colloid_id(:,c))) then
          print *, 'inconsistency2 found, c=head() does not return assigned colloidIDs', c, cell_i(:), colloid_id(:,c)
          stop
        endif
        c = list(c)
      enddo
    enddo
  enddo
enddo
return
end subroutine check_list
!-----------------------------------------------------------------------

!Deallocation of arrays

subroutine finalize_list()
implicit none
deallocate(list)
deallocate(head)
deallocate(colloid_id)
end subroutine finalize_list
!------------------------------------------------------------------------

!Subroutine to compute pore size
subroutine fun(rin,f)
implicit none

! Inputs
real(8),intent(in) :: rin(3)      ! center of the pore being calculated 

! Outputs
real(8), intent(out) :: f         ! the calculated fraction of pore volume

! Internal Variables
integer :: c                      ! the index for the particle tag in the current cell
integer :: i, j, k                ! indices for x,y,z cell coordinates
integer :: new_i, new_j, new_k    ! updated location of the cell
integer :: cell_i(3)              ! the current cell
real(8) :: dist                   ! scalar r_ij distance between the pore center and the nearest particles
real(8) :: rij(3)                 ! the center-center distance between the pore center and the nearest particles
real(8) :: ri(3)                  ! location of the center of the pore (AKA rc_in)
real(8) :: pore_size              ! current pore_size estimate (updated during the calculation)
real(8) :: max_size               ! maximum pore size (the size of the simulation box) AKA initial pore_size estimate

f=0.d0
! initialize the pore_size as it's max possible value (sim box volume)
max_size = dmax1(box_size(1),box_size(2),box_size(3))
pore_size = max_size

! read the current pore location being checked
ri = rin
! adjust the position to include interactions across the periodic boundaries
ri(:) = ri(:) - box_size(:) * dnint(ri(:)*inv_box_size(:))

! initialize the current cell at zero
cell_i(:) = 0

! find the cell that contains the current pore being checked
call list_cell_i(ri,cell_i)

! check/update cell information with neighbors
do i = cell_i(1)-1, cell_i(1)+1
  do j = cell_i(2)-1, cell_i(2)+1
    do k = cell_i(3)-1, cell_i(3)+1
      new_i = modulo(i,ncell(1))
      new_j = modulo(j,ncell(2))
      new_k = modulo(k,ncell(3))
      c = head(new_i,new_j,new_k)

      ! for each colloid in the cell
      do while(c .ne. 0)
        ! calculate the distance between the center of the pore and the colloid 
        rij = rpos(:,c) - ri
        ! update the position for interactions across the sim's periodic boundaries
        rij(:) = rij(:) - box_size(:)* dnint(rij(:)*inv_box_size(:))
        ! convert to scalar and subtract the colloid radius
        dist = dsqrt(sum(rij*rij)) - radii(c)
        ! if the probe is inside a colloid, then dist<0 (and there is no pore)
        if(dist .lt. 0.d0) then
          return
        ! else, if the distance is less than the current pore size (initialized at max), update the size
        else if (dist .lt. pore_size) then
          pore_size=dist     
        end if
        ! update c and continue checking through the remaining particles in the cell
        c = list(c)
      enddo
    enddo
  enddo
enddo
! if the pore size is still larger than the size of a cell, check against all colloids 
if(pore_size .ge. dcell_init) then
  do i = 0,ncolloid-1
    rij = rpos(:,i) - ri
    rij(:) = rij(:) - box_size(:) * dnint(rij(:)*inv_box_size(:))
    dist = dsqrt(sum(rij*rij)) - radii(i)
    if (dist .lt. pore_size) then
      pore_size = dist
    endif
  enddo
endif
! use pore_size to update f
if(pore_size .gt. 0) f = -(pore_size) * (pore_size)
return
end subroutine fun
!----------------------------------------------------------------------------

!Subroutine to compute constraint value
subroutine func(rin,f)
implicit none

! Inputs
real(8),intent(in) :: rin(3)  ! the center of the current pore

! Outputs
real(8), intent(out) :: f     ! the calculated fraction of pore volume

! Internal Values
real(8) :: f_local, dr2       ! the input fraction of pore volume, the constraint value

! get the Torquato's pore size
call fun(rin,f_local)

! check the constraint value (should be less than or equal to 0)
dr2 = sum((rp-rin)*(rp-rin)) + f_local
f = max(0.d0,dr2)
return

end subroutine func
!----------------------------------------------------------------------------


! Add a subroutine to deallocate arrays
subroutine deallocate_arrays()
  deallocate(rxi, ryi, rzi, box_length, radii)
end subroutine deallocate_arrays
!----------------------------------------------------------------------------

end module pore_size_calculation
                 
