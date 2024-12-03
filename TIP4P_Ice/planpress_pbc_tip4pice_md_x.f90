! ==========================================================================
! Calculation of planar pressure tensor for pure rigid molecules (eg TIP4P/ICE):
!     
!   - Irving-Kirkwood contour is used
!   - Periodic boundary conditions are applied
!   - Wolf potential (damped and cut) is applied to reproduce the Ewald method used in simulation
!   - Input file is xtc format from MD (GROMACS) simulation 
!   - x-axis is always normal to the planar interface
!   - full tensor calculations are made available (07/03/2023)
!  
! Author: Kaihang Shi
! Email: kaihangshi0@gmail.com
! ==========================================================================

PROGRAM virialpress_planar

! use an external module for analyzing the GROMACS trajectory
use gmxfort_trajectory

IMPLICIT NONE

! ---------------- Input Parameters ---------------------------!
! Control parameters 
Integer, Parameter ::  st_frame=119401
Integer, Parameter ::  nd_frame=120000

! xtc file PATH 
Character*256, Parameter :: xtcpath = "/vscratch/grp-kaihangs/kaihangs/traj_ice/pPrim_stable/&
                                       &/traj_50ns.xtc"  
! solid cluster center-of-mass file (evaluated using oxygen atom only)
Character*256, Parameter :: compath = "/vscratch/grp-kaihangs/kaihangs/traj_ice/pPrim_stable/&
                                       &/CM_coordinates.dat"

! System parameter 
Double Precision, Parameter :: temp = 270.0d0      ! system temperature [Kelvin]
Integer, Parameter :: n_mol_tot = 20736              ! total number of molecules in the system
Integer, Parameter :: n_sites_tot = 82944           ! Total number of atoms in the system
Integer, Parameter :: n_mol_sites = 4               ! number of sites in the water molecule    

! Intermolecular parameter 
Double Precision, Parameter :: sigma_o = 3.1668d0      ! sigma for oxygen atom [Angstrom] 
Double Precision, Parameter :: epsilonkb_o = 106.1d0    ! epsilon/k_B for oxygen atom [Kelvin]
Double Precision, Parameter :: q_m = -1.1794d0        ! point charge for M-site, TIP4P/ICE [e]
Double Precision, Parameter :: q_h = 0.5897d0        ! point charge for H, TIP4P/ICE model [e]
Double Precision, Parameter :: r_ljcut = 9.0d0      ! LJ cutoff radius [Angstrom]
Double Precision, Parameter :: alp = 0.195d0       ! Alpha Parameter for Wolf potential method [1/Angstrom]
Double Precision, Parameter :: r_coulcut = 20.0d0     ! cutoff radius for Wolf potential method [Angstrom]

! Mass 
Double Precision, Parameter :: mol_mass = 18.01528d0  ! molecular weight of water [g/mol]
Double Precision, Parameter :: mass_o = 15.999400
Double Precision, Parameter :: mass_h = 1.007940

! Density/pressure profile resolution
Double Precision, Parameter :: rden_cut = 120.0d0  ! Cutoff distance for density/pressure calculation [Angstrom]
Integer, Parameter :: n_bins = 240                  ! density profile resolution 


! ---------------- Define Variables -------------------------!
! ----- XTC variables
! Box dimension
Real, Dimension(3,3) :: mybox
Double Precision, Dimension(:,:), Allocatable :: boxsz
! Box origin (to be set to the center of a sphere/box)
Double Precision, Dimension(:), Allocatable :: ox, oy, oz
! Read in parameter
Integer :: frame_id, xtcframes, ns_frame, c_frame
Integer :: natom
! Trajectory data type
Type(Trajectory) :: trj
! Coordinates
Real, Dimension(3) :: myatom
Double Precision, Dimension(:,:,:), Allocatable :: rx_s, ry_s, rz_s
Double Precision, Dimension(:,:), Allocatable :: rx, ry, rz

! Statistics
! Only fluid-fluid contribution is considered here
Double Precision, Dimension(:), Allocatable :: virialpress_pxx
Double Precision, Dimension(:), Allocatable :: virialpress_pxy
Double Precision, Dimension(:), Allocatable :: virialpress_pxz
Double Precision, Dimension(:), Allocatable :: virialpress_pyx
Double Precision, Dimension(:), Allocatable :: virialpress_pyy
Double Precision, Dimension(:), Allocatable :: virialpress_pyz
Double Precision, Dimension(:), Allocatable :: virialpress_pzx
Double Precision, Dimension(:), Allocatable :: virialpress_pzy
Double Precision, Dimension(:), Allocatable :: virialpress_pzz
Double Precision, Dimension(:), Allocatable :: virialpxxavg
Double Precision, Dimension(:), Allocatable :: virialpxyavg
Double Precision, Dimension(:), Allocatable :: virialpxzavg
Double Precision, Dimension(:), Allocatable :: virialpyxavg
Double Precision, Dimension(:), Allocatable :: virialpyyavg
Double Precision, Dimension(:), Allocatable :: virialpyzavg
Double Precision, Dimension(:), Allocatable :: virialpzxavg
Double Precision, Dimension(:), Allocatable :: virialpzyavg
Double Precision, Dimension(:), Allocatable :: virialpzzavg
Double Precision, Dimension(:), Allocatable :: prkin

! bin width
Double Precision :: bin_wid, bin_dvol
Double Precision, Dimension(n_bins) :: bin_posp, bin_posn
! Statistics 
Double Precision, Dimension(:), Allocatable :: denavg
Double Precision, Dimension(:,:), Allocatable :: denblk

!Variables for pressure calculation
INTEGER :: imol, jmol, ibin, iframe, ierr
INTEGER :: isite, jsite

DOUBLE PRECISION :: rxi, ryi, rzi, rij, rijsq, rxj, ryj, rzj, xij, yij, zij, rxii, ryii, rzii
DOUBLE PRECISION :: rxis, ryis, rzis, rxjs, ryjs, rzjs
DOUBLE PRECISION :: rijs, rijssq,xijs, yijs, zijs
DOUBLE PRECISION :: xijssq, yijssq, zijssq, xijsq, yijsq, zijsq

DOUBLE PRECISION :: sr3, sr6, sr12, phitz, alppi, alpsq, f_lj, f_coul, f_coul_rc,f_tot
DOUBLE PRECISION :: pxx, pxy, pxz, pyx, pyy, pyz, pzx, pzy, pzz

Double Precision :: solk
Double Precision :: area
! CPU time
Double Precision :: cpu_st, cpu_nd

! Constants
Double Precision, Parameter :: Pi = 3.141592653589d0
Double Precision, parameter :: two_Pi = 2.0d0*Pi
Double Precision, Parameter :: Kb = 1.3806488d-23   ! [J/K] 
Double Precision, Parameter :: Na = 6.02214129d23  ! [1/mol] */
Double Precision, Parameter :: h = 6.62606957d-34   ! [J*s] = [kg*m^2/s]
Double Precision, Parameter :: R = Kb*Na            ! [J/(mol*K)] */
Double Precision, Parameter :: EETOK = 167100.9558738398d0 ! EETOK = 1/(4*pi*epsilon*kb) * e^2 * 10^10
Double Precision, Parameter :: PVTOK = 0.7242971565d-02   ! PVTOK = 1.0d5*1.0d-30/kb , [bar]*[A^3] to [K]
Double Precision, Parameter :: PCOEFF = Kb*1.0d25        ! [K]/[A]^3 to [bar]

! -------------------------------------------------------------------------------------------------! 

! -------- Precalculate some variabels -----!
! total number of frames for sampling purpose
ns_frame = nd_frame - st_frame + 1

! -------- Allocate variables ----------!
Allocate(boxsz(3,ns_frame), STAT=ierr)
Allocate(ox(ns_frame), STAT=ierr)
Allocate(oy(ns_frame), STAT=ierr)
Allocate(oz(ns_frame), STAT=ierr)

Allocate(rx_s(n_mol_sites,n_mol_tot,ns_frame), STAT=ierr)
Allocate(ry_s(n_mol_sites,n_mol_tot,ns_frame), STAT=ierr)
Allocate(rz_s(n_mol_sites,n_mol_tot,ns_frame), STAT=ierr)

Allocate(rx(n_mol_tot,ns_frame), STAT=ierr)
Allocate(ry(n_mol_tot,ns_frame), STAT=ierr)
Allocate(rz(n_mol_tot,ns_frame), STAT=ierr)

Allocate(virialpress_pxx(n_bins), STAT=ierr) 
Allocate(virialpress_pxy(n_bins), STAT=ierr) 
Allocate(virialpress_pxz(n_bins), STAT=ierr) 

Allocate(virialpress_pyx(n_bins), STAT=ierr) 
Allocate(virialpress_pyy(n_bins), STAT=ierr) 
Allocate(virialpress_pyz(n_bins), STAT=ierr) 

Allocate(virialpress_pzx(n_bins), STAT=ierr) 
Allocate(virialpress_pzy(n_bins), STAT=ierr) 
Allocate(virialpress_pzz(n_bins), STAT=ierr) 

Allocate(virialpxxavg(n_bins), STAT=ierr) 
Allocate(virialpxyavg(n_bins), STAT=ierr) 
Allocate(virialpxzavg(n_bins), STAT=ierr) 

Allocate(virialpyxavg(n_bins), STAT=ierr) 
Allocate(virialpyyavg(n_bins), STAT=ierr) 
Allocate(virialpyzavg(n_bins), STAT=ierr) 

Allocate(virialpzxavg(n_bins), STAT=ierr) 
Allocate(virialpzyavg(n_bins), STAT=ierr) 
Allocate(virialpzzavg(n_bins), STAT=ierr) 

Allocate(prkin(n_bins), STAT=ierr) 

Allocate(denavg(n_bins), STAT=ierr) 
Allocate(denblk(n_bins,ns_frame), STAT=ierr) 


! -- Initialize parameter -----!
virialpress_pxx = 0.0d0
virialpress_pxy = 0.0d0
virialpress_pxz = 0.0d0

virialpress_pyx = 0.0d0
virialpress_pyy = 0.0d0
virialpress_pyz = 0.0d0

virialpress_pzx = 0.0d0
virialpress_pzy = 0.0d0
virialpress_pzz = 0.0d0

virialpxxavg = 0.0d0
virialpxyavg = 0.0d0
virialpxzavg = 0.0d0

virialpyxavg = 0.0d0
virialpyyavg = 0.0d0
virialpyzavg = 0.0d0

virialpzxavg = 0.0d0
virialpzyavg = 0.0d0
virialpzzavg = 0.0d0

prkin = 0.0d0
denavg = 0.0d0
denblk = 0.0d0
bin_wid = rden_cut/DBLE(n_bins)

! Calculate constants in planar geometry
Do ibin = 1, n_bins
  ! Get y-distance of ibin 
  bin_posp(ibin) = (DBLE(ibin)-0.5d0)*bin_wid
  bin_posn(ibin) = -bin_posp(ibin)
End Do


! Calculate parameters in shifted-force Wolf potential
alpsq = alp**2
alppi = 2.0d0*alp/dSQRT(Pi)
f_coul_rc = dERFC(alp*r_coulcut)/r_coulcut**2 + alppi*dEXP(-alpsq*r_coulcut**2)/r_coulcut



! --------- Write header ------------!
Write(*,'(A)') '==========================================================  '
Write(*,'(A)') ' This program is used to calculate the pressure tensor for  ' 
Write(*,'(A)') ' planar geometry from virial route                     '
Write(*,'(A)') ' Only coded for TIP4P/Ice model                         '
Write(*,'(A)') ' Using the Irving-Kirkwood (IK) definition               '
Write(*,'(A)') '                                                            '  
Write(*,'(A)') ' Author: Kaihang Shi                                     '    
Write(*,'(A)') ' ========================================================== '
Write(*,'(A)') '                                                            '

! Creat a file recording running progress
Open(5,File='progress',action='write', status='replace')
! Open COM file to read
Open(7,File=compath,Status='Old',Access='Sequential',Action= 'Read')

! ---- Read in coordinates from GROMACS XTC format file ------!
! Open xtc file
Call trj%open(xtcpath)

! Get total number of atoms 
natom = trj%natoms()
! Get total number of frames 
xtcframes = trj%nframes

! Make sure this run is sensible
If (natom .NE. n_sites_tot) Then
  Write(*,'(A,I8)') 'True total number of atoms in the system is', natom
  Write(*,*) 'Reset n_sites_tot!'
  STOP
End If

! initialize counter
c_frame = 1

! Loop over pre-specified frames 
Do iframe = 1, nd_frame

  ! Read in useful frames
  If (iframe .GE. st_frame) Then

    ! Read in one frame
    frame_id = trj%read_next()

    ! Get box dimension for the current frame
    mybox = trj%box(1)
    ! Reassign box dimension to local variables (and convert from [nm] to [Angstrom])
    boxsz(1,c_frame) = mybox(1,1)*10.0
    boxsz(2,c_frame) = mybox(2,2)*10.0
    boxsz(3,c_frame) = mybox(3,3)*10.0

    ! the water bond is 1 A, to make sure both molecular COM and atom positon satisfy the minimum image convention
    If ((r_coulcut+2.0d0) .GT. Min(boxsz(1,c_frame),boxsz(2,c_frame),boxsz(3,c_frame))/2.0d0) Then
      Write(5,*) 'r_coulcut is too large for this box size!'
      STOP
    End If

    ! Read in COM (in nanometer)
    Read(7,*) ox(c_frame), oy(c_frame), oz(c_frame)

    ! Reassign coordinates into new array according to the file structure 
    jsite = 1
    ! loop over all molecules
    Do imol = 1, n_mol_tot
      ! loop over sites on each molecule
      Do isite = 1, n_mol_sites

        ! Get atom's coordinates in [nanometer]
        myatom = trj%x(1,jsite)

        ! translate coordinates using center of neucleus as the origin and convert units to [Angstrom]
        ! -------------------------------------------------------------------------
        ! site type sequence for TIP4P/ice: O, H, H, M 
        ! -------------------------------------------------------------------------
        rx_s(isite,imol,c_frame) = (myatom(1) - ox(c_frame))*10.0
        ry_s(isite,imol,c_frame) = (myatom(2) - oy(c_frame))*10.0
        rz_s(isite,imol,c_frame) = (myatom(3) - oz(c_frame))*10.0

        ! update site index
        jsite = jsite + 1 
      End do

      ! Calculation of COM of molecules
      rx(imol,c_frame) = (rx_s(1,imol,c_frame)*mass_o + (rx_s(2,imol,c_frame)+rx_s(3,imol,c_frame))*mass_h)/(mass_o+2.0*mass_h)
      ry(imol,c_frame) = (ry_s(1,imol,c_frame)*mass_o + (ry_s(2,imol,c_frame)+ry_s(3,imol,c_frame))*mass_h)/(mass_o+2.0*mass_h)
      rz(imol,c_frame) = (rz_s(1,imol,c_frame)*mass_o + (rz_s(2,imol,c_frame)+rz_s(3,imol,c_frame))*mass_h)/(mass_o+2.0*mass_h)

    End Do

    ! Put COM/sites into the central box with box centered at (0,0,0)
    Do imol = 1, n_mol_tot
      ! Put COM into the central box
      rx(imol,c_frame) = rx(imol,c_frame) - dNINT(rx(imol,c_frame)/boxsz(1,c_frame))*boxsz(1,c_frame) 
      ry(imol,c_frame) = ry(imol,c_frame) - dNINT(ry(imol,c_frame)/boxsz(2,c_frame))*boxsz(2,c_frame) 
      rz(imol,c_frame) = rz(imol,c_frame) - dNINT(rz(imol,c_frame)/boxsz(3,c_frame))*boxsz(3,c_frame) 

      ! this is redundent and for visualization only
      Do isite = 1, n_mol_sites

        ! Apply minimum image convention, move all atoms to the central box
        rx_s(isite,imol,c_frame) = rx_s(isite,imol,c_frame) - dNINT(rx_s(isite,imol,c_frame)/boxsz(1,c_frame))*boxsz(1,c_frame)
        ry_s(isite,imol,c_frame) = ry_s(isite,imol,c_frame) - dNINT(ry_s(isite,imol,c_frame)/boxsz(2,c_frame))*boxsz(2,c_frame)
        rz_s(isite,imol,c_frame) = rz_s(isite,imol,c_frame) - dNINT(rz_s(isite,imol,c_frame)/boxsz(3,c_frame))*boxsz(3,c_frame)
        
      End Do
    End Do

    ! Update counter
    c_frame = c_frame + 1 


  ! Skip useless frames
  Else
    frame_id = trj%read_next()
    Read(7,*) 
  ENDIF

! End reading specified frames
End Do 

! Close xtc file
Call trj%close()

Write(5,*) 'Finished reading coordinates and box information from XTC file!'
FLUSH(5)
Close(5)
Close(7)
! unit=6 is screen
FLUSH(6)



Call CPU_TIME(cpu_st)
! --------------- Start postprocessing ----------------------!
Do iframe = 1, ns_frame

  
  ! assuming x-direction is perpendicular to the planar interface
  area = boxsz(2,iframe)*boxsz(3,iframe)   ! Ly  * Lz
  bin_dvol = bin_wid*area

    
  ! ----------- Sampling x-density ------------!
  ! Loop over all molecules
  Do imol = 1, n_mol_tot

    if (dABS(rx(imol,iframe)) .LT. rden_cut) Then

      ibin = FLOOR(rx(imol,iframe)/bin_wid)

      if (ibin .LT. 0) Then
        ibin = - ibin
      Else if (ibin .GE. 0) Then
        ibin = ibin + 1
      End If

      ! Accumulate number 
      denblk(ibin,iframe) = denblk(ibin,iframe) + 1.0d0

    End If

  ! End loop over all molecules 
  End Do


  ! Loop over bins
  Do ibin = 1, n_bins
    
    ! Convert to number density (1/A^3)
    ! divided by 2 to average over two interfaces
    denavg(ibin) = denavg(ibin) + denblk(ibin,iframe)/(2.0d0 * bin_dvol)
        
  End Do  
  

  ! ----- Sampling pressure tensor -------!   
  ! Loop over all molecules
  Do imol = 1, n_mol_tot-1

    rxi = rx(imol,iframe) 
    ryi = ry(imol,iframe) 
    rzi = rz(imol,iframe)

    Do jmol = imol +1, n_mol_tot

      rxj = rx(jmol,iframe) 
      ryj = ry(jmol,iframe) 
      rzj = rz(jmol,iframe)

      !Calculate vector between sites
      xij = rxj - rxi
      yij = ryj - ryi
      zij = rzj - rzi

      ! Apply minimum image convention
      xij = xij - dNINT(xij/boxsz(1,iframe))*boxsz(1,iframe)
      yij = yij - dNINT(yij/boxsz(2,iframe))*boxsz(2,iframe)
      zij = zij - dNINT(zij/boxsz(3,iframe))*boxsz(3,iframe)

      !Square the values
      xijsq = xij**2
      yijsq = yij**2
      zijsq = zij**2
      rijsq = xijsq + yijsq + zijsq

      !Calculate distance between i and j molecule
      rij = dSQRT(rijsq)

      ! Check with force cutoff (H-H interaction requires rij<=2*r_OH+r_cut) 
      If (rij .GT. (r_coulcut + 2.0d0)) CYCLE


      ! initialize full tensor
      pxx = 0.0d0 
      pxy = 0.0d0 
      pxz = 0.0d0 
      pyx = 0.0d0 
      pyy = 0.0d0 
      pyz = 0.0d0 
      pzx = 0.0d0 
      pzy = 0.0d0 
      pzz = 0.0d0 

      ! loop over sites on imol
      DO isite = 1, n_mol_sites

        rxis = rx_s(isite,imol,iframe) 
        ryis = ry_s(isite,imol,iframe) 
        rzis = rz_s(isite,imol,iframe)


        DO jsite = 1, n_mol_sites

          rxjs = rx_s(jsite,jmol,iframe) 
          ryjs = ry_s(jsite,jmol,iframe) 
          rzjs = rz_s(jsite,jmol,iframe)

          !Calculate vector between sites
          xijs = rxjs - rxis
          yijs = ryjs - ryis
          zijs = rzjs - rzis

          ! Apply minimum image convention
          xijs = xijs - dNINT(xijs/boxsz(1,iframe))*boxsz(1,iframe)
          yijs = yijs - dNINT(yijs/boxsz(2,iframe))*boxsz(2,iframe)
          zijs = zijs - dNINT(zijs/boxsz(3,iframe))*boxsz(3,iframe)

          !Square the values
          xijssq = xijs*xijs
          yijssq = yijs*yijs
          zijssq = zijs*zijs

          !Calculate distance between i and j sites
          rijssq = xijssq + yijssq + zijssq
          rijs = dSQRT(rijssq)

          ! assuming r_coulcut > r_ljcut always
          If (rijs .GT. r_coulcut) CYCLE

          ! Determine interaction
          f_tot = 0.0d0

          ! O-O interaction (only LJ)
          If( (isite .eq. 1) .and. (jsite .eq. 1) ) Then

            If (rijs .LE. r_ljcut) Then

              !Calculate 12-6LJ force (based on sites) [K/A]
              sr3=(sigma_o/rijs)**3
              sr6=sr3*sr3
              sr12=sr6*sr6
              phitz = 24.0*epsilonkb_o/rijs
              f_lj = - phitz*(sr6 - 2.0d0*sr12)
              f_tot = f_lj

            End if

          ! H-H interaction (only Coulomb)
          Else if ( ((isite .eq. 2) .or. (isite .eq. 3)) .and. ((jsite .eq. 2) .or. (jsite .eq. 3)) ) Then

            ! Calculate Coulombic force using shifted-version of Wolf potential [K/A]
            f_coul = q_h*q_h*( (dERFC(alp*rijs)/rijssq + alppi*dEXP(-alpsq*rijssq)/rijs) - f_coul_rc )*EETOK
            f_tot = f_coul

          ! H-M interaction (only Coulomb)
          Else if ( ( ((isite .eq. 2) .or. (isite .eq. 3)) .and. (jsite .eq. 4) ) .or. &
            & ( ((jsite .eq. 2) .or. (jsite .eq. 3)) .and. (isite .eq. 4) ) ) Then

            ! Calculate Coulombic force using shifted-version of Wolf potential [K/A]
            f_coul = q_m*q_h*( (dERFC(alp*rijs)/rijssq + alppi*dEXP(-alpsq*rijssq)/rijs) - f_coul_rc )*EETOK
            f_tot = f_coul

          ! M-M interaction (only Coulomb)
          Else if ( (isite .eq. 4) .and. (jsite .eq. 4) ) Then

            ! Calculate Coulombic force using shifted-version of Wolf potential [K/A]
            f_coul = q_m*q_m*( (dERFC(alp*rijs)/rijssq + alppi*dEXP(-alpsq*rijssq)/rijs) - f_coul_rc )*EETOK
            f_tot = f_coul

          End If

          If (f_tot .eq. 0.0d0) CYCLE


          ! configurational part of the pressure
          pxx = pxx + xijs * xij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          pxy = pxy + xijs * yij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          pxz = pxz + xijs * zij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)

          pyx = pyx + yijs * xij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          pyy = pyy + yijs * yij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          pyz = pyz + yijs * zij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          
          pzx = pzx + zijs * xij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          pzy = pzy + zijs * yij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)
          pzz = pzz + zijs * zij * f_tot / (rijs * dABS(xij)) * (1.0d0/area)


        !End loop over jmol sites 
        ENDDO
        
      !End loop over imol sites
      ENDDO


      ! Special treatment of PBC (only account for y-direction which is perpendicular to the planar interface)
      If (xij .NE. (rxj -rxi)) Then

        ! -- Assume imol in the central box and ij pair contribue once
        ! == Irving-Kirkwood definition ==
        !Loop through the bins on the POSITIVE x side
        DO ibin= 1, n_bins
          
          ! Calculate intersection with planes ibin 
          solk = (bin_posp(ibin) - rxi)/xij

          if((solk .LT. 0.0d0) .OR. (solk .GT. 1.0d0)) CYCLE

          virialpress_pxx(ibin) = virialpress_pxx(ibin) + pxx
          virialpress_pxy(ibin) = virialpress_pxy(ibin) + pxy
          virialpress_pxz(ibin) = virialpress_pxz(ibin) + pxz

          virialpress_pyx(ibin) = virialpress_pyx(ibin) + pyx
          virialpress_pyy(ibin) = virialpress_pyy(ibin) + pyy
          virialpress_pyz(ibin) = virialpress_pyz(ibin) + pyz

          virialpress_pzx(ibin) = virialpress_pzx(ibin) + pzx
          virialpress_pzy(ibin) = virialpress_pzy(ibin) + pzy
          virialpress_pzz(ibin) = virialpress_pzz(ibin) + pzz

        !End bin cycle
        ENDDO

        !Loop through the bins on the NEGATIVE x side
        DO ibin= 1, n_bins
          
          ! Calculate intersection with planes ibin 
          solk = (bin_posn(ibin) - rxi)/xij

          if((solk .LT. 0.0d0) .OR. (solk .GT. 1.0d0)) CYCLE

          virialpress_pxx(ibin) = virialpress_pxx(ibin) + pxx
          virialpress_pxy(ibin) = virialpress_pxy(ibin) + pxy
          virialpress_pxz(ibin) = virialpress_pxz(ibin) + pxz

          virialpress_pyx(ibin) = virialpress_pyx(ibin) + pyx
          virialpress_pyy(ibin) = virialpress_pyy(ibin) + pyy
          virialpress_pyz(ibin) = virialpress_pyz(ibin) + pyz

          virialpress_pzx(ibin) = virialpress_pzx(ibin) + pzx
          virialpress_pzy(ibin) = virialpress_pzy(ibin) + pzy
          virialpress_pzz(ibin) = virialpress_pzz(ibin) + pzz

        !End bin cycle
        ENDDO


        ! -- assume jmol in the central box and ij pair contribute again
        rxii = rxj - xij
        ! == Irving-Kirkwood definition ==
        !Loop through the bins on the POSITIVE x side
        DO ibin= 1, n_bins
          
          ! Calculate intersection with planes ibin 
          solk = (bin_posp(ibin) - rxii)/xij

          if((solk .LT. 0.0d0) .OR. (solk .GT. 1.0d0)) CYCLE

          virialpress_pxx(ibin) = virialpress_pxx(ibin) + pxx
          virialpress_pxy(ibin) = virialpress_pxy(ibin) + pxy
          virialpress_pxz(ibin) = virialpress_pxz(ibin) + pxz

          virialpress_pyx(ibin) = virialpress_pyx(ibin) + pyx
          virialpress_pyy(ibin) = virialpress_pyy(ibin) + pyy
          virialpress_pyz(ibin) = virialpress_pyz(ibin) + pyz

          virialpress_pzx(ibin) = virialpress_pzx(ibin) + pzx
          virialpress_pzy(ibin) = virialpress_pzy(ibin) + pzy
          virialpress_pzz(ibin) = virialpress_pzz(ibin) + pzz

        !End bin cycle
        ENDDO

        !Loop through the bins on the NEGATIVE x side
        DO ibin= 1, n_bins
          
          ! Calculate intersection with planes ibin 
          solk = (bin_posn(ibin) - rxii)/xij

          if((solk .LT. 0.0d0) .OR. (solk .GT. 1.0d0)) CYCLE

          virialpress_pxx(ibin) = virialpress_pxx(ibin) + pxx
          virialpress_pxy(ibin) = virialpress_pxy(ibin) + pxy
          virialpress_pxz(ibin) = virialpress_pxz(ibin) + pxz

          virialpress_pyx(ibin) = virialpress_pyx(ibin) + pyx
          virialpress_pyy(ibin) = virialpress_pyy(ibin) + pyy
          virialpress_pyz(ibin) = virialpress_pyz(ibin) + pyz

          virialpress_pzx(ibin) = virialpress_pzx(ibin) + pzx
          virialpress_pzy(ibin) = virialpress_pzy(ibin) + pzy
          virialpress_pzz(ibin) = virialpress_pzz(ibin) + pzz

        !End bin cycle
        ENDDO


      ! if both imol and jmol are in the central box 
      Else
        ! == Irving-Kirkwood definition ==
        !Loop through the bins on the POSITIVE x side
        DO ibin= 1, n_bins
        
          ! Calculate intersection with planes ibin
          solk = (bin_posp(ibin) - rxi)/xij

          if((solk .LT. 0.0d0) .OR. (solk .GT. 1.0d0)) CYCLE

          virialpress_pxx(ibin) = virialpress_pxx(ibin) + pxx
          virialpress_pxy(ibin) = virialpress_pxy(ibin) + pxy
          virialpress_pxz(ibin) = virialpress_pxz(ibin) + pxz

          virialpress_pyx(ibin) = virialpress_pyx(ibin) + pyx
          virialpress_pyy(ibin) = virialpress_pyy(ibin) + pyy
          virialpress_pyz(ibin) = virialpress_pyz(ibin) + pyz

          virialpress_pzx(ibin) = virialpress_pzx(ibin) + pzx
          virialpress_pzy(ibin) = virialpress_pzy(ibin) + pzy
          virialpress_pzz(ibin) = virialpress_pzz(ibin) + pzz

        !End bin cycle
        ENDDO

        !Loop through the bins on the NEGATIVE x side
        DO ibin= 1, n_bins
        
          ! Calculate intersection with planes ibin
          solk = (bin_posn(ibin) - rxi)/xij

          if((solk .LT. 0.0d0) .OR. (solk .GT. 1.0d0)) CYCLE

          virialpress_pxx(ibin) = virialpress_pxx(ibin) + pxx
          virialpress_pxy(ibin) = virialpress_pxy(ibin) + pxy
          virialpress_pxz(ibin) = virialpress_pxz(ibin) + pxz

          virialpress_pyx(ibin) = virialpress_pyx(ibin) + pyx
          virialpress_pyy(ibin) = virialpress_pyy(ibin) + pyy
          virialpress_pyz(ibin) = virialpress_pyz(ibin) + pyz

          virialpress_pzx(ibin) = virialpress_pzx(ibin) + pzx
          virialpress_pzy(ibin) = virialpress_pzy(ibin) + pzy
          virialpress_pzz(ibin) = virialpress_pzz(ibin) + pzz

        !End bin cycle
        ENDDO

      ! End PBC treatment
      End If


    ! End loop over jmol
    ENDDO

  !End loop over imol 
  ENDDO  


  ! write progress
  If(MOD(iframe,1) .eq. 0) Then
    Call CPU_TIME(cpu_nd)
    Open(5,File='progress',action='write', status='replace')
    Write(5,'(A,I10,F15.7,F15.7,F15.7)') 'Finished sampling of frame #', st_frame+iframe-1, boxsz(1,iframe), &
      & boxsz(2,iframe), boxsz(3,iframe)
    Write(5,'(A,F7.3,A,2X,F7.3,A)') 'Time left:', (cpu_nd - cpu_st)/iframe*(ns_frame-iframe)/3600.0, 'hours', &
      & (cpu_nd - cpu_st)/iframe*(ns_frame-iframe)/3600.0/24.0, 'days'
    FLUSH(5)
    Close(5)
  ENDIF

! End loop over frames
End do

! ------------- Average statistics and write statistics to file ------------------!
! --------------- x-density ---------------!
! Open file
OPEN(3,FILE='x-density.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Loop over bins
Do ibin = 1, n_bins
  
  denavg(ibin) = denavg(ibin)/DBLE(ns_frame)

End Do

! Write data to file
! loop over molecule types
! Write molecule info
!Write(3,'(2A)') 'Molecule name: ', site_name
Write(3,*) ' x             x-rho [1/A^3]'

! Loop over bins
Do ibin = 1, n_bins

  ! Write to file in units of [1/A^3]
  Write(3,'(F8.4,8X,E15.7)')  bin_posp(ibin), denavg(ibin)
                                        
End Do

! Close file
CLOSE(3)


! ------------ planar pressure tensor ----------------!
Do ibin = 1, n_bins

  ! Kinetic part in unit of [K/A^3]
  prkin(ibin) = denavg(ibin)*temp

  ! Configurational part in unit of [K/A^3]
  ! average over two interfaces and over total number of frames
  virialpxxavg(ibin) = virialpress_pxx(ibin)/(2.0d0 * DBLE(ns_frame)) 
  virialpxyavg(ibin) = virialpress_pxy(ibin)/(2.0d0 * DBLE(ns_frame)) 
  virialpxzavg(ibin) = virialpress_pxz(ibin)/(2.0d0 * DBLE(ns_frame)) 

  virialpyxavg(ibin) = virialpress_pyx(ibin)/(2.0d0 * DBLE(ns_frame)) 
  virialpyyavg(ibin) = virialpress_pyy(ibin)/(2.0d0 * DBLE(ns_frame)) 
  virialpyzavg(ibin) = virialpress_pyz(ibin)/(2.0d0 * DBLE(ns_frame)) 

  virialpzxavg(ibin) = virialpress_pzx(ibin)/(2.0d0 * DBLE(ns_frame)) 
  virialpzyavg(ibin) = virialpress_pzy(ibin)/(2.0d0 * DBLE(ns_frame)) 
  virialpzzavg(ibin) = virialpress_pzz(ibin)/(2.0d0 * DBLE(ns_frame)) 


End Do

! Open file
OPEN(4,FILE='press_plan.txt', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Write file head
Write(4,*) 'Planar pressure tensor using IK route'
Write(4,*) 'Unit: pressure in [bar] and length in [Angstrom]'
Write(4,'(A)') 'x Pkin Pxx Pxy Pxz Pyx Pyy Pyz Pzx Pzy Pzz'

! Write data to file
! Loop over bins
Do ibin = 1, n_bins

  ! Write to file
  Write(4,'(F8.4,F16.4,F16.4,F16.4,F16.4,F16.4,F16.4,F16.4,F16.4,F16.4,F16.4)')  &
  & bin_posp(ibin), &
  & prkin(ibin)*PCOEFF, &
  & (virialpxxavg(ibin)+prkin(ibin))*PCOEFF, &
  & virialpxyavg(ibin)*PCOEFF, &
  & virialpxzavg(ibin)*PCOEFF, &

  & virialpyxavg(ibin)*PCOEFF, &
  & (virialpyyavg(ibin)+prkin(ibin))*PCOEFF, &
  & virialpyzavg(ibin)*PCOEFF, &

  & virialpzxavg(ibin)*PCOEFF, &
  & virialpzyavg(ibin)*PCOEFF, &
  & (virialpzzavg(ibin)+prkin(ibin))*PCOEFF

End Do

! Close file
CLOSE(4)




! ---------- Deallocate variables to free space -------------!
Deallocate(boxsz)
Deallocate(ox)
Deallocate(oy)
Deallocate(oz)
Deallocate(rx_s)
Deallocate(ry_s)
Deallocate(rz_s)
Deallocate(rx)
Deallocate(ry)
Deallocate(rz)

Deallocate(virialpress_pxx) 
Deallocate(virialpress_pxy) 
Deallocate(virialpress_pxz) 
Deallocate(virialpress_pyx) 
Deallocate(virialpress_pyy) 
Deallocate(virialpress_pyz) 
Deallocate(virialpress_pzx) 
Deallocate(virialpress_pzy) 
Deallocate(virialpress_pzz) 
Deallocate(virialpxxavg) 
Deallocate(virialpxyavg) 
Deallocate(virialpxzavg) 
Deallocate(virialpyxavg) 
Deallocate(virialpyyavg) 
Deallocate(virialpyzavg) 
Deallocate(virialpzxavg) 
Deallocate(virialpzyavg) 
Deallocate(virialpzzavg) 
Deallocate(prkin) 

Deallocate(denavg) 
Deallocate(denblk) 





END PROGRAM virialpress_planar






