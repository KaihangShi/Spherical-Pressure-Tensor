! ==========================================================================
! This version is for calculating IK contour local pressure tensor
!  in spherical coordinates
!  Input file is xtc format from MD (GROMACS) simulation
! --------------------------------------------------------------------------
! Note: 
! 1). Periodic boundary condition has been applied for pressure calculation
! 2). Hard sphere models are used here
! 
! Author: Kaihang Shi
! Last update: Feb 27, 2020
! ==========================================================================

PROGRAM virialpress_sphere

! use an external module for analyzing the GROMACS trajectory
use gmxfort_trajectory

IMPLICIT NONE

! ---------------- Input Parameters ---------------------------!
! Control parameters 
Integer, Parameter :: st_frame = 1
Integer, Parameter :: nd_frame = 1
! xtc file PATH 
Character*256, Parameter :: xtcpath = "/gpfs_backup/gubbins_data/kshi3/sphe_pressure_tensor/iv/traj_IV.xtc"  
! solid cluster center-of-mass file
Character*256, Parameter :: compath = "/gpfs_backup/gubbins_data/kshi3/sphe_pressure_tensor/iv/CM_true_IV.dat"
! System parameter 
Double Precision, Parameter :: temp = 179.805d0      ! [Kelvin]
Integer, Parameter :: n_mol_types = 1               ! number of molecule (hard sphere) types  
Integer, Parameter :: n_sites_tot = 48207            ! total number of hard spheres in the system
! Intermolecular parameter 
Double Precision, Parameter :: sigma = 3.405d0      ! [Angstrom] 
Double Precision, Parameter :: epsilonkb = 119.87d0  ! [Kelvin]
Double Precision, Parameter :: r_cut = 50.0/49.0*sigma    ! cutoff radius for hard sphere model (Jover, 2012 JCP) [Angstrom]
!Double Precision, Parameter :: mol_mass = 39.948d0  ! [g/mol]
! Parameters for radial density and pressure calculation
Double Precision, Parameter :: rden_cut = 62.5d0   ! Cutoff distance for density/pressure calculation [Angstrom]
Integer, Parameter :: rden_bins = 300             ! density profile resolution

! ---------------- Define Variables -------------------------!
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
Double Precision, Dimension(:,:), Allocatable :: rx_s, ry_s, rz_s
! Statistics
! Only fluid-fluid contribution is considered here
Double Precision, Dimension(:), Allocatable :: virialpress_pnr
Double Precision, Dimension(:), Allocatable :: virialpress_ptt
Double Precision, Dimension(:), Allocatable :: virialpress_ptp
Double Precision, Dimension(:), Allocatable :: virialpnravg
Double Precision, Dimension(:), Allocatable :: virialpttavg
Double Precision, Dimension(:), Allocatable :: virialptpavg
Double Precision, Dimension(:), Allocatable :: prkin
! r-density statistics (Density profile of each molecule type in the cylindrical system)
! Cutoff radius for density and cylindrical pressure tensor calculation
! limit in radial direction
Double Precision :: rden_lim
! r-density length in r direction for each bin
Double Precision :: rden_dr, rden_drsq, rden_drcb
Double Precision, Dimension(rden_bins) :: rden_dvol, posr
! Statistics 
Double Precision, Dimension(:,:), Allocatable :: rdenavg
Double Precision, Dimension(:,:,:), Allocatable :: rdenblk

!Variables for pressure calculation
INTEGER :: imol, jmol, ibin, iframe, idirec, ierr, isol
INTEGER :: isite, jsite, itype, jtype, isitetype, jsitetype
!DOUBLE PRECISION :: rxi, ryi, rzi, rij, rxj, ryj, rzj
DOUBLE PRECISION :: rxis, ryis, rzis, rxjs, ryjs, rzjs, rxiis, ryiis, rziis
DOUBLE PRECISION :: sr, sr25, phitz
DOUBLE PRECISION :: pnrc, pttc, dudr, ptpc
DOUBLE PRECISION :: rijs, xijs, yijs, zijs
DOUBLE PRECISION :: xijssq, yijssq, zijssq
Double Precision :: radistsq_i, radistsq_j, radist_ij
Double Precision, Dimension(:), Allocatable :: radist, radistsq
Double Precision :: aa, bb, cc, dcrt
Double Precision, Dimension(2) :: alpk
Double Precision :: rijer, rijetsq, rijepsq, dxysq, xl, yl, zl
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
Allocate(rx_s(n_sites_tot,ns_frame), STAT=ierr)
Allocate(ry_s(n_sites_tot,ns_frame), STAT=ierr)
Allocate(rz_s(n_sites_tot,ns_frame), STAT=ierr)
Allocate(virialpress_pnr(rden_bins), STAT=ierr) 
Allocate(virialpress_ptt(rden_bins), STAT=ierr) 
Allocate(virialpress_ptp(rden_bins), STAT=ierr) 
Allocate(virialpnravg(rden_bins), STAT=ierr) 
Allocate(virialpttavg(rden_bins), STAT=ierr) 
Allocate(virialptpavg(rden_bins), STAT=ierr) 
Allocate(prkin(rden_bins), STAT=ierr) 
Allocate(rdenavg(rden_bins,n_mol_types), STAT=ierr) 
Allocate(rdenblk(rden_bins,n_mol_types,ns_frame), STAT=ierr) 
Allocate(radist(n_sites_tot), STAT=ierr)
Allocate(radistsq(n_sites_tot), STAT=ierr)

! -- Initialize parameter -----!
virialpress_pnr = 0.0d0
virialpress_ptt = 0.0d0
virialpress_ptp = 0.0d0
virialpnravg = 0.0d0
virialpttavg = 0.0d0
virialptpavg = 0.0d0
prkin = 0.0d0
rdenavg = 0.0d0
rdenblk = 0.0d0
rden_dr = rden_cut/DBLE(rden_bins)
rden_drsq = rden_dr**2
rden_drcb = rden_dr**3
rden_lim = rden_cut + r_cut
! Calculate constants in spherical geometry
Do ibin = 1, rden_bins
    ! Calculate shperical shell volume
    rden_dvol(ibin) = 4.0/3.0*Pi*rden_drcb*DBLE(3*ibin**2-3*ibin+1)
    ! Get r-distance of ibin 
    posr(ibin) = (DBLE(ibin)-0.5d0)*rden_dr
End Do



! --------- Write header ------------!
Write(*,'(A)') '==========================================================  '
Write(*,'(A)') ' This program is used to calculate the pressure tensor for  ' 
Write(*,'(A)') ' spherical geometry from virial route                     '
Write(*,'(A)') ' using the Irving-Kirkwood (IK) definition of contour      '
Write(*,'(A)') '                                                            '  
Write(*,'(A)') ' This version is for continueous hard sphere model          '    
Write(*,'(A)') ' ========================================================== '
Write(*,'(A)') '                                                            '

! Creat a file recording running progress
Open(5,File='progress',action='write', status='replace')
! Open solid cluster COM file to read
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

    ! Read in solid cluster COM (in nanometers)
    Read(7,*) ox(c_frame), oy(c_frame), oz(c_frame)

    ! Loop over all hs particles
    Do isite = 1, n_sites_tot

        ! Get atom's coordinates in [nanometer]
        myatom = trj%x(1,isite)
        ! Reassign coordinates using center of box as the origin and convert units to [Angstrom]
        rx_s(isite,c_frame) = (myatom(1) - ox(c_frame))*10.0
        ry_s(isite,c_frame) = (myatom(2) - oy(c_frame))*10.0
        rz_s(isite,c_frame) = (myatom(3) - oz(c_frame))*10.0

        ! Apply minimum image convention
        rx_s(isite,c_frame) = rx_s(isite,c_frame) - dNINT(rx_s(isite,c_frame)/boxsz(1,c_frame))*boxsz(1,c_frame)
        ry_s(isite,c_frame) = ry_s(isite,c_frame) - dNINT(ry_s(isite,c_frame)/boxsz(2,c_frame))*boxsz(2,c_frame)
        rz_s(isite,c_frame) = rz_s(isite,c_frame) - dNINT(rz_s(isite,c_frame)/boxsz(3,c_frame))*boxsz(3,c_frame)

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

! ------------ TEST coordiantes processing
! OPEN(10,FILE='testatom.XYZ', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
! Write(10,*) n_sites_tot
! iframe = 1
! Write(10,*) (boxsz(idirec,iframe), idirec=1,3), ox(iframe), oy(iframe), oz(iframe)
! Do isite = 1, n_sites_tot
!     Write (10,'(A,F15.7,F15.7,F15.7)') 'Ar', rx_s(isite,iframe), ry_s(isite,iframe), rz_s(isite,iframe)
! End do
! ! Close file
! CLOSE(10)
! STOP


Call CPU_TIME(cpu_st)
! --------------- Start postprocessing ----------------------!
Do iframe = 1, ns_frame
    
  ! ----------- Sampling r-density ------------!
  ! Loop over all sites/molecules
  Do isite = 1, n_sites_tot

    ! Calculate r-distance of isite in spherical coordiantes
    radistsq(isite) = rx_s(isite,iframe)**2 + ry_s(isite,iframe)**2 + rz_s(isite,iframe)**2
    radist(isite) = dSQRT(radistsq(isite))

    If (radist(isite) .LT. rden_cut) Then
      ! Calculate ibin number
      ibin = FLOOR(radist(isite)/rden_dr) + 1

      ! Accumulate number 
      rdenblk(ibin,1,iframe) = rdenblk(ibin,1,iframe) + 1.0d0
    End If

  ! End loop over all molecules 
  End Do

  ! Loop over bins
  Do ibin = 1, rden_bins
    ! Loop over molecule types
    Do itype = 1, n_mol_types

      ! Convert to number density (1/A^3)
      rdenavg(ibin,itype) = rdenavg(ibin,itype) + rdenblk(ibin,itype,iframe)/rden_dvol(ibin)

    End Do          
  End Do  
    

  ! ----- Sampling spherical pressure tensor -------!   
  ! Loop over sites/molecules
  DO isite = 1, n_sites_tot

    rxis = rx_s(isite,iframe) 
    ryis = ry_s(isite,iframe) 
    rzis = rz_s(isite,iframe)

    ! Determine if carrying on (rden_lim = rden_cut + r_cut)
    If (radist(isite) .GT. rden_lim) CYCLE

    !Fluid-fluid interaction, avoiding double count of i,j and j,i so just
    !loop by using i<j     
    IF (isite .LT. n_sites_tot) THEN
      DO jsite = isite+1, n_sites_tot

        rxjs = rx_s(jsite,iframe) 
        ryjs = ry_s(jsite,iframe) 
        rzjs = rz_s(jsite,iframe)

        ! Determine if carrying on
        If (radist(jsite) .GT. rden_lim) CYCLE

        !Calculate vector between sites
        xijs = rxjs - rxis
        yijs = ryjs - ryis
        zijs = rzjs - rzis

        ! Apply minimum image convention
        xijs = xijs - dNINT(xijs/boxsz(1,iframe))*boxsz(1,iframe)
        yijs = yijs - dNINT(yijs/boxsz(2,iframe))*boxsz(2,iframe)
        zijs = zijs - dNINT(zijs/boxsz(3,iframe))*boxsz(3,iframe)

        !Square the values
        xijssq=xijs*xijs
        yijssq=yijs*yijs
        zijssq=zijs*zijs

        !Calculate distance between i and j sites
        rijs=dSQRT(xijssq+yijssq+zijssq)

        ! Check with force cutoff
        If (rijs .GE. r_cut) CYCLE

        !Calculate 50-49 continuous hard sphere force (based on sites)
        sr = sigma/rijs
        sr25 = sr**25
        phitz = 134.552662342121d0*epsilonkb/rijs
        dudr = phitz*(49.0*sr25*sr25/sr - 50.0*sr25*sr25)


        ! Special treatment of PBC
        If((xijs .NE. (rxjs -rxis)) .or. (yijs .NE. (ryjs -ryis)) .or. (zijs .NE. (rzjs-rzis))) Then

            ! -- Assume i site in the central box and ij pair contribue once
            ! Calculate coefficient in ax^2+bx+c form for R=R_l use
            aa = xijssq + yijssq + zijssq
            bb = 2.0*(rxis*xijs + ryis*yijs + rzis*zijs)

            ! == Irving-Kirkwood definition ==
            !Loop through the bins
            DO ibin= 1, rden_bins

                ! Calculate alpha_k
                cc = radistsq(isite) - posr(ibin)**2
                ! Discriminant 
                dcrt = bb**2 - 4.0*aa*cc
                if(dcrt .LE. 0.0d0) CYCLE
                ! two solutions for alpha
                alpk(1) = (-bb + dSQRT(dcrt))/(2.0*aa)
                alpk(2) = (-bb - dSQRT(dcrt))/(2.0*aa)

                if((alpk(1) .LT. 0.0d0 ) .and. (alpk(2) .GT. 1.0d0)) CYCLE
                if((alpk(2) .LT. 0.0d0 ) .and. (alpk(1) .GT. 1.0d0)) CYCLE

                ! Loop over solutions
                Do isol = 1,2 

                    ! contribute
                    if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then

                        ! - Normal pressure component
                        ! Calculate basic parameters
                        xl = rxis + alpk(isol)*xijs
                        yl = ryis + alpk(isol)*yijs
                        zl = rzis + alpk(isol)*zijs    
                        rijer = dABS(xijs * xl/posr(ibin) + yijs * yl/posr(ibin) + zijs * zl/posr(ibin))

                        ! Update normal pressure in radial direction
                        pnrc = rijer * dudr/rijs 
                        virialpress_pnr(ibin) = virialpress_pnr(ibin) + pnrc

                        ! - Tangential pressure in polar direction (theta)
                        ! Calculate basic parameters
                        dxysq = xl**2 + yl**2
                        rijetsq = (xl*(rzis*xijs - rxis*zijs) + yl*(rzis*yijs - zijs*ryis))**2/(posr(ibin)**2 * dxysq)

                        ! update tangential pressure in theta-theta direction
                        pttc = rijetsq/rijer * dudr/rijs 
                        virialpress_ptt(ibin) = virialpress_ptt(ibin) + pttc

                        ! - Tangential pressure in azimuthal direction (phi)
                        ! Calculate basic parameters
                        rijepsq = (xl*yijs - yl*xijs)**2/dxysq

                        ! update tangential pressure in phi-phi direction
                        ptpc = rijepsq/rijer * dudr/rijs 
                        virialpress_ptp(ibin) = virialpress_ptp(ibin) + ptpc


                    End if
                End Do

            !End bin cycle
            ENDDO


            ! -- assume j site in the central box and ij pair contribute again
            rxiis = rxjs - xijs
            ryiis = ryjs - yijs
            rziis = rzjs - zijs

            ! Calculate coefficient in ax^2+bx+c form for R=R_l use
            aa = xijssq + yijssq + zijssq
            bb = 2.0*(rxiis*xijs + ryiis*yijs + rziis*zijs)

            ! == Irving-Kirkwood definition ==
            !Loop through the bins
            DO ibin= 1, rden_bins

                ! Calculate alpha_k
                cc = rxiis**2 + ryiis**2 + rziis**2 - posr(ibin)**2
                ! Discriminant 
                dcrt = bb**2 - 4.0*aa*cc
                if(dcrt .LE. 0.0d0) CYCLE
                ! two solutions for alpha
                alpk(1) = (-bb + dSQRT(dcrt))/(2.0*aa)
                alpk(2) = (-bb - dSQRT(dcrt))/(2.0*aa)

                if((alpk(1) .LT. 0.0d0 ) .and. (alpk(2) .GT. 1.0d0)) CYCLE
                if((alpk(2) .LT. 0.0d0 ) .and. (alpk(1) .GT. 1.0d0)) CYCLE

                ! Loop over solutions
                Do isol = 1,2 

                    ! contribute
                    if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then

                        ! - Normal pressure component
                        ! Calculate basic parameters
                        xl = rxiis + alpk(isol)*xijs
                        yl = ryiis + alpk(isol)*yijs
                        zl = rziis + alpk(isol)*zijs    
                        rijer = dABS(xijs * xl/posr(ibin) + yijs * yl/posr(ibin) + zijs * zl/posr(ibin))

                        ! Update normal pressure in radial direction
                        pnrc = rijer * dudr/rijs 
                        virialpress_pnr(ibin) = virialpress_pnr(ibin) + pnrc

                        ! - Tangential pressure in polar direction (theta)
                        ! Calculate basic parameters
                        dxysq = xl**2 + yl**2
                        rijetsq = (xl*(rziis*xijs - rxiis*zijs) + yl*(rziis*yijs - zijs*ryiis))**2/(posr(ibin)**2 * dxysq)

                        ! update tangential pressure in theta-theta direction
                        pttc = rijetsq/rijer * dudr/rijs 
                        virialpress_ptt(ibin) = virialpress_ptt(ibin) + pttc

                        ! - Tangential pressure in azimuthal direction (phi)
                        ! Calculate basic parameters
                        rijepsq = (xl*yijs - yl*xijs)**2/dxysq

                        ! update tangential pressure in phi-phi direction
                        ptpc = rijepsq/rijer * dudr/rijs 
                        virialpress_ptp(ibin) = virialpress_ptp(ibin) + ptpc


                    End if
                End Do

            !End bin cycle
            ENDDO


        ! if both particle i and particle j are in the central box 
        Else

        	! Calculate coefficient in ax^2+bx+c form for R=R_l use
            aa = xijssq + yijssq + zijssq
            bb = 2.0*(rxis*xijs + ryis*yijs + rzis*zijs)

            ! == Irving-Kirkwood definition ==
            !Loop through the bins
            DO ibin= 1, rden_bins

                ! Calculate alpha_k
                cc = radistsq(isite) - posr(ibin)**2
                ! Discriminant 
                dcrt = bb**2 - 4.0*aa*cc
                if(dcrt .LE. 0.0d0) CYCLE
                ! two solutions for alpha
                alpk(1) = (-bb + dSQRT(dcrt))/(2.0*aa)
                alpk(2) = (-bb - dSQRT(dcrt))/(2.0*aa)

                if((alpk(1) .LT. 0.0d0 ) .and. (alpk(2) .GT. 1.0d0)) CYCLE
                if((alpk(2) .LT. 0.0d0 ) .and. (alpk(1) .GT. 1.0d0)) CYCLE

                ! Loop over solutions
                Do isol = 1,2 

                    ! contribute
                    if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then

                        ! - Normal pressure component
                        ! Calculate basic parameters
                        xl = rxis + alpk(isol)*xijs
                        yl = ryis + alpk(isol)*yijs
                        zl = rzis + alpk(isol)*zijs    
                        rijer = dABS(xijs * xl/posr(ibin) + yijs * yl/posr(ibin) + zijs * zl/posr(ibin))

                        ! Update normal pressure in radial direction
                        pnrc = rijer * dudr/rijs
                        virialpress_pnr(ibin) = virialpress_pnr(ibin) + pnrc

                        ! - Tangential pressure in polar direction (theta)
                        ! Calculate basic parameters
                        dxysq = xl**2 + yl**2
                        rijetsq = (xl*(rzis*xijs - rxis*zijs) + yl*(rzis*yijs - zijs*ryis))**2/(posr(ibin)**2 * dxysq)

                        ! update tangential pressure in theta-theta direction
                        pttc = rijetsq/rijer * dudr/rijs 
                        virialpress_ptt(ibin) = virialpress_ptt(ibin) + pttc

                        ! - Tangential pressure in azimuthal direction (phi)
                        ! Calculate basic parameters
                        rijepsq = (xl*yijs - yl*xijs)**2/dxysq

                        ! update tangential pressure in phi-phi direction
                        ptpc = rijepsq/rijer * dudr/rijs 
                        virialpress_ptp(ibin) = virialpress_ptp(ibin) + ptpc


                    End if
                End Do

            !End bin cycle
            ENDDO

           
        ! End PBC treatment
        End If

      !End loop over jmol sites 
      ENDDO       
    !End IF imol < n_mol_tot  
    ENDIF    
  !End loop over imol sites
  ENDDO  


  ! write progress
  If(MOD(iframe,10) .eq. 0) Then
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
! --------------- Radial density ---------------!
! Open file
OPEN(2,FILE='r-density.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Loop over bins
Do ibin = 1, rden_bins
    ! Loop over molecule types
    Do itype = 1, n_mol_types

        rdenavg(ibin,itype) = rdenavg(ibin,itype)/DBLE(ns_frame)
        ! Convert unit from [1/A^3] to [g/ml]
        !rdenavg(ibin,itype) = (mol_mass/(Na*1.0d-24))*rdenavg(ibin,itype)

    End Do
End Do

! Write data to file
! loop over molecule types
Do itype = 1, n_mol_types
  	! Write molecule info
    !Write(2,'(2A)') 'Molecule name: ', site_name
    Write(2,*) ' R             R-rho [1/A^3]'

    ! Loop over bins
    Do ibin = 1, rden_bins

        ! Write to file in units of [1/A^3]
        Write(2,'(F8.4,8X,E15.7)')  posr(ibin), rdenavg(ibin,itype)
                                            
    End Do
End Do

! Close file
CLOSE(2)

! ------------ Spherical pressure tensor ----------------!
Do ibin = 1, rden_bins

  ! Pressure (Kinetic part) in unit of [K/A^3]
  prkin(ibin) = rdenavg(ibin,1)*temp

  ! Irving-Kirkwood route
  virialpnravg(ibin) = virialpress_pnr(ibin)/DBLE(ns_frame)
  virialpttavg(ibin) = virialpress_ptt(ibin)/DBLE(ns_frame)
  virialptpavg(ibin) = virialpress_ptp(ibin)/DBLE(ns_frame)

  ! Calculate final pressure (configurational part) in unit of [K/A^3]
  virialpnravg(ibin) = -1.0d0/(2.0*two_Pi*posr(ibin)**2)*virialpnravg(ibin)
  virialpttavg(ibin) = -1.0d0/(2.0*two_Pi*posr(ibin)**2)*virialpttavg(ibin)
  virialptpavg(ibin) = -1.0d0/(2.0*two_Pi*posr(ibin)**2)*virialptpavg(ibin)


End Do

! Open file
OPEN(3,FILE='press_spheIK.txt', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Write file head
Write(3,*) 'Spherical pressure tensor from virial route using Irving-Kirkwood definition'
Write(3,*) 'Unit: pressure in [bar] and length in [Angstrom]'
Write(3,'(A)') &
                & '  R               Pkin          Pnr(ff)             Pnr(tot)          Ptt(ff)&
                &         Ptt(tot)         Ptp(ff)          Ptp(tot)'

! Write data to file
! Loop over bins
Do ibin = 1, rden_bins

    ! Write to file
    Write(3,'(F8.4,F16.4,F16.4,3X,F16.4,3X,F16.4,F16.4,F16.4,F16.4)')  &
        & posr(ibin), prkin(ibin)*PCOEFF, &
        & virialpnravg(ibin)*PCOEFF, &
        & (virialpnravg(ibin)+prkin(ibin))*PCOEFF, &
        & virialpttavg(ibin)*PCOEFF, &
        & (virialpttavg(ibin)+prkin(ibin))*PCOEFF, &
        & virialptpavg(ibin)*PCOEFF, &
        & (virialptpavg(ibin)+prkin(ibin))*PCOEFF

End Do

! Close file
CLOSE(3)




! ---------- Deallocate variables to free space -------------!
Deallocate(boxsz)
Deallocate(ox)
Deallocate(oy)
Deallocate(oz)
Deallocate(rx_s)
Deallocate(ry_s)
Deallocate(rz_s)
Deallocate(virialpress_pnr) 
Deallocate(virialpress_ptt) 
Deallocate(virialpress_ptp) 
Deallocate(virialpnravg) 
Deallocate(virialpttavg) 
Deallocate(virialptpavg) 
Deallocate(prkin) 
Deallocate(rdenavg) 
Deallocate(rdenblk) 
Deallocate(radist)
Deallocate(radistsq)



END PROGRAM virialpress_sphere






