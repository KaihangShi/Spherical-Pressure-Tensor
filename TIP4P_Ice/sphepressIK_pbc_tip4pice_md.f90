! ======================================================================================================
! Calculation of spherical pressure tensor for pure rigid molecules (eg TIP4P/ICE):
! 	- Irving-Kirkwood contour is used
! 	- Periodic boundary conditions are applied
!   - Wolf potential (damped and cut) is applied to reproduce the Ewald method used in simulation
!   - Input file is xtc format from MD (GROMACS) simulation
!  
! Author: Kaihang Shi
! Email: kaihangshi0@gmail.com
! Last update: March 30, 2021
! ======================================================================================================

PROGRAM virialpress_sphere

! use an external module for analyzing the GROMACS trajectory
use gmxfort_trajectory

IMPLICIT NONE


! ---------------- Input Parameters ---------------------------!
! Control parameters (starting/ending frame No. for sampling)
Integer, Parameter :: st_frame = 1
Integer, Parameter :: nd_frame = 1

! xtc file PATH 
Character*256, Parameter :: xtcpath = "/share/gubbins/kshi3/traj.xtc"  
! solid cluster center-of-mass file (evaluated using oxygen atom only)
Character*256, Parameter :: compath = "/gpfs_backup/gubbins_data/kshi3/sphe_pressure_tensor/tip4p_ice/CM_all.dat"

! System parameter 
Double Precision, Parameter :: temp = 255.0d0      ! system temperature [Kelvin]
Integer, Parameter :: n_mol_tot = 78856              ! total number of molecules in the system
Integer, Parameter :: n_sites_tot = 315424           ! Total number of atoms in the system
Integer, Parameter :: n_mol_sites = 4               ! number of sites in the water molecule    

! Intermolecular parameter 
Double Precision, Parameter :: sigma_o = 3.1668d0      ! sigma for oxygen atom [Angstrom] 
Double Precision, Parameter :: epsilonkb_o = 106.1d0    ! epsilon/k_B for oxygen atom [Kelvin]
Double Precision, Parameter :: q_m = -1.1794d0        ! point charge for M-site, TIP4P/ICE [e]
Double Precision, Parameter :: q_h = 0.5897d0        ! point charge for H, TIP4P/ICE model [e]
Double Precision, Parameter :: r_ljcut = 9.0d0      ! LJ cutoff radius [Angstrom]
Double Precision, Parameter :: alp = 0.04d0       ! Alpha Parameter for Wolf potential method [1/Angstrom]
Double Precision, Parameter :: r_coulcut = 65.0d0     ! cutoff radius for Wolf potential method [Angstrom]

! Mass 
Double Precision, Parameter :: mol_mass = 18.01528d0  ! molecular weight of water [g/mol]
Double Precision, Parameter :: mass_o = 15.999400
Double Precision, Parameter :: mass_h = 1.007940

! Parameters for radial density and pressure calculation
Double Precision, Parameter :: rden_cut = 10.5d0   ! Cutoff distance for density/pressure calculation [Angstrom]
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
Double Precision, Dimension(:,:,:), Allocatable :: rx_s, ry_s, rz_s
Double Precision, Dimension(:,:), Allocatable :: rx, ry, rz
Character*8, Dimension(n_mol_sites) :: atom_name

! Statistics
! Only fluid-fluid contribution is considered here
Double Precision, Dimension(:), Allocatable :: virialpress_pnr
Double Precision, Dimension(:), Allocatable :: virialpress_ptt
Double Precision, Dimension(:), Allocatable :: virialpress_ptp
Double Precision, Dimension(:), Allocatable :: virialpnravg
Double Precision, Dimension(:), Allocatable :: virialpttavg
Double Precision, Dimension(:), Allocatable :: virialptpavg
Double Precision, Dimension(:), Allocatable :: prkin


! r-density statistics 
! r-density length in r direction for each bin
Double Precision :: rden_dr, rden_drsq, rden_drcb
Double Precision, Dimension(rden_bins) :: rden_dvol, posr
! Statistics 
Double Precision, Dimension(:), Allocatable :: rdenavg
Double Precision, Dimension(:,:), Allocatable :: rdenblk


!Variables for pressure calculation
INTEGER :: imol, jmol, ibin, iframe, ierr, isol
INTEGER :: isite, jsite
DOUBLE PRECISION :: rxi, ryi, rzi, rij, rijsq, rxj, ryj, rzj, xij, yij, zij, rxii, ryii, rzii
DOUBLE PRECISION :: rxis, ryis, rzis, rxjs, ryjs, rzjs
DOUBLE PRECISION :: rijs, rijssq,xijs, yijs, zijs
DOUBLE PRECISION :: xijssq, yijssq, zijssq, xijsq, yijsq, zijsq

DOUBLE PRECISION :: sr3, sr6, sr12, phitz, alppi, alpsq, f_lj, f_coul, f_coul_rc,f_tot
DOUBLE PRECISION :: pnrc, pttc, ptpc

Double Precision, Dimension(:), Allocatable :: radist, radistsq
Double Precision :: aa, bb, cc, dcrt
Double Precision, Dimension(2) :: alpk

Double Precision :: rijer, rijser, xl, yl, zl, dxy, rijet, rijset, rijep, rijsep

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


Allocate(virialpress_pnr(rden_bins), STAT=ierr) 
Allocate(virialpress_ptt(rden_bins), STAT=ierr) 
Allocate(virialpress_ptp(rden_bins), STAT=ierr) 
Allocate(virialpnravg(rden_bins), STAT=ierr) 
Allocate(virialpttavg(rden_bins), STAT=ierr) 
Allocate(virialptpavg(rden_bins), STAT=ierr) 
Allocate(prkin(rden_bins), STAT=ierr) 

Allocate(rdenavg(rden_bins), STAT=ierr) 
Allocate(rdenblk(rden_bins,ns_frame), STAT=ierr) 

Allocate(radist(n_mol_tot), STAT=ierr)
Allocate(radistsq(n_mol_tot), STAT=ierr)



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

! Calculate constants in spherical geometry
Do ibin = 1, rden_bins
	! Calculate shperical shell volume
	rden_dvol(ibin) = 4.0/3.0*Pi*rden_drcb*DBLE(3*ibin**2-3*ibin+1)
	! Get r-distance of ibin 
	posr(ibin) = (DBLE(ibin)-0.5d0)*rden_dr
End Do

! Calculate parameters in shifted-force Wolf potential
alpsq = alp**2
alppi = 2.0d0*alp/dSQRT(Pi)
f_coul_rc = dERFC(alp*r_coulcut)/r_coulcut**2 + alppi*dEXP(-alpsq*r_coulcut**2)/r_coulcut



! --------- Write header ------------!
Write(*,'(A)') '==========================================================  '
Write(*,'(A)') ' Calculation of the molecular pressure tensor in         ' 
Write(*,'(A)') ' spherical geometry from virial route                     '
Write(*,'(A)') ' using the Irving-Kirkwood (IK) definition of contour      '
Write(*,'(A)') '                                                            '  
Write(*,'(A)') ' This version is for TIP4P/ice model                       '    
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

		! the water bond is 1 A, to make sure both molecular COM and atom positon satisfy the minimum image convention
		If ((r_coulcut+2.0d0) .GT. Min(boxsz(1,c_frame),boxsz(2,c_frame),boxsz(3,c_frame))/2.0d0) Then
			Write(5,*) 'Cutoff is too large for this box size!'
			STOP
		End If

		! Read in solid cluster COM (in nanometers)
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

! ! TEST coordiantes processing
! OPEN(2,FILE='testatom.XYZ', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
! !OPEN(3,FILE='testmol.XYZ', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
! atom_name(1) = 'O'
! atom_name(2) = 'H'
! atom_name(3) = 'H'
! atom_name(4) = 'M'
! Write(2,*) natom
! Do iframe = 1, 1
! 	Write(2,*) frame_id, boxsz(1,iframe), boxsz(2,iframe), boxsz(3,iframe)
! 	Do imol = 1, n_mol_tot
! 		!Write (3,'(A,F15.7,F15.7,F15.7)') 'O', rx(imol,iframe), ry(imol,iframe), rz(imol,iframe)
		
! 		Do isite = 1, n_mol_sites
! 			! Write to file
! 			Write(2,*) atom_name(isite), rx_s(isite,imol,iframe), ry_s(isite,imol,iframe), rz_s(isite,imol,iframe)
! 		End Do
! 	End do
! End Do



Call CPU_TIME(cpu_st)
! --------------- Start postprocessing ----------------------!
Do iframe = 1, ns_frame
		
	! ----------- Sampling molecular r-density ------------!
	! Loop over all molecules
	Do imol = 1, n_mol_tot

		! Calculate r-distance of imol in spherical coordiantes
		radistsq(imol) = rx(imol,iframe)**2 + ry(imol,iframe)**2 + rz(imol,iframe)**2
		radist(imol) = dSQRT(radistsq(imol))

		If (radist(imol) .LT. rden_cut) Then
			! Calculate ibin number
			ibin = FLOOR(radist(imol)/rden_dr) + 1

			! Accumulate number 
			rdenblk(ibin,iframe) = rdenblk(ibin,iframe) + 1.0d0
		End If

	! End loop over all molecules 
	End Do

	! Loop over bins
	Do ibin = 1, rden_bins
		! Convert to number density (1/A^3)
		rdenavg(ibin) = rdenavg(ibin) + rdenblk(ibin,iframe)/rden_dvol(ibin)		 
	End Do  
		

	! ----- Sampling molecular spherical pressure tensor -------!   
	! Loop over all molecules
	DO imol = 1, n_mol_tot - 1

		rxi = rx(imol,iframe) 
		ryi = ry(imol,iframe) 
		rzi = rz(imol,iframe)

		Do jmol = imol + 1, n_mol_tot
			
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


			! loop over sites on imol
			Do isite = 1, n_mol_sites

				rxis = rx_s(isite,imol,iframe) 
				ryis = ry_s(isite,imol,iframe) 
				rzis = rz_s(isite,imol,iframe)

				! loop over sites on jmol
				Do jsite = 1, n_mol_sites 

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

					! -- Assume imol in the central box 
					! Calculate coefficient in ax^2+bx+c form for R=R_l use
					aa = rijsq
					bb = 2.0*(rxi*xij + ryi*yij + rzi*zij)

					!Loop through the bins
					DO ibin= 1, rden_bins

						! Calculate alpha_k
						cc = radistsq(imol) - posr(ibin)**2
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
								xl = rxi + alpk(isol)*xij
								yl = ryi + alpk(isol)*yij
								zl = rzi + alpk(isol)*zij    
								rijer = xij * xl/posr(ibin) + yij * yl/posr(ibin) + zij * zl/posr(ibin)
								rijser = xijs * xl/posr(ibin) + yijs * yl/posr(ibin) + zijs * zl/posr(ibin)

								! Update normal pressure in radial direction
								pnrc = rijer * rijser * f_tot/(dABS(rijer)*rijs) 
								virialpress_pnr(ibin) = virialpress_pnr(ibin) + pnrc

								! - Tangential pressure in polar direction (theta)
								! Calculate basic parameters
								dxy = dSQRT(xl**2 + yl**2)
								rijet = xij * (xl*zl)/(posr(ibin)*dxy) + yij * (yl*zl)/(posr(ibin)*dxy) - zij * dxy/posr(ibin)
								rijset = xijs * (xl*zl)/(posr(ibin)*dxy) + yijs * (yl*zl)/(posr(ibin)*dxy) - zijs * dxy/posr(ibin)

								! update tangential pressure in theta-theta direction
								pttc = rijet * rijset/dABS(rijer) * f_tot/rijs 
								virialpress_ptt(ibin) = virialpress_ptt(ibin) + pttc

								! - Tangential pressure in azimuthal direction (phi)
								! Calculate basic parameters
								rijep = yij * xl/dxy - xij * yl/dxy
								rijsep = yijs * xl/dxy - xijs * yl/dxy

								! update tangential pressure in phi-phi direction
								ptpc = rijep * rijsep/dABS(rijer) * f_tot/rijs 
								virialpress_ptp(ibin) = virialpress_ptp(ibin) + ptpc

							! end check alpk
							End if
						! End loop over solutions
						End Do
					!End bin cycle
					ENDDO

					! Special treatment of Periodic boundary condition
					If((xij .NE. (rxj -rxi)) .or. (yij .NE. (ryj -ryi)) .or. (zij .NE. (rzj-rzi))) Then

						! -- Assume j site in the central box and ij pair contribute again
						rxii = rxj - xij
						ryii = ryj - yij
						rzii = rzj - zij

						! Calculate coefficient in ax^2+bx+c form for R=R_l use
						aa = rijsq
						bb = 2.0*(rxii*xij + ryii*yij + rzii*zij)

						!Loop through the bins
						DO ibin= 1, rden_bins

							! Calculate alpha_k
							cc = rxii**2 + ryii**2 + rzii**2 - posr(ibin)**2
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

								! The following code is slightly different from above with rxi/ryi/rzi replaced by rxii/ryii/rzii
								if((alpk(isol) .ge. 0.0d0) .and. (alpk(isol) .le. 1.0d0)) Then

									! - Normal pressure component
									! Calculate basic parameters
									xl = rxii + alpk(isol)*xij
									yl = ryii + alpk(isol)*yij
									zl = rzii + alpk(isol)*zij    
									rijer = xij * xl/posr(ibin) + yij * yl/posr(ibin) + zij * zl/posr(ibin)
									rijser = xijs * xl/posr(ibin) + yijs * yl/posr(ibin) + zijs * zl/posr(ibin)

									! Update normal pressure in radial direction
									pnrc = rijer * rijser * f_tot/(dABS(rijer)*rijs) 
									virialpress_pnr(ibin) = virialpress_pnr(ibin) + pnrc

									! - Tangential pressure in polar direction (theta)
									! Calculate basic parameters
									dxy = dSQRT(xl**2 + yl**2)
									rijet = xij * (xl*zl)/(posr(ibin)*dxy) + yij * (yl*zl)/(posr(ibin)*dxy) - zij * dxy/posr(ibin)
									rijset = xijs * (xl*zl)/(posr(ibin)*dxy) + yijs * (yl*zl)/(posr(ibin)*dxy) - zijs * dxy/posr(ibin)

									! update tangential pressure in theta-theta direction
									pttc = rijet * rijset/dABS(rijer) * f_tot/rijs 
									virialpress_ptt(ibin) = virialpress_ptt(ibin) + pttc

									! - Tangential pressure in azimuthal direction (phi)
									! Calculate basic parameters
									rijep = yij * xl/dxy - xij * yl/dxy
									rijsep = yijs * xl/dxy - xijs * yl/dxy

									! update tangential pressure in phi-phi direction
									ptpc = rijep * rijsep/dABS(rijer) * f_tot/rijs 
									virialpress_ptp(ibin) = virialpress_ptp(ibin) + ptpc

								! end check alpk
								End if

							! end loop over solutions		
							End Do

						!End bin cycle
						ENDDO

					! End PBC treatment
					End If


				! End loop over jsite
				End Do

			! End loop over isite
			End Do

		!End loop over jmol sites 
		ENDDO     

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
		
	rdenavg(ibin) = rdenavg(ibin)/DBLE(ns_frame)
	! Convert unit from [1/A^3] to [g/ml]
	!rdenavg(ibin,itype) = (mol_mass/(Na*1.0d-24))*rdenavg(ibin,itype)

End Do

! Write data to file
! Write molecule info
!Write(2,'(2A)') 'Molecule name: ', site_name
Write(2,*) ' R             R-rho [1/A^3]'

! Loop over bins
Do ibin = 1, rden_bins

		! Write to file in units of [1/A^3]
		Write(2,'(F8.4,8X,E15.7)')  posr(ibin), rdenavg(ibin)
																				
End Do

! Close file
CLOSE(2)

! ------------ Spherical molecular pressure tensor ----------------!
Do ibin = 1, rden_bins

	! Pressure (Kinetic part) in unit of [K/A^3]
	prkin(ibin) = rdenavg(ibin)*temp

	! Irving-Kirkwood route
	virialpnravg(ibin) = virialpress_pnr(ibin)/DBLE(ns_frame)
	virialpttavg(ibin) = virialpress_ptt(ibin)/DBLE(ns_frame)
	virialptpavg(ibin) = virialpress_ptp(ibin)/DBLE(ns_frame)

	! Calculate final pressure (configurational part) in unit of [K/A^3]
	virialpnravg(ibin) = 1.0d0/(2.0*two_Pi*posr(ibin)**2)*virialpnravg(ibin)
	virialpttavg(ibin) = 1.0d0/(2.0*two_Pi*posr(ibin)**2)*virialpttavg(ibin)
	virialptpavg(ibin) = 1.0d0/(2.0*two_Pi*posr(ibin)**2)*virialptpavg(ibin)

End Do

! Open file
OPEN(3,FILE='press_spheIK.txt', STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

! Write file head
Write(3,*) 'Spherical molecular pressure tensor from virial route using Irving-Kirkwood definition'
Write(3,*) 'Unit: pressure in [bar] and length in [Angstrom]'
Write(3,'(A)') &
								& '  R               Pkin          Pnr(c)             Pnr(tot)          Ptt(c)&
								&         Ptt(tot)         Ptp(c)          Ptp(tot)'

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
Deallocate(rx)
Deallocate(ry)
Deallocate(rz)
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






