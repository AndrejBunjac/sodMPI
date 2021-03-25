!     This is the first parallel working version. Only the subroutine part is parallelized.
!     Next version should parallelize the time evolution part as well.
!     Debugging commands are removed.

      program MPISOD
      include 'mpif.h'

!     =================================================================
!     Second Order Difference Wave Packet Propagation 
!     Method in Cylindrical Coordinates (rho,z)
!     =================================================================

!     =================================================================
!     Declaration of SOD variables and constants
!     =================================================================
!     Variables
      real*8 Lq,Lz,norm,Vecp,t,omf,ef
      real*8 q,z
      complex acf
!     For automating file generation
      integer numfiles,tdiv,file1,file2
!     Sequences
      integer, dimension(2) :: nn
      real*8, dimension(1024,2048) :: psir,vabs,vpot
      COMPLEX*16, dimension(1024,2048) :: psiN3L0, psi, psi0, psi1
      complex, dimension(1024,2048) :: wf
!     Constants
      COMPLEX*16, parameter :: IU = (0,1)
      real*8, parameter :: pi = 3.1415926535898d0
      real*8 qmin,zmin,dq,dz,qc
!     Hellman potential parameters A*e^(-al r)/r
      real*8, parameter :: A = 21.d0 
      real*8, parameter :: al = 2.54920d0
!     Miscalaneous
      character*16 fname1,fname2
      character*4 tnamesuf

!     =================================================================
!     Declaration of MPI variables and parameters
!     =================================================================
!     message tags and descriptors
      integer,parameter :: BEGIN=1
      integer,parameter :: LNGHBOR=2
      integer,parameter :: RNGHBOR=3
      integer,parameter :: DONE=4
      integer,parameter :: NONE=0
      integer,parameter :: MASTER=0
      integer,parameter :: MINWORKER=3
!     set maxworker to be less than nn(2)/5
      integer,parameter :: MAXWORKER=200
      integer taskid,numtasks,nworkers,avecol,cols,offset,extra,check
      integer dest,source,nleft,nright,msgtype,rc,start,endd,ierr
      integer,dimension(MPI_STATUS_SIZE) :: status1


!     =================================================================
!     Initialize MPI, get ranks, check nuumber of tasks
!     =================================================================

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, ierr )
      nworkers = numtasks - 1

!     =================================================================
!     ************************* MASTER CODE ***************************
!     =================================================================
      if (taskid.eq.MASTER) then


!     check workers
      if ((nworkers.lt.MINWORKER).or.(nworkers.gt.MAXWORKER)) then
        print *,'MP_PROCS needs to be between',MINWORKER+1,'and',MAXWORKER+1,'for this code'
        call MPI_FINALIZE( ierr )
        stop
      endif

!     =================================================================
!     ------------------SERIAL PART OF MASTER CODE---------------------
!     =================================================================


!     =================================================================
!     User input
!     =================================================================
!     Data log file
      open(2,file='out.dat',status='unknown')
!
!     Laser parameters
!
      write(*,*) 'Laser parameters:'
      write(*,*) 'Laser field strength F<recomended - 0.0158351>:'
      read(*,*) ef
      write(*,*) 'Laser field frequency om<recomended - 0.0569538>:'
      read(*,*) omf
      write(2,*) 'F =',real(ef)
      write(2,*) 'omf =',real(omf)
!
!     Grid parameters
! 
      write(*,*) 'Grid parameters:'
      write(*,*) 'Grid size Lz<1000>, Nz<2048>:'
      read(*,*) Lz,nn(2)
      write(*,*) 'Grid start z_min<-500>:'
      read(*,*) zmin
!     NOTE: All z values of grid points will be shifted by dz/2 along the axis (to avoid the singularity at z=0)
      Lq = Lz/2
      nn(1) = nn(2)/2
!
!     Time parameters
!
      write(*,*) 'Time parameters:'
      write(*,*) 'tp<2356.437>, tmax<3000>, dt<.002>:'
      read(*,*) tp,tmax,dt
      write(*,*) 'tdiv - regulates the otput file count, set to 100/dt <50000>'
      read(*,*) tdiv
      write(*,*) 'rho_c<5>, r0<450>, V0abs<.03-.1>:'
      read(*,*) qc,r0,v0a
      write(*,*) 'nstep<1>:'
      read(*,*) nst
      dq = Lq/nn(1)
      dz = Lz/nn(2)
      qmin = 0.d0
      numfiles = Floor(tmax/100)

!     Update values in the log file
      write(2,*) 'Lz =',real(Lz),' zmin =',real(zmin)
      write(2,*) 'Lrho =',real(Lq),' rho_min =',real(qmin)
      write(2,*) 'Nz =',nn(2)
      write(2,*) 'Nrho =',nn(1)
      write(2,*) 'Tmax =',real(tmax),' dt =',real(dt)
      write(2,*) 'rho_c =',real(qc),' r0 =',real(r0),' V0a =',real(v0a)

!     =================================================================
!     Send initial data to workers
!     =================================================================

      call MPI_Barrier( MPI_COMM_WORLD )
      call MPI_BCAST ( qmin, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )
      call MPI_BCAST ( zmin, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( dq, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( dz, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( qc, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )
      call MPI_BCAST ( nn(1), 2, MPI_INTEGER, MASTER, MPI_COMM_WORLD, ierr )        

!     =================================================================
!     Preparing output files
!     =================================================================
      print *,'Preparing output files...'
!      
      OPEN(11,FILE='acf.dat',STATUS='UNKNOWN')
      OPEN(12,FILE='reacf.dat',STATUS='UNKNOWN')
      OPEN(13,FILE='cf2N3L0.dat',STATUS='UNKNOWN')
!
      do i=1,numfiles
!
        write(tnamesuf, 50) i*100
50      format (I4)
        fname1 = trim('wf'//trim(adjustl(tnamesuf))//'.dat')
        file1 = 500+i
        fname2 = trim('wf2-'//trim(adjustl(tnamesuf))//'.dat')
        file2 = 500+numfiles+i
!
        open(file1,FILE=fname1,status='UNKNOWN')
        open(file2,FILE=fname2,status='UNKNOWN')
!
      enddo

!     =================================================================
!     Opening input file and preparing the initial wave-packet
!     =================================================================
!
      print *,'Preparing the initial wave packet...'
      OPEN(21,FILE='N3L0NaSOD2048.dat',status='old')

!     Initial wave packet
	
100	Read (21,*,END=150) i,j,psir(i,j)
	Go TO 100
150	close(21)
!
!     *** psir - loads raw data from a file
!     *** psiN3L0 = calculates the regularized function by multiplying with sqrt(q/2)
!     *** keeps psiN3L0 in memory for later calculation
!     *** copies psiN3L0 to psi to be used as a variable
!
!     Initial wave packet and absorbing pot.
!  
      DO j=1,nn(2)
        z = zmin + dz/2 + (j-1)*dz
        DO i=1,nn(1)
          q = qmin + (i-1)*dq
          r = dsqrt(q*q+z*z)	
!
          psiN3L0(i,j) = dsqrt(q/2)*psir(i,j)
!
          psi(i,j) = psiN3L0(i,j)

          if (r.gt.r0) vabs(i,j) = v0a*(r-r0)**2
        enddo
      enddo
!     ================================================================
!     Initialize MPI data to workers
!     ================================================================
!     avecol - average number of columns taken per process
!     extra the remaining columns redistributed to the processes
      avecol = nn(2)/nworkers
      extra = mod(nn(2),nworkers)
      offset = 1
!
      do i=1,nworkers
        if (i.le.extra) then
          cols = avecol + 1
        else
          cols = avecol
        endif
!     figure out neighbours for each worker for the purpose of exchaning border data
      if (i.eq.1) then
        nleft = NONE
      else
        nleft = i - 1
      endif
      if (i.eq.nworkers) then
        nright = NONE
      else
        nright = i + 1
      endif
!     send to workers
      dest = i
      call MPI_SEND( offset, 1, MPI_INTEGER, dest, BEGIN, MPI_COMM_WORLD, ierr )
      call MPI_SEND( cols, 1, MPI_INTEGER, dest, BEGIN, MPI_COMM_WORLD, ierr )
      call MPI_SEND( nleft, 1, MPI_INTEGER, dest, BEGIN, MPI_COMM_WORLD, ierr )
      call MPI_SEND( nright, 1, MPI_INTEGER, dest, BEGIN, MPI_COMM_WORLD, ierr )
!     update offset
      offset = offset + cols
      enddo

!     ================================================================
!     THE TIME PROPAGATION - MAIN PART OF THE CODE
!     ================================================================
!
!     t - time in atomic units
!     it - number of time counts (iterations)
!     icount - counts the file adress to write to
!
      print *,'Initializing time propagation...'
      t = 0.d0
      it = 0
      icount = 501
      acf = (1.d0,0.d0)
!     0-th iteration data in the files
      write(11,*) real(t),acf
      write(12,*) real(t),real(acf)
!
!     RETURN POINT FOR TIME PROPAGATION!!!
!
   10 CONTINUE
      it = it+1
      t = it*dt
      call MPI_Barrier( MPI_COMM_WORLD )
      call MPI_BCAST ( t, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )

!     ===============================================================
!     Potential - calculation of the potential on the grid
!     ===============================================================

!     Set border values to be 0
      do i=1,nn(1)
        vpot(i,1) = 0.d0
        vpot(i,nn(2)) = 0.d0
      enddo
      do j=1,nn(2)
        vpot(1,j) = 0.d0
        vpot(nn(1),j) = 0.d0
      enddo

!
!     Calculate the potential (except borders) using OMP parallelization
!
      !$OMP PARALLEL DO PRIVATE(q,z,r,vecp,vfld) SHARED(dq,dz,qmin,zmin,nst,nn,vpot)
      do i=2,nn(1)-1,nst
        q = qmin + (i-1)*dq
        do j=2,nn(2)-1,nst
          z = zmin + dz/2 + (j-1)*dz
          r = dsqrt(q*q+z*z)
          Vecp = -1.d0/r + A*dexp(-al*r)/r
!      Seting the potential to be 0 around the singularity for q
          if (q.le.5d-7) then
              vpot(i,j) = 0.
            else
!      The pulse turns off at t = tp       
            if (t.lt.tp) then
              vpot(i,j) = Vecp - 1.d0/(8*q*q) - z*ef*(cos(pi*(t-tp/2)/tp))**2*cos(omf*(t-tp/2)+phi)
              else
              vpot(i,j) = Vecp - 1.d0/(8*q*q)
            endif
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO

!     ===============================================================
!     ---------------PARALLEL PART OF MASTER CODE--------------------
!     ===============================================================
!                Send psi to workers to update to psi1
!     ===============================================================

!     work distribution to workers
!     send the data chunk of psi to each worker
      offset = 1
      do i=1,nworkers
        if (i.le.extra) then
          cols = avecol + 1
        else
          cols = avecol
        endif
        dest = i
        call MPI_SEND( psi(1,offset), cols*nn(1), MPI_DOUBLE_COMPLEX, dest, BEGIN, MPI_COMM_WORLD, ierr )
        offset = offset + cols
      enddo

!     wait for worker data results
      do i=1,nworkers
        source = i
        msgtype = DONE
        call MPI_RECV( offset, 1, MPI_INTEGER, source, msgtype, MPI_COMM_WORLD, status1, ierr )
        call MPI_RECV( cols, 1, MPI_INTEGER, source, msgtype, MPI_COMM_WORLD, status1, ierr )
        call MPI_RECV( psi1(1,offset), cols*nn(1), MPI_DOUBLE_COMPLEX, source, msgtype, MPI_COMM_WORLD, status1, ierr )
      enddo

!     ===============================================================
!     -------------- SERIAL PART OF THE MASTER CODE -----------------
!     ===============================================================
!        SOD - The wave function is updated to the next iteration
!     ===============================================================

!
!     OMP parallelization to calculate the new psi
!
      !$OMP PARALLEL DO SHARED(nn,psi,psi0,psi1,vpot,vabs)
      do j=1,nn(2)
        do i=1,nn(1)
          if (it.eq.1) then
            psi1(i,j) = psi(i,j) - iu*dt*psi1(i,j) - iu*dt*vpot(i,j)*psi(i,j)      
          else
            psi1(i,j) = psi0(i,j) - 2*iu*dt*psi1(i,j) - 2*iu*dt*vpot(i,j)*psi(i,j)
          endif
          psi0(i,j) = psi(i,j)
          psi(i,j) = psi1(i,j)
!         Absorbing pot.
          psi(i,j) = (1.d0-vabs(i,j)*dt)*psi(i,j)
        enddo
      enddo
      !$OMP END PARALLEL DO

!     ==============================================================
!     Auto-correlation function
!     ==============================================================

      acf = 0.d0
      do j=1,nn(2)
        z = zmin + dz/2 + (j-1)*dz
        do i=1,nn(1)
          q = qmin + (i-1)*dq
          r = dsqrt(q*q+z*z)
          acf = acf + psiN3L0(i,j)*psi(i,j)
        enddo
      enddo
      acf = acf*dq*dz*4*pi
!
!     Division by 10 is to make sure that every 10th value is written in the file
!
      if (dfloat(it)/10-it/10.eq.0.d0) then
        write(*,200) t,real(acf)
  200	  format(' t =',g10.4,' Re(acf) =',g10.4)
        write(11,*) real(t),acf
        write(12,*) real(t),real(acf)
        write(13,*) real(t),abs(acf)**2
      endif

!     ==============================================================
!     Wave Function Frames
!     ==============================================================

      if (icount.le.500+numfiles) then
        if (float(it)/tdiv-it/tdiv.eq.0.) then
          do j=1,nn(2)
            do i=1,nn(1)
              if (i.gt.1) then
                q = qmin + (i-1)*dq
                wf(i,j) = psi(i,j)/dsqrt(2*pi*q)
              endif
            enddo
          enddo
!         estimation:
          do j=1,nn(2)
            wf(1,j) = (4*wf(2,j)-wf(3,j))/3
          enddo
!
          do i=1,nn(1)
              do j=1,nn(2)
                write(icount+numfiles,*) cabs(wf(i,j))**2
                write(icount,*) i,j,wf(i,j)
 	      enddo
          enddo
          icount = icount + 1
        endif
      endif
!
!     reiterate time 
      if (t.lt.tmax) then
        check = 0
      else
        check = 1
      endif
      call MPI_Barrier( MPI_COMM_WORLD )
      call MPI_BCAST ( check, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, ierr )

  
      if (check.eq.0) go to 10


!     ==============================================================
!     Norm
!     ==============================================================
!
      spsi2 = 0.d0
      do j=1,nn(2)
        do i=1,nn(1)
          spsi2 = spsi2 + cdabs(psi(i,j))**2
        enddo
      enddo
      psinorm = dsqrt(spsi2*dq*dz*4*pi)
      write(*,*) 't=',real(t),' |Psi|=',real(psinorm)
      write(*,*) 'F=',ef,' Om=',omf
      write(2,*) 't=',real(t),' |Psi|=',real(psinorm)
!
!     end of MASTER CODE
      call MPI_FINALIZE ( ierr )
      endif


!     ==============================================================
!     ************************ WORKER CODE *************************
!     ==============================================================

!     ================================================================
!     Recieve initial data from MASTER
!     ================================================================

      if (taskid.ne.MASTER) then

      call MPI_Barrier( MPI_COMM_WORLD )
      call MPI_BCAST ( qmin, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )
      call MPI_BCAST ( zmin, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( dq, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( dz, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( qc, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )  
      call MPI_BCAST ( nn(1), 2, MPI_INTEGER, MASTER, MPI_COMM_WORLD, ierr )

!     ================================================================
!     Recieve data from MASTER
!     ================================================================
      source = MASTER
      msgtype = BEGIN
      call MPI_RECV( offset, 1, MPI_INTEGER, source, msgtype, MPI_COMM_WORLD, status1, ierr )
      call MPI_RECV( cols, 1, MPI_INTEGER, source, msgtype, MPI_COMM_WORLD, status1, ierr )
      call MPI_RECV( nleft, 1, MPI_INTEGER, source, msgtype, MPI_COMM_WORLD, status1, ierr )
      call MPI_RECV( nright, 1, MPI_INTEGER, source, msgtype, MPI_COMM_WORLD, status1, ierr )

!     Determine border elements. First and last columns for updating 
      start = offset
      endd = offset + cols - 1
    
!     Return point for time iteration of the worker code
   20 CONTINUE
!     Initialize grid to 0
!
      call MPI_Barrier( MPI_COMM_WORLD )
      call MPI_BCAST ( t, 1, MPI_REAL8, MASTER, MPI_COMM_WORLD, ierr )
!  
      DO j=1,nn(2)
        DO i=1,nn(1)	
          psi(i,j) = 0.0
          psi1(i,j) = 0.0
        enddo
      enddo

!     Recieve grid partition from the master
      source = MASTER
      msgtype = BEGIN
      call MPI_RECV( psi(1,offset), cols*nn(1), MPI_DOUBLE_COMPLEX, source, msgtype, MPI_COMM_WORLD, status1, ierr )

!     ==================================================================================================
!     Communicate border elements with neighbors. Need 4 border elements per neighbor for the subroutine
!     ==================================================================================================

!     Send to left
      if (nleft .ne. NONE) then
        call MPI_SEND( psi(1,offset), 4*nn(1), MPI_DOUBLE_COMPLEX, nleft, RNGHBOR, MPI_COMM_WORLD, ierr )
      endif

!     Receive from right
      if (nright .ne. NONE) then
        source = nright
        msgtype = RNGHBOR
        call MPI_RECV( psi(1,offset+cols), 4*nn(1), MPI_DOUBLE_COMPLEX, source, msgtype, MPI_COMM_WORLD, status1, ierr )
      endif

!     Send to right
      if (nright .ne. NONE) then
        call MPI_SEND( psi(1,offset+cols-4), 4*nn(1), MPI_DOUBLE_COMPLEX, nright, LNGHBOR, MPI_COMM_WORLD, ierr)
      endif

!     Recieve from left
      if (nleft .ne. NONE) then
        source = nleft
        msgtype = LNGHBOR
        call MPI_RECV( psi(1,offset-4), 4*nn(1), MPI_DOUBLE_COMPLEX, source, msgtype, MPI_COMM_WORLD, status1, ierr )
      endif

!     Call the subroutine to update psi to psi1
      call HamFun(start,endd,nn,psi,psi1,qmin,zmin,dq,dz,qc) 

!     Send my portion of the results back to master process
      call MPI_SEND( offset, 1, MPI_INTEGER, MASTER, DONE, MPI_COMM_WORLD, ierr )
      call MPI_SEND( cols, 1, MPI_INTEGER, MASTER, DONE, MPI_COMM_WORLD,  ierr )
      call MPI_SEND( psi1(1,offset), cols*nn(1), MPI_DOUBLE_COMPLEX, MASTER, DONE, MPI_COMM_WORLD, ierr )

!     Recieve time count from the MASTER process
      call MPI_Barrier( MPI_COMM_WORLD )
      call MPI_BCAST ( check, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, ierr )

!     end of WORKER CODE


      if (check.eq.0) go to 20
      call MPI_FINALIZE ( ierr )
      endif


      stop
      end


!     =================================================================
!     *****************************************************************
!     Subroutine that updates the wave function
!     *****************************************************************
!     =================================================================

      SUBROUTINE HamFun(start,endd,nn,f,hf,qmin,zmin,dq,dz,qc)
!
!     Calculates the operation T*f(q,z),
!     T=.5*d^2/dt^2. 
!
!     =================================================================
!     Declaration of SOD variables and constants
!     =================================================================

      integer,dimension(2) :: nn
      complex*16, parameter :: IU = (0,1)
      complex*16, dimension(1024,2048) :: f,Tf,Hf,u
      complex*16  uq,uqq,fqq,fzz
      real*8 qmin,zmin,dq,dz,qc
      real*8 q,z,sqrq
      integer start,endd,rstart,rendd
        

!     =================================================================
!     Main part of the subroutine
!     =================================================================

!
!     The kinetic energy part
!
!     regularized function u(q,z) = f(q,z)/sqrt(q)
!
!     rstart and rendd make sure that the extra elements (beyond ones being computed) are also regularized
!
      if (start.eq.1) then
        rstart = 1
      else
        rstart = start - 4
      endif
!
      if (endd.eq.nn(2)) then
        rendd = nn(2)
      else
        rendd = endd + 4
      endif
!
      do j=rstart,rendd
        do i=2,nn(1)
          q = qmin + (i-1)*dq
          u(i,j) = f(i,j)/dsqrt(q)
        enddo
      enddo
!     estimation:
      do j=rstart,rendd
        u(1,j) = (4*u(2,j)-u(3,j))/3
      enddo
!
!     derivatives
!     rstart and rendd make sure that the processes handling the first and last column skip the first and last elements
!
!     
      if (start.eq.1) then
        rstart = 2
      else
        rstart = start
      endif
      if (endd.eq.nn(2)) then
        rendd = nn(2) - 1
      else
        rendd = endd
      endif
!
!
      !$OMP PARALLEL DO PRIVATE(z,q,uq,uqq,sqrq,fqq,fzz) SHARED(u,f,Tf,zmin,qmin,dq,dz,rstart,rendd)
      do j=rstart,rendd
        do i=2,nn(1)-1
          q = qmin + (i-1)*dq
          if (q.le.5d-7) then
            fqq = 0.d0
            fzz = 0.d0
          else
            if (q.le.qc) then
!     ---     calc. fqq in terms of derivatives of u
              if ((i.eq.2).or.(i.eq.nn(1)-1)) then
                uq = (u(i+1,j)-u(i-1,j))/(2*dq)
                uqq = (u(i+1,j)-2*u(i,j)+u(i-1,j))/(dq*dq)
              endif
              if ((i.eq.3).or.(i.eq.nn(1)-2)) then
                uq = (u(i-2,j)-8*u(i-1,j)+8*u(i+1,j)-u(i+2,j))/(12*dq)
                uqq = (-u(i-2,j)+16*u(i-1,j)-30*u(i,j)+16*u(i+1,j)-u(i+2,j))/(12*dq*dq)
              endif
              if ((i.eq.4).and.(i.le.nn(1)-3)) then
                uq = (-u(i-3,j)+9*u(i-2,j)-45*u(i-1,j)+45*u(i+1,j)-9*u(i+2,j)+u(i+3,j))/(60*dq)
                uqq = (2*u(i-3,j)-27*u(i-2,j)+270*u(i-1,j)-490*u(i,j)+270*u(i+1,j)-27*u(i+2,j)+2*u(i+3,j))/(180*dq*dq)
              endif
              if ((i.ge.5).and.(i.le.nn(1)-4)) then
                uq = (u(i-4,j)/280-4*u(i-3,j)/105+u(i-2,j)/5-4*u(i-1,j)/5+4*u(i+1,j)/5-u(i+2,j)/5+4*u(i+3,j)/105-u(i+4,j)/280)/dq
                uqq = (-u(i-4,j)/560+8*u(i-3,j)/315-u(i-2,j)/5+8*u(i-1,j)/5-205*u(i,j)/72+8*u(i+1,j)/5-  &
                      u(i+2,j)/5+8*u(i+3,j)/315-u(i+4,j)/560)/(dq*dq)
              endif
              sqrq = dsqrt(q)
              fqq = -u(i,j)/(4*sqrq**3) + uq/sqrq + sqrq*uqq
!     ---
            else
!     ---     calc. fqq directly
              if ((i.eq.2).or.(i.eq.nn(1)-1)) then
                fqq = (f(i+1,j)-2*f(i,j)+f(i-1,j))/(dq*dq)
              endif
              if ((i.eq.3).or.(i.eq.nn(1)-2)) then
                fqq = (-f(i-2,j)+16*f(i-1,j)-30*f(i,j)+16*f(i+1,j)-f(i+2,j))/(12*dq*dq)
              endif
              if ((i.eq.4).or.(i.eq.nn(1)-3)) then
                fqq = (2*f(i-3,j)-27*f(i-2,j)+270*f(i-1,j)-490*f(i,j)+270*f(i+1,j)-27*f(i+2,j)+2*f(i+3,j))/(180*dq*dq)
              endif
              if ((i.ge.5).and.(i.le.nn(1)-4)) then
                fqq = (-f(i-4,j)/560+8*f(i-3,j)/315-f(i-2,j)/5+8*f(i-1,j)/5-205*f(i,j)/72+8*f(i+1,j)/5-  &
                      f(i+2,j)/5+8*f(i+3,j)/315-f(i+4,j)/560)/(dq*dq)
              endif
!     ---
            endif
!     ---   calc. fzz
            if ((j.eq.2).or.(j.eq.nn(2)-1)) then
              fzz = (f(i,j+1)-2*f(i,j)+f(i,j-1))/(dz*dz)
            endif
            if ((j.eq.3).or.(j.eq.nn(2)-2)) then
              fzz = (-f(i,j-2)+16*f(i,j-1)-30*f(i,j)+16*f(i,j+1)-f(i,j+2))/(12*dz*dz)
            endif
            if ((j.eq.4).or.(j.eq.nn(2)-3)) then
              fzz = (2*f(i,j-3)-27*f(i,j-2)+270*f(i,j-1)-490*f(i,j)+270*f(i,j+1)-27*f(i,j+2)+2*f(i,j+3))/(180*dz*dz)
            endif
            if ((j.ge.5).and.(j.le.nn(2)-4)) then
              fzz = (-f(i,j-4)/560+8*f(i,j-3)/315-f(i,j-2)/5+8*f(i,j-1)/5-205*f(i,j)/72+8*f(i,j+1)/5-  &
                    f(i,j+2)/5+8*f(i,j+3)/315-f(i,j+4)/560)/(dz*dz)
            endif
!     ---
          endif
          Tf(i,j) = -(fqq + fzz)/2
        enddo
      enddo
      !$OMP END PARALLEL DO
!
!     The potential energy part and H*f=(T+V)*f
!
      if (start.eq.1) then
        do i=1,nn(1)
          Hf(i,1) = 0.d0
        enddo
      endif
      if (endd.eq.nn(2)) then
        do i=1,nn(1)
          Hf(i,nn(2)) = 0.d0
        enddo
      endif
!
      do j=start,endd
        Hf(1,j) = 0.d0
        Hf(nn(1),j) = 0.d0
      enddo
!
      !$OMP PARALLEL DO PRIVATE(z,q,uq,uqq,sqrq,fqq,fzz) SHARED(nn,Hf,Tf,qmin,zmin,dq,dz,rstart,rendd)
      do j=rstart,rendd
        do i=2,nn(1)-1
          q = qmin + (i-1)*dq
          if (q.le.5d-7) then
	      Hf(i,j) = 0.
          else
	      Hf(i,j) = Tf(i,j)
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO
      RETURN
      END
