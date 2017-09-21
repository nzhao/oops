*----------------------------------------------------------------------|
!      subroutine ZGEXPV( n, m, t, v, w, tol, anorm,
!     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

      SUBROUTINE MPI_HZCSRKEXPV(n, m, v, w, tol, anorm, 
     .  wsp, lwsp, iwsp, liwsp, itrace, iflag, 
     .  nspin, h, prefactor, nterm, tlst, tn, alst, n_op, oplst, 
     .  reqs, status, cache, nchunk)
      implicit none
      
      INCLUDE "mpif.h"
      
      INTEGER           nspin, nterm, tn, nchunk, ierr, n_op
      DOUBLE PRECISION  tlst(tn)
      COMPLEX*16        oplst(tn), cache(2 * n)
      INTEGER           rank, nprocess, idx, idx1
      INTEGER           reqs(2 * nchunk)
      INTEGER           status(MPI_Status_size * 2 * nchunk)
      INTEGER           tt
      COMPLEX*16        cz
      DOUBLE PRECISION  cr
      COMPLEX*16        ONEIM
      PARAMETER( ONEIM=(0.0d0,-1.0d0) )
      INTEGER*8         h
      INTEGER*8         alst(n_op * 20)
      COMPLEX*16        prefactor
      
      integer          n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm
      complex*16       v(n), w(n), wsp(lwsp)
!      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZGEXPV computes w = exp(t*A)*v
*     for a Zomplex (i.e., complex double precision) matrix A 
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        complex*16 x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|

      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hij, op

      intrinsic AINT,ABS,DCMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
      
*      write(*,'(a,I16)'),'[mpi_zcsrkexpv] n     = ',n    
*      write(*,'(a,I16)'),'[mpi_zcsrkexpv] nspin = ',nspin    
*      write(*,'(a,I16)'),'[mpi_zcsrkexpv] lwsp  = ',lwsp    
*      write(*,'(a,I16)'),'[mpi_zcsrkexpv] nterm = ',nterm    
*      write(*,'(a,Z16)'),'[mpi_zcsrkexpv] n     = ',loc(n)    
*      write(*,'(a,Z16)'),'[mpi_zcsrkexpv] nspin = ',loc(nspin)    
*      write(*,'(a,Z16)'),'[mpi_zcsrkexpv] lwsp  = ',loc(lwsp)    
*      write(*,'(a,Z16)'),'[mpi_zcsrkexpv] nterm = ',loc(nterm)
*      stop
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      CALL MPI_ALLREDUCE(beta, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
     .  MPI_COMM_WORLD, ierr)
      beta = cr
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
      tlst_output: DO tt = 1, tn
      t_out = tlst(tt)
 100  if ( t_now.ge.t_out ) goto 500
      
*---  initialize wsp before mpi_hamvec
!$OMP PARALLEL DO PRIVATE(j)
      DO i = 1, n
        DO j = 1, m + 2
          wsp(i + (j - 1) * n) = ZERO
        END DO
      END DO
      
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         !call matvec( wsp(j1v-n), wsp(j1v) )
!         CALL mpi_hamvec(nspin, coeff_lst, wsp(j1v-n), wsp(j1v), reqs, 
!     .     status, cache, nchunk)
!         CALL ZSCAL(n, ONEIM, wsp(j1v), 1)! -i * h;
         CALL mpi_zcsrkhmv(nspin, h, prefactor, wsp(j1v-n), wsp(j1v), 
     .     nchunk, reqs, status, cache)
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            CALL MPI_ALLREDUCE(hij, cz, 2, MPI_DOUBLE_PRECISION, 
     .        MPI_SUM, MPI_COMM_WORLD, ierr)
            hij = cz
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DZNRM2( n, wsp(j1v),1 )
         CALL MPI_ALLREDUCE(hj1j, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
     .     MPI_COMM_WORLD, ierr)
         hj1j = cr
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = DCMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      !call matvec( wsp(j1v-n), wsp(j1v) )
!      CALL mpi_hamvec(nspin, coeff_lst, wsp(j1v-n), wsp(j1v), reqs, 
!     .  status, cache, nchunk)
!      CALL ZSCAL(n, ONEIM, wsp(j1v), 1)! -i * h;
      CALL mpi_zcsrkhmv(nspin, h, prefactor, wsp(j1v-n), wsp(j1v), 
     .  nchunk, reqs, status, cache)
      avnorm = DZNRM2( n, wsp(j1v),1 )
      CALL MPI_ALLREDUCE(avnorm, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
     .  MPI_COMM_WORLD, ierr)
      avnorm = cr
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = ONE
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         enddo
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue
* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZGEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hij = DCMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      CALL MPI_ALLREDUCE(beta, cr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
     .  MPI_COMM_WORLD, ierr)
      beta = cr
      hump = MAX( hump, beta )

!*initialize wsp before mpi_hamvec
!!$OMP PARALLEL DO PRIVATE(j)
!      do i=1,n
!          do j=1,m+2
!              wsp( i+(j-1)*n )= ZERO
!          enddo
!      enddo
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue
      
      do j=1,n_op
!$OMP PARALLEL DO
        do i=1,n
           wsp( i )= ZERO
        enddo
        
        CALL mpi_zcsrkhmv(nspin, alst((j - 1) * 20 + 1), ONE, w,wsp,
     .    nchunk, reqs, status, cache)
        op = ZDOTC( n, w, 1, wsp,1 )
        CALL MPI_ALLREDUCE(op, cz, 2, MPI_DOUBLE_PRECISION, 
     .       MPI_SUM, MPI_COMM_WORLD, ierr)
        op = cz
        oplst(tt + (j-1)*tn) = op
      enddo
      
      END DO tlst_output

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = DCMPLX( step_min )
      wsp(2)  = DCMPLX( step_max )
      wsp(3)  = DCMPLX( 0.0d0 )
      wsp(4)  = DCMPLX( 0.0d0 )
      wsp(5)  = DCMPLX( x_error )
      wsp(6)  = DCMPLX( s_error )
      wsp(7)  = DCMPLX( tbrkdwn )
      wsp(8)  = DCMPLX( sgn*t_now )
      wsp(9)  = DCMPLX( hump/vnorm )
      wsp(10) = DCMPLX( beta/vnorm )
      END
*----------------------------------------------------------------------|
