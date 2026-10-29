c Downloaded from https://netlib.org/cgi-bin/netlibfiles.pl?filename=/quadpack/qawse.f (public domain) on 29.10.2026
c Changelog:
c - Added d1mach from https://people.math.sc.edu/Burkardt/f77_src/machine/machine.f (GNU LGPL)

      subroutine dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,
     *   result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqawse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1
c***keywords  automatic integrator, special-purpose,
c             algebraico-logarithmic end point singularities,
c             clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f*w over (a,b),
c            (where w shows a singular behaviour at the end points,
c            see parameter integr).
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration of functions having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        double precision version
c
c        parameters
c         on entry
c            f      - double precision
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - double precision
c                     lower limit of integration
c
c            b      - double precision
c                     upper limit of integration, b.gt.a
c                     if b.le.a, the routine will end with ier = 6.
c
c            alfa   - double precision
c                     parameter in the weight function, alfa.gt.(-1)
c                     if alfa.le.(-1), the routine will end with
c                     ier = 6.
c
c            beta   - double precision
c                     parameter in the weight function, beta.gt.(-1)
c                     if beta.le.(-1), the routine will end with
c                     ier = 6.
c
c            integr - integer
c                     indicates which weight function is to be used
c                     = 1  (x-a)**alfa*(b-x)**beta
c                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c                     if integr.lt.1 or integr.gt.4, the routine
c                     will end with ier = 6.
c
c            epsabs - double precision
c                     absolute accuracy requested
c            epsrel - double precision
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.2
c                     if limit.lt.2, the routine will end with ier = 6.
c
c         on return
c            result - double precision
c                     approximation to the integral
c
c            abserr - double precision
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for the integral and error
c                             are less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement, it is advised to analyze the
c                             integrand in order to determine the
c                             integration difficulties which prevent the
c                             requested tolerance from being achieved.
c                             in case of a jump discontinuity or a local
c                             singularity of algebraico-logarithmic type
c                             at one or more interior points of the
c                             integration range, one should proceed by
c                             splitting up the interval at these
c                             points and calling the integrator on the
c                             subranges.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             b.le.a or alfa.le.(-1) or beta.le.(-1), or
c                             integr.lt.1 or integr.gt.4, or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             or limit.lt.2.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            blist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the given integration range (a,b)
c
c            rlist  - double precision
c                     vector of dimension at least limit,the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - double precision
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     of which are pointers to the error
c                     estimates over the subintervals, so that
c                     elist(iord(1)), ..., elist(iord(k)) with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise form a decreasing sequence
c
c            last   - integer
c                     number of subintervals actually produced in
c                     the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqc25s,dqmomo,dqpsrt
c***end prologue  dqawse
c
      double precision a,abserr,alfa,alist,area,area1,area12,area2,a1,
     *  a2,b,beta,blist,b1,b2,centre,dabs,dmax1,d1mach,elist,epmach,
     *  epsabs,epsrel,errbnd,errmax,error1,erro12,error2,errsum,f,
     *  resas1,resas2,result,rg,rh,ri,rj,rlist,uflow
      integer ier,integr,iord,iroff1,iroff2,k,last,limit,maxerr,nev,
     *  neval,nrmax
c
      external f
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit),ri(25),rj(25),rh(25),rg(25)
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqawse
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(b.le.a.or.(epsabs.eq.0.0d+00.and.
     *  epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28)).or.alfa.le.(-0.1d+01).
     *  or.beta.le.(-0.1d+01).or.integr.lt.1.or.integr.gt.4.or.
     *  limit.lt.2) go to 999
      ier = 0
c
c           compute the modified chebyshev moments.
c
      call dqmomo(alfa,beta,ri,rj,rg,rh,integr)
c
c           integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
c
      centre = 0.5d+00*(b+a)
      call dqc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
      neval = nev
      call dqc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
      last = 2
      neval = neval+nev
      result = area1+area2
      abserr = error1+error2
c
c           test on accuracy.
c
      errbnd = dmax1(epsabs,epsrel*dabs(result))
c
c           initialization
c           --------------
c
      if(error2.gt.error1) go to 10
      alist(1) = a
      alist(2) = centre
      blist(1) = centre
      blist(2) = b
      rlist(1) = area1
      rlist(2) = area2
      elist(1) = error1
      elist(2) = error2
      go to 20
   10 alist(1) = centre
      alist(2) = a
      blist(1) = b
      blist(2) = centre
      rlist(1) = area2
      rlist(2) = area1
      elist(1) = error2
      elist(2) = error1
   20 iord(1) = 1
      iord(2) = 2
      if(limit.eq.2) ier = 1
      if(abserr.le.errbnd.or.ier.eq.1) go to 999
      errmax = elist(1)
      maxerr = 1
      nrmax = 1
      area = result
      errsum = abserr
      iroff1 = 0
      iroff2 = 0
c
c            main do-loop
c            ------------
c
      do 60 last = 3,limit
c
c           bisect the subinterval with largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
c
        call dqc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,area1,
     *  error1,resas1,integr,nev)
        neval = neval+nev
        call dqc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,area2,
     *  error2,resas2,integr,nev)
        neval = neval+nev
c
c           improve previous approximations integral and error
c           and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(a.eq.a1.or.b.eq.b2) go to 30
        if(resas1.eq.error1.or.resas2.eq.error2) go to 30
c
c           test for roundoff error.
c
        if(dabs(rlist(maxerr)-area12).lt.0.1d-04*dabs(area12)
     *  .and.erro12.ge.0.99d+00*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
   30   rlist(maxerr) = area1
        rlist(last) = area2
c
c           test on accuracy.
c
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 35
c
c           set error flag in the case that the number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c
c           set error flag in the case of roundoff error.
c
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
c
c           set error flag in the case of bad integrand behaviour
c           at interior points of integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (dabs(a2)+0.1d+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   35   if(error2.gt.error1) go to 40
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 50
   40   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with largest error estimate (to be bisected next).
c
   50   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if (ier.ne.0.or.errsum.le.errbnd) go to 70
   60 continue
c
c           compute final result.
c           ---------------------
c
   70 result = 0.0d+00
      do 80 k=1,last
        result = result+rlist(k)
   80 continue
      abserr = errsum
  999 return
      end
      subroutine dqc25s(f,a,b,bl,br,alfa,beta,ri,rj,rg,rh,result,
     *   abserr,resasc,integr,nev)
c***begin prologue  dqc25s
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (bl,br), with error
c            estimate, where the weight function w has a singular
c            behaviour of algebraico-logarithmic type at the points
c            a and/or b. (bl,br) is a part of (a,b).
c***description
c
c        integration rules for integrands having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           f      - double precision
c                    function subprogram defining the integrand
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - double precision
c                    left end point of the original interval
c
c           b      - double precision
c                    right end point of the original interval, b.gt.a
c
c           bl     - double precision
c                    lower limit of integration, bl.ge.a
c
c           br     - double precision
c                    upper limit of integration, br.le.b
c
c           alfa   - double precision
c                    parameter in the weight function
c
c           beta   - double precision
c                    parameter in the weight function
c
c           ri,rj,rg,rh - double precision
c                    modified chebyshev moments for the application
c                    of the generalized clenshaw-curtis
c                    method (computed in subroutine dqmomo)
c
c           result - double precision
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if b1 = a or br = b.
c                    in all other cases the 15-point kronrod
c                    rule is applied, obtained by optimal addition of
c                    abscissae to the 7-point gauss rule.
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           resasc - double precision
c                    approximation to the integral of abs(f*w-i/(b-a))
c
c           integr - integer
c                    which determines the weight function
c                    = 1   w(x) = (x-a)**alfa*(b-x)**beta
c                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
c                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
c                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
c                                 log(b-x)
c
c           nev    - integer
c                    number of integrand evaluations
c***references  (none)
c***routines called  dqcheb,dqk15w
c***end prologue  dqc25s
c
      double precision a,abserr,alfa,b,beta,bl,br,centr,cheb12,cheb24,
     *  dabs,dc,dlog,f,factor,fix,fval,hlgth,resabs,resasc,result,res12,
     *  res24,rg,rh,ri,rj,u,dqwgts,x
      integer i,integr,isym,nev
c
      dimension cheb12(13),cheb24(25),fval(25),rg(25),rh(25),ri(25),
     *  rj(25),x(11)
c
      external f,dqwgts
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ..., 11, to be used for the computation of the
c           chebyshev series expansion of f.
c
      data x(1) / 0.9914448613 7381041114 4557526928 563d0 /
      data x(2) / 0.9659258262 8906828674 9743199728 897d0 /
      data x(3) / 0.9238795325 1128675612 8183189396 788d0 /
      data x(4) / 0.8660254037 8443864676 3723170752 936d0 /
      data x(5) / 0.7933533402 9123516457 9776961501 299d0 /
      data x(6) / 0.7071067811 8654752440 0844362104 849d0 /
      data x(7) / 0.6087614290 0872063941 6097542898 164d0 /
      data x(8) / 0.5000000000 0000000000 0000000000 000d0 /
      data x(9) / 0.3826834323 6508977172 8459984030 399d0 /
      data x(10) / 0.2588190451 0252076234 8898837624 048d0 /
      data x(11) / 0.1305261922 2005159154 8406227895 489d0 /
c
c           list of major variables
c           -----------------------
c
c           fval   - value of the function f at the points
c                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
c                    k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (bl,br)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (bl,br)
c           res12  - approximation to the integral obtained from cheb12
c           res24  - approximation to the integral obtained from cheb24
c           dqwgts - external function subprogram defining
c                    the four possible weight functions
c           hlgth  - half-length of the interval (bl,br)
c           centr  - mid point of the interval (bl,br)
c
c***first executable statement  dqc25s
      nev = 25
      if(bl.eq.a.and.(alfa.ne.0.0d+00.or.integr.eq.2.or.integr.eq.4))
     * go to 10
      if(br.eq.b.and.(beta.ne.0.0d+00.or.integr.eq.3.or.integr.eq.4))
     * go to 140
c
c           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod
c           scheme.
c
c
      call dqk15w(f,dqwgts,a,b,alfa,beta,integr,bl,br,
     *    result,abserr,resabs,resasc)
      nev = 15
      go to 270
c
c           this part of the program is executed only if a = bl.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
c                  *f(0.5*(br-a)*x+0.5*(br+a))
c
   10 hlgth = 0.5d+00*(br-bl)
      centr = 0.5d+00*(br+bl)
      fix = b-centr
      fval(1) = 0.5d+00*f(hlgth+centr)*(fix-hlgth)**beta
      fval(13) = f(centr)*(fix**beta)
      fval(25) = 0.5d+00*f(centr-hlgth)*(fix+hlgth)**beta
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix-u)**beta
        fval(isym) = f(centr-u)*(fix+u)**beta
   20 continue
      factor = hlgth**(alfa+0.1d+01)
      result = 0.0d+00
      abserr = 0.0d+00
      res12 = 0.0d+00
      res24 = 0.0d+00
      if(integr.gt.2) go to 70
      call dqcheb(x,fval,cheb12,cheb24)
c
c           integr = 1  (or 2)
c
      do 30 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   30 continue
      do 40 i=14,25
        res24 = res24+cheb24(i)*ri(i)
   40 continue
      if(integr.eq.1) go to 130
c
c           integr = 2
c
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
      do 50 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res12+cheb24(i)*rg(i)
   50 continue
      do 60 i=14,25
        res24 = res24+cheb24(i)*rg(i)
   60 continue
      go to 130
c
c           compute the chebyshev series expansion of the
c           following function
c           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
c
   70 fval(1) = fval(1)*dlog(fix-hlgth)
      fval(13) = fval(13)*dlog(fix)
      fval(25) = fval(25)*dlog(fix+hlgth)
      do 80 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*dlog(fix-u)
        fval(isym) = fval(isym)*dlog(fix+u)
   80 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           integr = 3  (or 4)
c
      do 90 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   90 continue
      do 100 i=14,25
        res24 = res24+cheb24(i)*ri(i)
  100 continue
      if(integr.eq.3) go to 130
c
c           integr = 4
c
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
      do 110 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res24+cheb24(i)*rg(i)
  110 continue
      do 120 i=14,25
        res24 = res24+cheb24(i)*rg(i)
  120 continue
  130 result = (result+res24)*factor
      abserr = (abserr+dabs(res24-res12))*factor
      go to 270
c
c           this part of the program is executed only if b = br.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
c                *f(0.5*(b-bl)*x+0.5*(b+bl))
c
  140 hlgth = 0.5d+00*(br-bl)
      centr = 0.5d+00*(br+bl)
      fix = centr-a
      fval(1) = 0.5d+00*f(hlgth+centr)*(fix+hlgth)**alfa
      fval(13) = f(centr)*(fix**alfa)
      fval(25) = 0.5d+00*f(centr-hlgth)*(fix-hlgth)**alfa
      do 150 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix+u)**alfa
        fval(isym) = f(centr-u)*(fix-u)**alfa
  150 continue
      factor = hlgth**(beta+0.1d+01)
      result = 0.0d+00
      abserr = 0.0d+00
      res12 = 0.0d+00
      res24 = 0.0d+00
      if(integr.eq.2.or.integr.eq.4) go to 200
c
c           integr = 1  (or 3)
c
      call dqcheb(x,fval,cheb12,cheb24)
      do 160 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  160 continue
      do 170 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  170 continue
      if(integr.eq.1) go to 260
c
c           integr = 3
c
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
      do 180 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  180 continue
      do 190 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  190 continue
      go to 260
c
c           compute the chebyshev series expansion of the
c           following function
c           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
c
  200 fval(1) = fval(1)*dlog(hlgth+fix)
      fval(13) = fval(13)*dlog(fix)
      fval(25) = fval(25)*dlog(fix-hlgth)
      do 210 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*dlog(u+fix)
        fval(isym) = fval(isym)*dlog(fix-u)
  210 continue
      call dqcheb(x,fval,cheb12,cheb24)
c
c           integr = 2  (or 4)
c
      do 220 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  220 continue
      do 230 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  230 continue
      if(integr.eq.2) go to 260
      dc = dlog(br-bl)
      result = res24*dc
      abserr = dabs((res24-res12)*dc)
      res12 = 0.0d+00
      res24 = 0.0d+00
c
c           integr = 4
c
      do 240 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  240 continue
      do 250 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  250 continue
  260 result = (result+res24)*factor
      abserr = (abserr+dabs(res24-res12))*factor
  270 return
      end
      subroutine dqcheb(x,fval,cheb12,cheb24)
c***begin prologue  dqcheb
c***refer to  dqc25c,dqc25f,dqc25s
c***routines called  (none)
c***revision date  830518   (yymmdd)
c***keywords  chebyshev series expansion, fast fourier transform
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes the chebyshev series expansion
c            of degrees 12 and 24 of a function using a
c            fast fourier transform method
c            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
c            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
c            where t(k,x) is the chebyshev polynomial of degree k.
c***description
c
c        chebyshev series expansion
c        standard fortran subroutine
c        double precision version
c
c        parameters
c          on entry
c           x      - double precision
c                    vector of dimension 11 containing the
c                    values cos(k*pi/24), k = 1, ..., 11
c
c           fval   - double precision
c                    vector of dimension 25 containing the
c                    function values at the points
c                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
c                    where (a,b) is the approximation interval.
c                    fval(1) and fval(25) are divided by two
c                    (these values are destroyed at output).
c
c          on return
c           cheb12 - double precision
c                    vector of dimension 13 containing the
c                    chebyshev coefficients for degree 12
c
c           cheb24 - double precision
c                    vector of dimension 25 containing the
c                    chebyshev coefficients for degree 24
c
c***end prologue  dqcheb
c
      double precision alam,alam1,alam2,cheb12,cheb24,fval,part1,part2,
     *  part3,v,x
      integer i,j
c
      dimension cheb12(13),cheb24(25),fval(25),v(12),x(11)
c
c***first executable statement  dqcheb
      do 10 i=1,12
        j = 26-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   10 continue
      alam1 = v(1)-v(9)
      alam2 = x(6)*(v(3)-v(7)-v(11))
      cheb12(4) = alam1+alam2
      cheb12(10) = alam1-alam2
      alam1 = v(2)-v(8)-v(10)
      alam2 = v(4)-v(6)-v(12)
      alam = x(3)*alam1+x(9)*alam2
      cheb24(4) = cheb12(4)+alam
      cheb24(22) = cheb12(4)-alam
      alam = x(9)*alam1-x(3)*alam2
      cheb24(10) = cheb12(10)+alam
      cheb24(16) = cheb12(10)-alam
      part1 = x(4)*v(5)
      part2 = x(8)*v(9)
      part3 = x(6)*v(7)
      alam1 = v(1)+part1+part2
      alam2 = x(2)*v(3)+part3+x(10)*v(11)
      cheb12(2) = alam1+alam2
      cheb12(12) = alam1-alam2
      alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8)
     *  +x(9)*v(10)+x(11)*v(12)
      cheb24(2) = cheb12(2)+alam
      cheb24(24) = cheb12(2)-alam
      alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8)
     *  +x(3)*v(10)-x(1)*v(12)
      cheb24(12) = cheb12(12)+alam
      cheb24(14) = cheb12(12)-alam
      alam1 = v(1)-part1+part2
      alam2 = x(10)*v(3)-part3+x(2)*v(11)
      cheb12(6) = alam1+alam2
      cheb12(8) = alam1-alam2
      alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6)
     *  -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
      cheb24(6) = cheb12(6)+alam
      cheb24(20) = cheb12(6)-alam
      alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8)
     *  -x(9)*v(10)-x(5)*v(12)
      cheb24(8) = cheb12(8)+alam
      cheb24(18) = cheb12(8)-alam
      do 20 i=1,6
        j = 14-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   20 continue
      alam1 = v(1)+x(8)*v(5)
      alam2 = x(4)*v(3)
      cheb12(3) = alam1+alam2
      cheb12(11) = alam1-alam2
      cheb12(7) = v(1)-v(5)
      alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
      cheb24(3) = cheb12(3)+alam
      cheb24(23) = cheb12(3)-alam
      alam = x(6)*(v(2)-v(4)-v(6))
      cheb24(7) = cheb12(7)+alam
      cheb24(19) = cheb12(7)-alam
      alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
      cheb24(11) = cheb12(11)+alam
      cheb24(15) = cheb12(11)-alam
      do 30 i=1,3
        j = 8-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   30 continue
      cheb12(5) = v(1)+x(8)*v(3)
      cheb12(9) = fval(1)-x(8)*fval(3)
      alam = x(4)*v(2)
      cheb24(5) = cheb12(5)+alam
      cheb24(21) = cheb12(5)-alam
      alam = x(8)*fval(2)-fval(4)
      cheb24(9) = cheb12(9)+alam
      cheb24(17) = cheb12(9)-alam
      cheb12(1) = fval(1)+fval(3)
      alam = fval(2)+fval(4)
      cheb24(1) = cheb12(1)+alam
      cheb24(25) = cheb12(1)-alam
      cheb12(13) = v(1)-v(3)
      cheb24(13) = cheb12(13)
      alam = 0.1d+01/0.6d+01
      do 40 i=2,12
        cheb12(i) = cheb12(i)*alam
   40 continue
      alam = 0.5d+00*alam
      cheb12(1) = cheb12(1)*alam
      cheb12(13) = cheb12(13)*alam
      do 50 i=2,24
        cheb24(i) = cheb24(i)*alam
   50 continue
      cheb24(1) = 0.5d+00*alam*cheb24(1)
      cheb24(25) = 0.5d+00*alam*cheb24(25)
      return
      end
      subroutine dqk15w(f,w,p1,p2,p3,p4,kp,a,b,result,abserr,
     *   resabs,resasc)
c***begin prologue  dqk15w
c***date written   810101   (yymmdd)
c***revision date  830518   (mmddyy)
c***category no.  h2a2a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b), with error
c                           estimate
c                       j = integral of abs(f*w) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           double precision version
c
c           parameters
c             on entry
c              f      - double precision
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              w      - double precision
c                       function subprogram defining the integrand
c                       weight function w(x). the actual name for w
c                       needs to be declared e x t e r n a l in the
c                       calling program.
c
c              p1, p2, p3, p4 - double precision
c                       parameters in the weight function
c
c              kp     - integer
c                       key for indicating the type of weight function
c
c              a      - double precision
c                       lower limit of integration
c
c              b      - double precision
c                       upper limit of integration
c
c            on return
c              result - double precision
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule (resg).
c
c              abserr - double precision
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - double precision
c                       approximation to the integral of abs(f)
c
c              resasc - double precision
c                       approximation to the integral of abs(f-i/(b-a))
c
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15w
c
      double precision a,absc,absc1,absc2,abserr,b,centr,dabs,dhlgth,
     *  dmax1,dmin1,d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  p1,p2,p3,p4,resabs,resasc,resg,resk,reskh,result,uflow,w,wg,wgk,
     *  xgk
      integer j,jtw,jtwm1,kp
      external f,w
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(4)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point gauss-kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ... abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point gauss-kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/
     *     0.9914553711208126d+00,     0.9491079123427585d+00,
     *     0.8648644233597691d+00,     0.7415311855993944d+00,
     *     0.5860872354676911d+00,     0.4058451513773972d+00,
     *     0.2077849550078985d+00,     0.0000000000000000d+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/
     *     0.2293532201052922d-01,     0.6309209262997855d-01,
     *     0.1047900103222502d+00,     0.1406532597155259d+00,
     *     0.1690047266392679d+00,     0.1903505780647854d+00,
     *     0.2044329400752989d+00,     0.2094821410847278d+00/
c
      data wg(1),wg(2),wg(3),wg(4)/
     *     0.1294849661688697d+00,    0.2797053914892767d+00,
     *     0.3818300505051189d+00,    0.4179591836734694d+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f*w over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15w
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 15-point kronrod approximation to the
c           integral, and estimate the error.
c
      fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j=1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*
     *  0.5d+02)*resabs,abserr)
      return
      end
      subroutine dqmomo(alfa,beta,ri,rj,rg,rh,integr)
c***begin prologue  dqmomo
c***date written   820101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,c3a2
c***keywords  modified chebyshev moments
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes modified chebsyshev moments. the k-th
c            modified chebyshev moment is defined as the integral over
c            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
c            polynomial of degree k.
c***description
c
c        modified chebyshev moments
c        standard fortran subroutine
c        double precision version
c
c        parameters
c           alfa   - double precision
c                    parameter in the weight function w(x), alfa.gt.(-1)
c
c           beta   - double precision
c                    parameter in the weight function w(x), beta.gt.(-1)
c
c           ri     - double precision
c                    vector of dimension 25
c                    ri(k) is the integral over (-1,1) of
c                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
c
c           rj     - double precision
c                    vector of dimension 25
c                    rj(k) is the integral over (-1,1) of
c                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
c
c           rg     - double precision
c                    vector of dimension 25
c                    rg(k) is the integral over (-1,1) of
c                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           rh     - double precision
c                    vector of dimension 25
c                    rh(k) is the integral over (-1,1) of
c                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
c
c           integr - integer
c                    input parameter indicating the modified
c                    moments to be computed
c                    integr = 1 compute ri, rj
c                           = 2 compute ri, rj, rg
c                           = 3 compute ri, rj, rh
c                           = 4 compute ri, rj, rg, rh
c
c***references  (none)
c***routines called  (none)
c***end prologue  dqmomo
c
      double precision alfa,alfp1,alfp2,an,anm1,beta,betp1,betp2,ralf,
     *  rbet,rg,rh,ri,rj
      integer i,im1,integr
c
      dimension rg(25),rh(25),ri(25),rj(25)
c
c
c***first executable statement  dqmomo
      alfp1 = alfa+0.1d+01
      betp1 = beta+0.1d+01
      alfp2 = alfa+0.2d+01
      betp2 = beta+0.2d+01
      ralf = 0.2d+01**alfp1
      rbet = 0.2d+01**betp1
c
c           compute ri, rj using a forward recurrence relation.
c
      ri(1) = ralf/alfp1
      rj(1) = rbet/betp1
      ri(2) = ri(1)*alfa/alfp2
      rj(2) = rj(1)*beta/betp2
      an = 0.2d+01
      anm1 = 0.1d+01
      do 20 i=3,25
        ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
        rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1d+01
   20 continue
      if(integr.eq.1) go to 70
      if(integr.eq.3) go to 40
c
c           compute rg using a forward recurrence relation.
c
      rg(1) = -ri(1)/alfp1
      rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
      an = 0.2d+01
      anm1 = 0.1d+01
      im1 = 2
      do 30 i=3,25
        rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/
     *  (anm1*(an+alfp1))
        anm1 = an
        an = an+0.1d+01
        im1 = i
   30 continue
      if(integr.eq.2) go to 70
c
c           compute rh using a forward recurrence relation.
c
   40 rh(1) = -rj(1)/betp1
      rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
      an = 0.2d+01
      anm1 = 0.1d+01
      im1 = 2
      do 50 i=3,25
        rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+
     *  anm1*rj(i))/(anm1*(an+betp1))
        anm1 = an
        an = an+0.1d+01
        im1 = i
   50 continue
      do 60 i=2,25,2
        rh(i) = -rh(i)
   60 continue
   70 do 80 i=2,25,2
        rj(i) = -rj(i)
   80 continue
   90 return
      end
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c***begin prologue  dqpsrt
c***refer to  dqage,dqagie,dqagpe,dqawse
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  sequential sorting
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine maintains the descending ordering in the
c            list of the local error estimated resulting from the
c            interval subdivision process. at each call two error
c            estimates are inserted using the sequential search
c            method, top-down for the largest error estimate and
c            bottom-up for the smallest error estimate.
c***description
c
c           ordering routine
c           standard fortran subroutine
c           double precision version
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - double precision
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - double precision
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k elements
c                       of which contain pointers to the error
c                       estimates, such that
c                       elist(iord(1)),...,  elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c
c***end prologue  dqpsrt
c
      double precision elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     *  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  dqpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed if, due to a
c           difficult integrand, subdivision increased the error
c           estimate. in the normal case the insert procedure should
c           start after the nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to be maintained
c           in descending order. this number depends on the number of
c           subdivisions still allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
      double precision function dqwgts(x,a,b,alfa,beta,integr)
c***begin prologue  dqwgts
c***refer to dqk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, algebraico-logarithmic
c             end-point singularities
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine dqaws and defines the weight function.
c***end prologue  dqwgts
c
      double precision a,alfa,b,beta,bmx,dlog,x,xma
      integer integr
c***first executable statement  dqwgts
      xma = x-a
      bmx = b-x
      dqwgts = xma**alfa*bmx**beta
      go to (40,10,20,30),integr
   10 dqwgts = dqwgts*dlog(xma)
      go to 40
   20 dqwgts = dqwgts*dlog(bmx)
      go to 40
   30 dqwgts = dqwgts*dlog(xma)*dlog(bmx)
   40 return
      end
      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
c    D1MACH ( 3) = B^(-T), the smallest relative spacing.
c    D1MACH ( 4) = B^(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer i

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      else if ( i == 1 ) then
        d1mach = 4.450147717014403D-308
      else if ( i == 2 ) then
        d1mach = 8.988465674311579D+307
      else if ( i == 3 ) then
        d1mach = 1.110223024625157D-016
      else if ( i == 4 ) then
        d1mach = 2.220446049250313D-016
      else if ( i == 5 ) then
        d1mach = 0.301029995663981D+000
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      end if

      return
      end
      function i1mach ( i )

c*********************************************************************72
c
cc I1MACH returns integer machine dependent constants.
c
c  Discussion:
c
c    Input/output unit numbers.
c
c      I1MACH(1) = the standard input unit.
c      I1MACH(2) = the standard output unit.
c      I1MACH(3) = the standard punch unit.
c      I1MACH(4) = the standard error message unit.
c
c    Words.
c
c      I1MACH(5) = the number of bits per integer storage unit.
c      I1MACH(6) = the number of characters per integer storage unit.
c
c    Integers.
c
c    Assume integers are represented in the S digit base A form:
c
c      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
c
c    where 0 <= X(1:S-1) < A.
c
c      I1MACH(7) = A, the base.
c      I1MACH(8) = S, the number of base A digits.
c      I1MACH(9) = A^S-1, the largest integer.
c
c    Floating point numbers
c
c    Assume floating point numbers are represented in the T digit 
c    base B form:
c
c      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
c
c    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
c
c      I1MACH(10) = B, the base.
c
c    Single precision
c
c      I1MACH(11) = T, the number of base B digits.
c      I1MACH(12) = EMIN, the smallest exponent E.
c      I1MACH(13) = EMAX, the largest exponent E.
c
c    Double precision
c
c      I1MACH(14) = T, the number of base B digits.
c      I1MACH(15) = EMIN, the smallest exponent E.
c      I1MACH(16) = EMAX, the largest exponent E.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 16.
c
c    Output, integer I1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      integer i1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      else if ( i == 1 ) then
        i1mach = 5
      else if ( i == 2 ) then
        i1mach = 6
      else if ( i == 3 ) then
        i1mach = 7
      else if ( i == 4 ) then
        i1mach = 6
      else if ( i == 5 ) then
        i1mach = 32
      else if ( i == 6 ) then
        i1mach = 4
      else if ( i == 7 ) then
        i1mach = 2
      else if ( i == 8 ) then
        i1mach = 31
      else if ( i == 9 ) then
        i1mach = 2147483647
      else if ( i == 10 ) then
        i1mach = 2
      else if ( i == 11 ) then
        i1mach = 24
      else if ( i == 12 ) then
        i1mach = -125
      else if ( i == 13 ) then
        i1mach = 128
      else if ( i == 14 ) then
        i1mach = 53
      else if ( i == 15 ) then
        i1mach = -1021
      else if ( i == 16 ) then
        i1mach = 1024
      else if ( 16 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      end if

      return
      end