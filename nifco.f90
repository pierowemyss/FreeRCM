! Author: Piero Wemyss
! Date: July 8, 2024
! 
! __________________________________________________________________________ 
!                  Non-Ideal Fortran-computed COefficients
! __________________________________________________________________________ 
!    Fortran subroutines for computing coefficients for non-ideal mixtures
! --------------------------------------------------------------------------
! 
! HYSYS-modified NRTL | INPUT:  liquid mole fractions  |  OUTPUT: activity coefficients
!                     |         temperature (celsius)  |
!                     |         binary coefficients    |
!                     |         number of components   |
subroutine NRTL(x, Tcel, a, b, c, NComps, gamma)
    implicit none

    integer, intent(in) :: NComps 
    real(8), intent(in), dimension(0:NComps-1) :: x
    real(8), intent(in), dimension(0:NComps-1,0:NComps-1) :: a,b,c
    real(8), intent(in) :: Tcel
    real(8) :: tauji, Gji, tauki, Gki, term1, tauij, Gij, taukj, Gkj, taumj, Gmj, lngam
    real(8), dimension(0:NComps-1) :: term1j, term1k, term2k, term2m, term2j
    real(8), intent(out), dimension(0:NComps-1) :: gamma
    integer :: i, j, k, m
    real(8) :: T
    real(4), parameter :: R = 8.314

    T = Tcel + 273.15

    do i = 0, NComps-1 

        do j = 0, NComps-1 
            tauji = (a(j,i)+b(j,i)/T)
            Gji = exp(-c(j,i)*tauji)
            term1j(j) = x(j)*tauji*Gji
        enddo
        do k = 0, NComps-1 
            tauki = (a(k,i)+b(k,i)/T)
            Gki = exp(-c(k,i)*tauki)
            term1k(k) = x(k)*Gki
        end do

        term1 = sum(term1j)/sum(term1k)

        do j = 0, NComps-1 
            tauij = (a(i,j)+b(i,j)/T)
            Gij = exp(-c(i,j)*tauij)

            do k = 0, NComps-1 
                taukj = (a(k,j)+b(k,j)/T)
                Gkj = exp(-c(k,j)*taukj)
                term2k(k) = x(k)*Gkj
            enddo
            do m = 0, NComps-1 
                taumj = (a(m,j)+b(m,j)/T)
                Gmj = exp(-c(m,j)*taumj)
                term2m(m) = x(m)*taumj*Gmj
            enddo

            term2j(j) = (x(j)*Gij/(sum(term2k)))*(tauij-(sum(term2m)/sum(term2k)))
        enddo

        lngam = term1 + sum(term2j)

        gamma(i) = exp(lngam)
    enddo
end subroutine

!
! SRK EOS  | INPUT: vapor mole fractions  | OUTPUT: fugacity coefficient
!          |        temperature (celsius) |
!          |        pressure (bar)        |
!          |        crit temp (celsius)   |
!          |        crit press (bar)      |
!          |        acentricity           |
!          |        number of components  |
subroutine SRK(x,T,P,TcCel,Pc,omega,NComps,phi)
    
    use minpack_module, only: wp, hybrd1, dpmpar, enorm
    use iso_fortran_env, only: nwrite => output_unit
    
    implicit none

    common /shared_data/ beta, q

    integer, intent(in) :: NComps
    real(4), intent(in) :: P, T
    real(8), intent(in), dimension(0:NComps-1) :: x, TcCel, Pc, omega
    real(8), intent(out), dimension(0:NComps-1) :: phi
    real(8) :: aij, am, bm, beta, q, A, B, Z0
    real(8), dimension(0:NComps-1) :: Tc, Tr, Pr, alph, ai, bi, bvec, ajvec, aivec, AAi, BBi, logphi
    integer :: i, j, info

    integer,parameter :: n = 1
    integer,parameter :: lwa = (n*(3*n+13))/2

    external :: RKfunc
    real(8) :: Z(n), fvec(n), wa(lwa)
    real(8) :: tol

    real(4), parameter :: R = 8.14
    
    Tc = TcCel + 273.15
    Tr = (T+273.15)/Tc
    Pr = P/Pc
    
    alph = (1.0+(0.480+1.574*omega-0.176*omega**(2.0))*(1.0-Tr**(0.5)))**2.0
    ai = 0.42748*(alph*R**2*Tc**2)/Pc
    bi = 0.08664*R*Tc/Pc
    
    do i = 0, NComps-1
    
        bvec(i) = x(i)*bi(i)
    
        do j = 0, NComps-1
    
            aij = (ai(i)*ai(j))**(0.5)
            ajvec(j) = x(i)*x(j)*aij
    
        enddo
    
        aivec(i) = sum(ajvec)
    
    enddo
    
    am = sum(aivec)
    bm = sum(bvec)
    
    beta = bm*P/(R*(T+273.15))
    q = am/(bm*R*(T+273.15))
    
    Z0 = 1.0d0
    Z(1) = Z0

    tol = 1.0d-8

    call hybrd1(RKfunc, n, Z, fvec, tol, info, wa, lwa)

    if (info /= 1) then
        write(*, '(A, i2)') "NIFCO.SRK ERROR: Unsuccessful root find for compressibility."
    end if

    AAi = (ai/(R**(2.0)*(T+273.15)**(2.5)))**(0.5)
    A = (am/(R**(2.0)*(T+273.15)**(2.5)))**(0.5)
    
    BBi = bi/(R*(T+273.15))
    B = bm/(R*(T+273.15))

    do i = 0, NComps-1
        logphi(i) =  0.4343*(Z(1) - 1.0)*(BBi(i)/B) - log10(Z(1) - B*P) - (A**(2.0)/B)*((2.0*AAi(i)/A - BBi(i)/B))*log10(1.0 + B*P/Z(1))
    end do
    
    phi = 10**logphi

end subroutine

!
! SRK Z-fact  | INPUT: vapor mole fractions  | OUTPUT: compressibility factor
!             |        temperature (celsius) |
!             |        pressure (bar)        |
!             |        crit temp (celsius)   |
!             |        crit press (bar)      |
!             |        acentricity           |
!             |        number of components  |
subroutine Zfact(x,T,P,TcCel,Pc,omega,NComps,Z)
    
    use minpack_module, only: wp, hybrd1, dpmpar, enorm
    use iso_fortran_env, only: nwrite => output_unit
    
    implicit none

    common /shared_data/ beta, q

    integer, intent(in) :: NComps
    real(4), intent(in) :: P, T
    real(8), intent(in), dimension(0:NComps-1) :: x, TcCel, Pc, omega
    real(8) :: aij, am, bm, beta, q, Z0
    real(8), dimension(0:NComps-1) :: Tc, Tr, Pr, alph, ai, bi, bvec, ajvec, aivec
    integer :: i, j, info

    integer,parameter :: n = 1
    integer,parameter :: lwa = (n*(3*n+13))/2

    external :: RKfunc
    real(8) :: fvec(n), wa(lwa)
    real(8), intent(out) :: Z(n)
    real(8) :: tol

    real(4), parameter :: R = 8.14
    
    Tc = TcCel + 273.15
    Tr = (T+273.15)/Tc
    Pr = P/Pc
    
    alph = (1.0+(0.480+1.574*omega-0.176*omega**(2.0))*(1.0-Tr**(0.5)))**2.0
    ai = 0.42748*(alph*R**2*Tc**2)/Pc
    bi = 0.08664*R*Tc/Pc
    
    do i = 0, NComps-1
    
        bvec(i) = x(i)*bi(i)
    
        do j = 0, NComps-1
    
            aij = (ai(i)*ai(j))**(0.5)
            ajvec(j) = x(i)*x(j)*aij
    
        enddo
    
        aivec(i) = sum(ajvec)
    
    enddo
    
    am = sum(aivec)
    bm = sum(bvec)
    
    beta = bm*P/(R*(T+273.15))
    q = am/(bm*R*(T+273.15))
    
    Z0 = 1.0d0
    Z(1) = Z0

    tol = 1.0d-8

    call hybrd1(RKfunc, n, Z, fvec, tol, info, wa, lwa)

    if (info /= 1) then
        write(*, '(A, i2)') "NIFCO.SRK ERROR: Unsuccessful root find for compressibility."
    end if

end subroutine

subroutine RKfunc(n, Z, fvec, iflag)
    
    implicit none

    common /shared_data/ beta, q

    integer, intent(in) :: n
    real(8), intent(in) :: Z(n)
    real(8), intent(out) :: fvec(n)
    integer, intent(inout) :: iflag
!f2py intent(inout) :: iflag
    real(8) :: beta, q

    fvec(1) = 1.0d0 + beta - (q*beta*(Z(1)-beta))/(Z(1)+beta) - Z(1)

end subroutine
