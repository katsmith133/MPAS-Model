module ocn_adc_mixing_fused

  use mpas_kind_types
  use ocn_constants
  use ocn_turbulence
  use mpas_log
  use ocn_config
  use mpas_constants

  implicit none

  logical :: defineFirst, stopflag

  integer :: i1,i2

  contains

  subroutine swap_time_levels

    i1 = mod(iterCount - 1, 2) + 1
    i2 = mod(iterCount, 2) + 1

  end subroutine swap_time_levels

  subroutine dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
		integer,intent(in) :: nVertLevels, nCells
    real (kind=RKIND), dimension(2,nVertLevels,nCells), intent(in) :: activeTracers
    real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: alphaT, betaS
    integer :: iCell, k, ij, i

    real (kind=RKIND), dimension(nVertLevels) :: B, Bup, Bdn
    real (kind=RKIND), dimension(nVertLevels+1) :: KEsps, tke, BupEdge, BdnEdge, BEdge
    real (kind=RKIND) :: sigav, tumdav, sumdav, Tup, Tdn, Sup, Sdn
    real (kind=RKIND) :: sumv, sumv0
    real (kind=RKIND), parameter :: refT = 15.0_RKIND, refS = 35.0_RKIND, minlen = 0.55_RKIND

    !NOTE: will need to convert to some form of displaced density in the mpas framework soon
    !possibly go back to a more traditional length scale formulation
    do iCell = 1,nCells
       do k = 1,nVertLevels+1
          KEsps(k) = areaFraction(k,iCell) * KspsD(i1,k,iCell) + &
                     (1.0_RKIND - areaFraction(k,iCell)) * KspsU(i1,k,iCell)
          tke(k) = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell)) ! + KEsps(k)
       enddo

       do k = 1,nVertLevels
          B(k) = gravity * (alphaT(k,iCell) * (activeTracers(1,k,iCell) - refT) - &
                            betaS(k,iCell)  * (activeTracers(2,k,iCell) - refS))

          sigav  = 0.5_RKIND * (areaFraction(k,iCell) + areaFraction(k+1,iCell))
          tumdav = 0.5_RKIND * (tumd(k,iCell) + tumd(k+1,iCell))
          sumdav = 0.5_RKIND * (sumd(k,iCell) + sumd(k+1,iCell))

          Tup = activeTracers(1,k,iCell) + sigav*tumdav
          Tdn = activeTracers(1,k,iCell) - (1.0_RKIND - sigav)*tumdav
          Sup = activeTracers(2,k,iCell) + sigav*sumdav
          Sdn = activeTracers(2,k,iCell) - (1.0_RKIND - sigav)*sumdav

          Bup(k) = gravity * (alphaT(k,iCell) * (Tup - refT) - betaS(k,iCell) * (Sup - refS))
          Bdn(k) = gravity * (alphaT(k,iCell) * (Tdn - refT) - betaS(k,iCell) * (Sdn - refS))

          if(k>1) THEN
             BupEdge(k) = 0.5_RKIND * (Bup(k-1) + Bup(k))
             BdnEdge(k) = 0.5_RKIND * (Bdn(k-1) + Bdn(k))
             BEdge(k)   = 0.5_RKIND * (B(k-1) + B(k))
          endif
       enddo

       BdnEdge(nVertLevels+1) = BdnEdge(nVertLevels)
       BupEdge(nVertLevels+1) = BupEdge(nVertLevels)
       BEdge(nVertLevels+1)   = BEdge(nVertLevels)

       BdnEdge(1) = BdnEdge(2)
       BupEdge(1) = BupEdge(2)
       BEdge(1)   = BEdge(2)

       do k = 2,nVertLevels
          ! updraft length scale
          sumv  = 0.0_RKIND
          sumv0 = 0.0_RKIND
          ij=k
          lenup(k,iCell) = 0.0_RKIND
          do while(sumv0 <= tke(k) .and. ij > 1)
             sumv = sumv0 + (BEdge(ij) - BupEdge(k)) * (ze(ij-1,iCell)-ze(ij,iCell))
             if (sumv > tke(k)) THEN
                lenup(k,iCell) = max(minlen, &
                                     lenup(k,iCell) + &
                                     abs(ze(ij-1,iCell)-ze(ij,iCell))*(tke(k)-sumv0)/(sumv-sumv0))
                exit
             else
                lenup(k,iCell) = lenup(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
             endif
             sumv0 = sumv
             ij = ij - 1
          end do

          ! downdraft length scale
          sumv  = 0.0_RKIND
          sumv0 = 0.0_RKIND
          ij=k
          lendn(k,iCell) = 0.0_RKIND
          do while(sumv0 <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv0 + (BdnEdge(k) - BEdge(ij)) * (ze(ij-1,iCell)-ze(ij,iCell))
             if (sumv > tke(k)) THEN
                lendn(k,iCell) = max(minlen, &
                                     lendn(k,iCell) + &
                                     abs(ze(ij-1,iCell)-ze(ij,iCell))*(tke(k)-sumv0)/(sumv-sumv0))
                exit
             else
                lendn(k,iCell) = lendn(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
             endif
             sumv0 = sumv
             ij = ij + 1
          end do

          length(k,iCell) = 2.0_RKIND * lenup(k,iCell) * lendn(k,iCell)  &
                            / (lenup(k,iCell) + lendn(k,iCell))
          ! write(*,*) k, iCell, lenup(k,iCell), lendn(k,iCell)
       enddo
       ! write(*,*) '----'
    enddo

    length(1,iCell) = 0.55_RKIND
    length(nVertLevels+1,iCell) = 0.55_RKIND

  end subroutine dissipation_lengths2

  subroutine compute_ADC_tends(nCells,nVertLevels, nTracers, dt,activeTracers, uvel, vvel, BVF,  &
                  uwsfc, vwsfc, wtsfc, wssfc, alphaT, betaS, fCell, boundaryLayerDepth)

      integer,intent(in) :: nCells, nVertLevels, nTracers
      real,intent(in) :: dt

      real,dimension(nTracers,nVertLevels,nCells),intent(inout) :: activeTracers
      real,dimension(nVertLevels,nCells),intent(inout) :: uvel, vvel, alphaT, betaS
      real,dimension(nCells),intent(in) :: uwsfc, vwsfc, wtsfc, wssfc, fCell
      real,dimension(nCells),intent(inout) :: boundaryLayerDepth
      real,dimension(nVertLevels,nCells),intent(inout) :: BVF
      integer :: niter, iIter,iCell,k, supercycle

      real :: dt_small

      real,dimension(nVertLevels+1,nCells) :: w2tend
      real,dimension(nVertLevels,nCells) :: w3tend

      real :: Sw, St, Ss, Eav, Dav, sigav, sigavp1, wumdAv, tumdAv, sumdAv, wumdAvp1, tumdAvp1, sumdAvp1
      real :: Swup, KspsUav, KspsDav, KspsUavp1, KspsDavp1, KE, Mcav, lenav,u2av,v2av,w2av
      real :: w3temp, w3check, mval, KEsps, Uz, Vz, dz

      real :: invLen, l, len1, len2, lenmax, integrandTop, integrandBot
      real :: len0, len2_1, len2_2, sfcBuoy, lengthT, bvfT
      real :: KEm1, KEp1, tauUP, tauDN, tomUP, tomDN
      real :: tauM1, tau, tauP1, tauAV, utemp, vtemp
      real :: B, Cval, diff, wtav, dzmid, Ksps, Sz, Tz, w4k, w4kp1, w2k, w2kp1
      real :: lareaFraction, wstar, Q, w3av, tempMoment
      real :: sfcFrictionVelocitySquared, wtSumUp, wtSumDn, wsSumUp, wsSumDn

      real,dimension(nVertLevels,nCells) :: Swumd
      real,dimension(nVertLevels,nCells) :: tauw3, tauTemp, tauSalt, tauVel, tauvVel
      real,dimension(nVertLevels+1,nCells) :: wttend, wstend, uwtend, vwtend, u2tend, v2tend, &
        uttend, vttend, ustend, vstend, uvtend
      real,dimension(nVertLevels,nCells) :: epstend, areaFractionMid
      real,dimension(nVertLevels+1,nCells) :: KspsUtend, KspsDtend

      dt_small = config_adc_timestep
      niter = dt / dt_small

      call swap_time_levels

      !on further examination build_diagnostics array can live outside the iter loop
      do iCell=1,nCells
        Q = grav*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))* &
            boundaryLayerDepth(iCell)
        if(Q > 0) then
          wstar = abs(Q)**(1.0/3.0)
        else
          wstar = 0.0
        endif

        tumd(1,iCell) = 0.0_RKIND
        wumd(1,iCell) = 0.0_RKIND
        areaFraction(1,iCell) = 0.5_RKIND
        Mc(1,iCell) = 0.0_RKIND
        w2t(1,iCell) = -0.3_RKIND*wstar * wtsfc(iCell)
        w2s(1,iCell) = 0.3_RKIND*wstar * wssfc(iCell)

        sfcFrictionVelocitySquared = sqrt(uwsfc(iCell)**2 + vwsfc(iCell)**2)
        do k=1,2
          u2(k,1,iCell) = 0.0_RKIND! 2.0*sfcFrictionVelocitySquared + 0.3*wstar**2.0
          v2(k,1,iCell) = 0.0_RKIND!2.0*sfcFrictionVelocitySquared + 0.3*wstar**2.0
          uw(k,1,iCell) = -uwsfc(iCell)
          vw(k,1,iCell) = -vwsfc(iCell)
          wt(k,1,iCell) = wtsfc(iCell)
          ws(k,1,iCell) = wssfc(iCell)
          KE = 0.5*(u2(i1,1,iCell) + v2(i1,1,iCell))
          eps(k,1,iCell) = KE**1.5/(0.5*(ze(1,iCell) - ze(2,iCell) + 1.0E-10_RKIND))
        enddo
      enddo

      do iIter=1,niter
        !Kernel 1 inlined versions of the base arrays, needed for later to make them collapsible
        do iCell=1,nCells
          do k=2,nVertLevels
            w3av = 0.5_RKIND*(w3(i1,k-1,iCell) + w3(i1,k,iCell))

            Sw = w3av / (w2(i1,k,iCell)**1.5_RKIND + 1.0E-15_RKIND)
            lareaFraction = 0.5_RKIND + 0.5_RKIND*Sw / sqrt(4.0_RKIND + Sw**2)

            if(lareaFraction < 0.01_RKIND) lareaFraction = 0.01_RKIND
            if(lareaFraction > 0.99_RKIND) lareaFraction = 0.99_RKIND

            areaFraction(k,iCell) = lareaFraction


            wumd(k,iCell) = sqrt(w2(i1,k,iCell) / (areaFraction(k,iCell) * &
               (1.0_RKIND - areaFraction(k,iCell))))
            Mc(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - &
               areaFraction(k,iCell)) * wumd(k,iCell)

            tumd(k,iCell) = wt(i1,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))
            sumd(k,iCell) = ws(i1,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))

!            if(w2(i1,k,iCell) <= epsilon + 1.0e-12) then
!               areaFraction(k,iCell) = 0.5_RKIND
!               wumd(k,iCell) = 0.0_RKIND
!               tumd(k,iCell) = 0.0_RKIND
!               sumd(k,iCell) = 0.0_RKIND
!            endif
         enddo
        enddo

        do iCell = 1,nCells
          do k=2,nVertLevels
            sigav = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
            tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
            sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
            wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
            w2t(k,iCell) = -sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*tumdav
            w2s(k,iCell) = -sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*sumdav
            !also use this loop to reset the cliptends for the step
            u2cliptend(k,iCell) = 0.0_RKIND
            v2cliptend(k,iCell) = 0.0_RKIND
            w2cliptend(k,iCell) = 0.0_RKIND
          enddo
        enddo
          if(config_adc_use_old_length_scale)then
					  call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
          endif


        !This cell loop computes w3Tend and TOMs for later tendencies.  Should be okay to collapse loops here
        do iCell = 1,nCells
            do k=1,nVertLevels
              Eav = 0.5*(Entrainment(k+1,iCell) + Entrainment(k,iCell))
              Dav = 0.5*(Detrainment(k+1,iCell) + Detrainment(k,iCell))
              u2av = 0.5*(u2(i1,k,iCell) + u2(i1,k+1,iCell))
              v2av = 0.5*(v2(i1,k,iCell) + v2(i1,k+1,iCell))
              w2av = 0.5*(w2(i1,k,iCell) + w2(i1,k+1,iCell))

              sigav = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
              wumdav = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
              tumdav = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
              sumdav = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
              Mcav = 0.5*(Mc(k,iCell) + Mc(k+1,iCell))
              lenav = 0.5*(length(k,iCell) + length(k+1,iCell))
              KEsps = sigav*KspsDav+ (1.0 - sigav)*KspsUav
              KE = sqrt((u2av+v2av+w2av) + 0.0*KEsps)

              dz = ze(k,iCell) - ze(k+1,iCell)
              !KE = sqrt(areaFraction(k,iCell)*KspsUav + (1.0 - areaFraction(k,iCell))*KspsDav)
              Swumd(k,iCell) = - 2.0/3.0*(1.0_RKIND/sigav*(areaFraction(k,iCell)*KspsD(i1,k,iCell) - &
                areaFraction(k+1,iCell)*KspsD(i1,k+1,iCell)) / dz - 1.0_RKIND / (1.0_RKIND - &
                sigav)*((1.0_RKIND - areaFraction(k,iCell))*KspsU(i1,k,iCell) - (1.0_RKIND - &
                areaFraction(k+1,iCell))*KspsU(i1,k+1,iCell)) / dz)

              tauw3(k,iCell) = C_mom_w3*KE / (1.0E-15_RKIND + sqrt(2.0_RKIND)*lenAv)

              w3tend1(k,iCell) = -wumdav**3.0*(Eav*(3.0*sigav - 2.0) + Dav*(3.0*sigav - 1.0))
              w3tend2(k,iCell) = -wumdav**3.0*(6.0*sigav**2.0 - 6.0*sigav + 1)*             &
                    (areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))*  &
                    wumd(k,iCell) - areaFraction(k+1,iCell)*(1.0_RKIND - &
                    areaFraction(k+1,iCell))*wumd(k+1,iCell))/dz
              w3tend3(k,iCell) = -1.5_RKIND*sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)* &
                    wumdav**2.0_RKIND*((1.0_RKIND - 2.0_RKIND*areaFraction(k,iCell))*      &
                    wumd(k,iCell)**2.0_RKIND - (1.0_RKIND -                                &
                    2.0_RKIND*areaFraction(k+1,iCell))*wumd(k+1,iCell)**2.0_RKIND) / dz
              w3tend4(k,iCell) = -3.0_RKIND*(1.0_RKIND - 2.0_RKIND*sigav)*Mcav*wumdav*Swumd(k,iCell)
              w3tend5(k,iCell) = 3.0_RKIND*(1.0_RKIND - c11)*grav*(alphaT(k,iCell)*w2t(k,iCell) - &
                    betaS(k,iCell)*w2s(k,iCell))

              w3tend(k,iCell) = w3tend1(k,iCell) + w3tend2(k,iCell) + w3tend3(k,iCell) + &
                w3tend4(k,iCell) + w3tend5(k,iCell)

              if(k>1 .and. k < nVertLevels .and. kappa_w3 > 0.0) then
                 w3tend(k,iCell) = w3tend(k,iCell) + kappa_w3*(w3(i1,k-1,iCell) &
                    - w3(i1,k+1,iCell)) / (zm(k-1,iCell) - zm(k+1,iCell))**2.0
              endif

            ! now get all the downgradient TOMs
              KE = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
              KEp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
              lenav = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))
              diff = C_mom * sqrt(0.5_RKIND*(KE + KEp1)) * lenav
              dz = ze(k,iCell) - ze(k+1,iCell)
              uw2(k,iCell) = -diff*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) / dz
              vw2(k,iCell) = -diff*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) / dz
              u2w(k,iCell) = -diff*(u2(i1,k,iCell) - u2(i1,k+1,iCell)) / dz
              v2w(k,iCell) = -diff*(v2(i1,k,iCell) - v2(i1,k+1,iCell)) / dz
              uvw(k,iCell) = -diff*(uv(i1,k,iCell) - uv(i1,k+1,iCell)) / dz

              diff = C_therm*sqrt(0.5*(KE + KEp1)) * lenav
              uwt(k,iCell) = -diff*(ut(i1,k,iCell) - ut(i1,k+1,iCell)) / dz
              vwt(k,iCell) = -diff*(vt(i1,k,iCell) - vt(i1,k+1,iCell)) / dz
              uws(k,iCell) = -diff*(us(i1,k,iCell) - us(i1,k+1,iCell)) / dz
              vws(k,iCell) = -diff*(vs(i1,k,iCell) - vs(i1,k+1,iCell)) / dz

            enddo !nVertLevels
          enddo !nCells

!          if(config_adc_use_old_length_scale)then
!					  call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
!          endif

          !next all second order moment tendencies also apply tendencies
          do iCell=1,nCells
            do k=2,nVertLevels
              dzmid = (zm(k-1,iCell) - zm(k,iCell))
              dz = ze(k,iCell) - ze(k+1,iCell)
              B = grav*(alphaT(k,iCell)*areaFraction(k,iCell)*(1.0_RKIND - &
                areaFraction(k,iCell))*wumd(k,iCell)*tumd(k,iCell) -      &
                betaS(k,iCell)*areaFraction(k,iCell)*(1.0_RKIND -          &
                areaFraction(k,iCell))*wumd(k,iCell)*sumd(k,iCell))

              Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dzmid
              Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dzmid
              Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dzmid
              Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dzmid

              KE = sqrt(0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + areaFraction(k,iCell)* &
                        (1.0_RKIND - areaFraction(k,iCell))*wumd(k,iCell)))

              tauTemp(k,iCell) = c_Temp*KE / (1.0E-15_RKIND + length(k,iCell))
              tauSalt(k,iCell) = c_Salt*KE / (1.0E-15_RKIND + length(k,iCell))
              tauVel(k,iCell) = c_Vel*KE / (1.0E-15_RKIND + length(k,iCell))
              tauvVel(k,iCell) = c_vVel*KE / (1.0E-15_RKIND + length(k,iCell))

              w2tend1(k,iCell) = -wumd(k,iCell)**2.0_RKIND*(         &
                Entrainment(k,iCell) + Detrainment(k,iCell))
              w2tend2(k,iCell) = (Mc(k-1,iCell)*(1.0_RKIND - 2.0_RKIND*    &
                areaFraction(k-1,iCell))*wumd(k-1,iCell)**2.0 - Mc(k+1,iCell)* &
                (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell))*wumd(k+1,iCell)**2.0) &
                / (ze(k-1,iCell) - ze(k+1,iCell))
              w2tend2(k,iCell) = -(w3(i1,k-1,iCell) - w3(i1,k,iCell)) / dzmid
              w2tend3(k,iCell) = tauvVel(k,iCell)*(u2(i1,k,iCell) + v2(i1,k,iCell))/3.0_RKIND
              w2tend4(k,iCell) = (2.0_RKIND - 4.0_RKIND/3.0_RKIND*C_2)*Mc(k,iCell)* &
                (grav*alphaT(k,iCell)*tumd(k,iCell) - grav*betaS(k,iCell)*sumd(k,iCell))
              w2tend5(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 -            &
                alpha2)*(uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz) +         &
                Mc(k,iCell)*(Swumd(k-1,iCell) + Swumd(k,iCell))

              w2tend(k,iCell) = w2tend1(k,iCell) + w2tend2(k,iCell) + &
                w2tend3(k,iCell) + w2tend4(k,iCell) + w2tend5(k,iCell)

              wttend1(k,iCell) = -(Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                wumd(k,iCell)*tumd(k,iCell)
              wttend2(k,iCell) = ((1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell))* &
                wumd(k-1,iCell)*tumd(k-1,iCell)*Mc(k-1,iCell) - (1.0_RKIND - 2.0_RKIND* &
                areaFraction(k+1,iCell))*wumd(k+1,iCell)*tumd(k+1,iCell)*Mc(k+1,iCell)) &
                / (ze(k-1,iCell) - ze(k+1,iCell)) - Mc(k,iCell)*wumd(k,iCell)*Tz
              wttend2(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / (zm(k-1,iCell) - zm(k,iCell)) - &
                  Mc(k,iCell)*wumd(k,iCell)*Tz
              wttend3(k,iCell) = (1.0_RKIND - gamma1)*areaFraction(k,iCell)*  &
                (1.0_RKIND - areaFraction(k,iCell))*grav*(alphaT(k,iCell) &
                *tumd(k,iCell)**2.0 - betaS(k,iCell)*tumd(k,iCell)*sumd(k,iCell))
              wttend4(k,iCell) = -alpha3/4.0*(ut(i1,k,iCell)*Uz + &
                vt(i1,k,iCell)*Vz)
              wttend5(k,iCell) = kappa_FL*(wt(i1,k-1,iCell) -       &
                    wt(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0
              wttend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))* &
                    tumd(k,iCell)*0.5_RKIND*(Swumd(k-1,iCell) + Swumd(k,iCell)) - Mc(k,iCell)* &
                    (1.0_RKIND / areaFraction(k,iCell) * (areaFraction(k-1,iCell)*wt_spsU(k-1,iCell) - &
                    areaFraction(k+1,iCell)*wt_spsU(k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell)) - &
                    1.0_RKIND/(1.0_RKIND - areaFraction(k,iCell))*((1.0_RKIND - areaFraction(k-1,iCell))* &
                    wt_spsD(k-1,iCell) - (1.0_RKIND - areaFraction(k+1,iCell))*wt_spsD(k+1,iCell)) / &
                    (ze(k-1,iCell) - ze(k+1,iCell)))

              wttend(k,iCell) = wttend1(k,iCell) + wttend2(k,iCell) +        &
                wttend3(k,iCell) + wttend4(k,iCell) + wttend5(k,iCell) + wttend6(k,iCell)

              wstend1(k,iCell) = -(Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                wumd(k,iCell)*sumd(k,iCell)
              wstend2(k,iCell) = ((1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell))* &
                wumd(k-1,iCell)*sumd(k-1,iCell)*Mc(k-1,iCell) - (1.0_RKIND - 2.0_RKIND* &
                areaFraction(k+1,iCell))*wumd(k+1,iCell)*sumd(k+1,iCell)*Mc(k+1,iCell)) &
                / (ze(k-1,iCell) - ze(k+1,iCell)) - Mc(k,iCell)*wumd(k,iCell)*Sz
              wstend3(k,iCell) = (1.0_RKIND - gamma1)*grav*areaFraction(k,iCell)* &
                (1.0_RKIND - areaFraction(k,iCell))*(alphaT(k,iCell) &
                *tumd(k,iCell)*sumd(k,iCell) - betaS(k,iCell)*sumd(k,iCell)*sumd(k,iCell))
              wstend4(k,iCell) = -alpha3/4.0*(us(i1,k,iCell)*Uz + &
                vs(i1,k,iCell)*Vz)
              wstend5(k,iCell) = kappa_FL*(ws(i1,k-1,iCell) -       &
                    ws(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0
              wstend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))* &
                    sumd(k,iCell)*0.5_RKIND*(Swumd(k-1,iCell) + Swumd(k,iCell)) - Mc(k,iCell)* &
                    (1.0_RKIND / areaFraction(k,iCell) * (areaFraction(k-1,iCell)*ws_spsU(k-1,iCell) - &
                    areaFraction(k+1,iCell)*ws_spsU(k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell)) - &
                    1.0_RKIND/(1.0_RKIND - areaFraction(k,iCell))*((1.0_RKIND - areaFraction(k-1,iCell))* &
                    ws_spsD(k-1,iCell) - (1.0_RKIND - areaFraction(k+1,iCell))*ws_spsD(k+1,iCell)) / &
                    (ze(k-1,iCell) - ze(k+1,iCell)))

              wstend(k,iCell) = wstend1(k,iCell) + wstend2(k,iCell) + &
                wstend3(k,iCell) + wstend4(k,iCell) + wstend5(k,iCell) + wstend6(k,iCell)

              uwtend1(k,iCell) = -(uw2(k-1,iCell) - uw2(k,iCell)) / dzmid
              uwtend2(k,iCell) = 0.5*((0.8 - 4.0*alpha1/3.0)*0.5*KE**2.0 +  &
                (alpha1 - alpha2)*u2(i1,k,iCell) + (alpha1 +  &
                alpha2 - 2.0_RKIND)*w2(i1,k,iCell))*Uz
              uwtend3(k,iCell) = 0.5_RKIND*(alpha1 - alpha2)*    &
                uv(i1,k,iCell)*Vz
              uwtend4(k,iCell) = beta5*grav*(alphaT(k,iCell)*   &
                ut(i1,k,iCell) - betaS(k,iCell)*us(i1,k,iCell))
              uwtend5(k,iCell) = -tauVel(k,iCell)*uw(i1,k,iCell) + &
                kappa_FL*(uw(i1,k-1,iCell) - uw(i1,k+1,iCell)) / &
                (ze(k-1,iCell) - ze(k+1,iCell))**2.0

              uwtend(k,iCell) = uwtend1(k,iCell) + uwtend2(k,iCell) + &
                uwtend3(k,iCell) + uwtend4(k,iCell) + uwtend5(k,iCell)

              vwtend(k,iCell) = (-(vw2(k-1,iCell) - vw2(k,iCell)) / dzmid +   &
                0.5_RKIND*((0.8_RKIND - 4.0*alpha1/3.0)*0.5_RKIND*KE**2.0_RKIND +  &
                (alpha1 - alpha2)*v2(i1,k,iCell) + (alpha1 +   &
                alpha2 - 2.0_RKIND)*w2(i1,k,iCell))*Vz + 0.5_RKIND*(alpha1 &
                - alpha2)*uv(i1,k,iCell)*Uz + beta5*grav*       &
                (alphaT(k,iCell)*vt(i1,k,iCell) - betaS(k,iCell)*vs(i1,k,iCell))) -            &
                tauVel(k,iCell)*vw(i1,k,iCell) + kappa_FL*(vw(i1,k-1,iCell) &
                - vw(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0

              uvtend(k,iCell) = (-(uvw(k-1,iCell) - uvw(k,iCell)) / dz - &
                (1.0_RKIND - 0.5_RKIND*(alpha1+alpha2))*(uw(i1,k,iCell)*Vz &
                + vw(i1,k,iCell)*Uz)) - tauVel(k,iCell)*uv(i1,k,iCell) +           &
                kappa_VAR*(uv(i1,k-1,iCell) - uv(i1,k+1,iCell)) /       &
                (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

              u2tend1(k,iCell) = -(u2w(k-1,iCell) - u2w(k,iCell)) / dzmid
              u2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - &
                2.0_RKIND)*uw(i1,k,iCell)*Uz
              u2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*vw(i1,k,iCell)*Vz
              u2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1.0_RKIND-beta5)*B
!              u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
!                tauVel(k,iCell)*(KE**2.0_RKIND/3.0_RKIND - u2(i1,k,iCell))
              u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
                tauVel(k,iCell)*(v2(i1,k,iCell) + areaFraction(k,iCell)*( &
                1.0_RKIND - areaFraction(k,iCell))*wumd(k,iCell)**2)/3.0_RKIND

              u2tend(k,iCell) = u2tend1(k,iCell) + u2tend2(k,iCell) + &
                u2tend3(k,iCell) + u2tend4(k,iCell) + u2tend5(k,iCell)

              v2tend1(k,iCell) = -(v2w(k-1,iCell) - v2w(k,iCell)) / dzmid
              v2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - &
                2.0_RKIND)*vw(i1,k,iCell)*Vz
              v2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*uw(i1,k,iCell)*Uz
              v2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1-beta5)*B
!              v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + tauVel(k,iCell)* &
!                (KE**2/3. - v2(i1,k,iCell))
              v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
                tauVel(k,iCell)*(u2(i1,k,iCell) + areaFraction(k,iCell)*( &
                1.0_RKIND - areaFraction(k,iCell))*wumd(k,iCell)**2)/3.0_RKIND

              v2tend(k,iCell) = v2tend1(k,iCell) + v2tend2(k,iCell) + &
                v2tend3(k,iCell) + v2tend4(k,iCell) + v2tend5(k,iCell)

              uttend(k,iCell) = (-(uwt(k-1,iCell) - uwt(k,iCell))/dz -  &
                uw(i1,k,iCell)*Tz - (1.0_RKIND - alpha3)*wt(i1,k,iCell) &
                *Uz) - ut(i1,k,iCell)*tauTemp(k,iCell)

              vttend(k,iCell) = (-(vwt(k-1,iCell) - vwt(k,iCell))/dz -  &
                vw(i1,k,iCell)*Tz - (1.0_RKIND - alpha3)*wt(i1,k,iCell) &
                *Vz) - vt(i1,k,iCell)*tauTemp(k,iCell)

              ustend(k,iCell) = (-(uws(k-1,iCell) - uws(k,iCell))/dz -  &
                uw(i1,k,iCell)*Sz - (1.0_RKIND - alpha3)*ws(i1,k,iCell) &
                *Uz) - us(i1,k,iCell)*tauSalt(k,ICell)

              vstend(k,iCell) = (-(vws(k-1,iCell) - vws(k,iCell))/dz -  &
                vw(i1,k,iCell)*Sz - (1.0_RKIND - alpha3)*ws(i1,k,iCell)* &
                Vz) - vs(i1,k,iCell)*tauSalt(k,iCell)

              t2(i2,k,iCell) = tumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
                (1.0_RKIND-areaFraction(k,iCell))
              s2(i2,k,iCell) = sumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
                (1.0_RKIND-areaFraction(k,iCell))
              ts(i2,k,iCell) = tumd(k,iCell)*sumd(k,iCell)* &
                areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell))

              if(.not. config_adc_use_old_length_scale)then
              !epsilon tendency is next
              KEm1 = 0.5_RKIND*(u2(i1,k-1,iCell) + v2(i1,k-1,iCell) + w2(i1,k-1,iCell))
              KE = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
              KEp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
              tauM1 = 2.0_RKIND*KEm1 / (eps(i1,k-1,iCell) + 1.0E-15_RKIND)
              tau = 2.0_RKIND*KE / (eps(i1,k,iCell) + 1.0E-15_RKIND)
              tauP1 = 2.0_RKIND*KEp1 / (eps(i1,k+1,iCell) + 1.0E-15_RKIND)
              tauAV = 0.5_RKIND*(tauUP + tauDN)

              tomUP = 0.5577/1.3_RKIND*(0.5_RKIND*(KEm1+KE))**2/(0.5_RKIND*(eps(i1,k-1,iCell) + &
                eps(i1,k,iCell)) + 1.0E-15_RKIND)*(eps(i1,k-1,iCell) - &
                eps(i1,k,iCell)) / (ze(k-1,iCell) - ze(k,iCell) + 1.0E-10_RKIND)
              tomDN = 0.5577/1.3_RKIND*(0.5_RKIND*(KE+KEp1))**2/(0.5_RKIND*(eps(i1,k,iCell) + &
                eps(i1,k+1,iCell)) + 1.0E-15_RKIND)*(eps(i1,k,iCell) - &
                eps(i1,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell) + 1.0E-10_RKIND)

              epstend(k,iCell) = (tomUP - tomDN) / (zm(k-1,iCell) - zm(k,iCell) + 1.0E-10_RKIND) - &
              !-0.8_RKIND*((tomUP / (0.5_RKIND*(tauM1 + tau) + 1.0E-10_RKIND) - tomDN / &
               ! (0.5_RKIND*(tau + tauP1) + 1.0E-10_RKIND)) / (zm(k-1,iCell) - zm(k,iCell))) - &
                2.88/(tau + 1.0E-10_RKIND)*uw(i1,k,iCell)*(uvel(k-1,iCell) - uvel(k,iCell)) / &
                (zm(k-1,iCell) - zm(k,iCell)) - &
                2.88/(tau + 1.0E-10_RKIND)*vw(i1,k,iCell)*(vvel(k-1,iCell) - vvel(k,iCell)) / &
                (zm(k-1,iCell) - zm(k,iCell)) + &
                2.88/(tau + 1.0E-10_RKIND)*grav*(alphaT(k,iCell)*wt(i1,k,iCell) - betaS(k,iCell)* &
                ws(i1,k,iCell)) - 3.84*eps(i1,k,iCell)/(tau + 1.0E-10_RKIND)
              endif

              !finally update all subplume fluxes
              if(BVF(k,iCell) <= 0.0_RKIND) then
                lenspsU(k,iCell) = dzmid
                lenspsD(k,iCell) = dzmid
              else
                lenspsU(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsU(i1,k,iCell)/BVF(k,iCell)))
                lenspsD(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsD(i1,k,iCell)/BVF(k,iCell)))
              endif

              KmU(k,iCell) = 0.1_RKIND*lenspsU(k,iCell)*sqrt( KspsU(i1,k,iCell) )
              KhU(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsU(k,iCell)/dz)*KmU(k,iCell)
              wt_spsU(k,iCell) =  -KhU(k,iCell)*Tz
              ws_spsU(k,iCell) =  -KhU(k,iCell)*Sz

              KmD(k,iCell) = 0.1_RKIND*lenspsD(k,iCell)*sqrt( KspsD(i1,k,iCell) )
              KhD(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsD(k,iCell)/dz)*KmD(k,iCell)
              wt_spsD(k,iCell) = -KhD(k,iCell)*Tz
              ws_spsD(k,iCell) = -KhD(k,iCell)*Sz

              !change length scale to lenup and lendown
              Entrainment(k,iCell) = Cww_E*areaFraction(k,iCell)*(1.0_RKIND- &
                areaFraction(k,iCell))*Mc(k,iCell) / ( lenup(k,iCell) + EPSILON )
              Detrainment(k,iCell) = Cww_D*areaFraction(k,iCell)*(1.0_RKIND- &
                areaFraction(k,iCell))*Mc(k,iCell) / ( lendn(k,iCell) + EPSILON )

              if(k==2) then
                Cval = 3.96
              else
                Cval = (0.19_RKIND+0.51_RKIND*lenspsU(k,iCell)/dzmid)
              endif

              KspsUtend(k,iCell) = -grav*(alphaT(k,iCell)*wt_spsU(k,iCell) - &
                betaS(k,iCell)*ws_spsU(k,ICell)) + ((KmU(k-1,iCell) +                 &
                KmU(k,iCell))*(KspsU(i1,k-1,iCell) - KspsU(i1,k,iCell)) /         &
                (ze(k-1,iCell) - ze(k,iCell)) - (KmU(k,iCell) +             &
                KmU(k+1,iCell)) * (KspsU(i1,k,iCell) - KspsU(i1,k+1,iCell)) /     &
                (ze(k,iCell) - ze(k+1,iCell))) / dzmid -                                 &
                Cval*KspsU(i1,k,iCell)**1.5_RKIND/(1.0E-15 + lenspsU(k,iCell)) +                     &
                eps(i1,k,iCell) / (2.0_RKIND*(1.0_RKIND - areaFraction(k,iCell))) + KmU(k,iCell)* &
                (Uz**2.0 + Vz**2.0)

              if(k==2) then
                Cval = 3.96_RKIND
              else
                Cval = (0.19_RKIND+0.51_RKIND*lenspsD(k,iCell)/dz)
              endif

              KspsDtend(k,iCell) = -grav*(alphaT(k,iCell)*wt_spsD(k,iCell) - &
                betaS(k,iCell)*ws_spsD(k,iCell)) + ((KmD(k-1,iCell) +                 &
                KmD(k,iCell))*(KspsD(i1,k-1,iCell) - KspsD(i1,k,iCell)) /         &
                (ze(k-1,iCell) - ze(k,iCell)) - (KmD(k,iCell) +             &
                KmD(k+1,iCell)) * (KspsD(i1,k,iCell) - KspsD(i1,k+1,iCell)) /     &
                (ze(k,iCell) - ze(k+1,iCell))) / dz -                                 &
                Cval*KspsD(i1,k,iCell)**1.5_RKIND /(1.0E-15_RKIND + lenspsD(k,iCell)) +                   &
                eps(i1,k,iCell) / (2.0_RKIND*(areaFraction(k,iCell))) + &
                KmD(k,iCell)*(Uz**2.0 + Vz**2.0)

            enddo ! nVertLevels
          enddo ! nCells

          !In this step we update second moments except w3 which needs updated w2
          do iCell = 1,nCells
            do k=2,nVertLevels
                !update second order moment tendency here
                w2(i2,k,iCell) = (w2(i1,k,iCell) + dt_small*w2tend(k,iCell)) / &
                  (1.0_RKIND + dt_small*tauvVel(k,iCell)*2.0_RKIND/3.0_RKIND)
                if(w2(i2,k,iCell) < epsilon) then
                  w2cliptend(k,iCell) = epsilon - w2(i2,k,iCell)
                  w2(i2,k,iCell) = epsilon
                endif

                if(abs(w2(i2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: w2 out of range, w2 = $r, location = $i, $i", &
                    MPAS_LOG_CRIT,realArgs=(/w2(i2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                u2(i2,k,iCell) = (u2(i1,k,iCell) + dt_small*u2tend(k,iCell)) / &
                  (1.0_RKIND + dt_small*tauVel(k,iCell)*2.0_RKIND/3.0_RKIND)
                if(u2(i2,k,iCell) < epsilon) then
                   u2cliptend(k,iCell) = epsilon - u2(i2,k,iCell)
                   u2(i2,k,iCell) = epsilon
                endif

                v2(i2,k,iCell) = (v2(i1,k,iCell) + dt_small*v2tend(k,iCell)) / &
                  (1.0_RKIND + dt_small*tauVel(k,iCell)*2.0_RKIND/3.0_RKIND)
                if(v2(i2,k,iCell) < epsilon) then
                   v2cliptend(k,iCell) = epsilon - v2(i2,k,iCell)
                   v2(i2,k,iCell) = epsilon
                endif

                uw(i2,k,iCell) = uw(i1,k,iCell) + dt_small*uwtend(k,iCell)
                vw(i2,k,iCell) = vw(i1,k,iCell) + dt_small*vwtend(k,iCell)
                uv(i2,k,iCell) = uv(i1,k,iCell) + dt_small*uvtend(k,iCell)
                ut(i2,k,iCell) = ut(i1,k,iCell) + dt_small*uttend(k,iCell)
                wt(i2,k,iCell) = (wt(i1,k,iCell) + dt_small*wttend(k,iCell)) / &
                  (1.0_RKIND + dt_small*tauTemp(k,iCell))
                vt(i2,k,iCell) = vt(i1,k,iCell) + dt_small*vttend(k,iCell)
                us(i2,k,iCell) = us(i1,k,iCell) + dt_small*ustend(k,iCell)
                vs(i2,k,iCell) = vs(i1,k,iCell) + dt_small*vstend(k,iCell)
                ws(i2,k,iCell) = (ws(i1,k,iCell) + dt_small*wstend(k,iCell)) / &
                   (1.0_RKIND + dt_small*tauSalt(k,iCell))
                if(abs(wt(i2,k,iCell)) > 1.0_RKIND) then
                  call mpas_log_write("ERROR: wt out of range, wt = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/wt(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif

                if(abs(ws(i2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: ws out of range, ws = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/ws(i2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                if(abs(u2(i2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: u2 out of range, u2 = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/u2(i2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                if(abs(v2(i2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: v2 out of range, v2 = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/v2(i2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                if(.not. config_adc_use_old_length_scale)then
                   eps(i2,k,iCell) = eps(i1,k,iCell) + dt_small*epstend(k,iCell)
                   if (eps(i2,k,iCell) < 1.0E-10_RKIND) then
                      eps(i2,k,iCell) = 1.0E-10_RKIND
                   endif
                   !Amrapalli: 16.6 here can be tuned if needed, would show up in u2 and v2 primarily
                   ! also in KEspsU/D
                   KE = 0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
                   length(k,iCell) = KE**1.5 / (c_epsilon * eps(i2,k,iCell))
                 else
                   KE = 0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
                   eps(:,k,iCell) = KE**1.5 / (c_epsilon * length(k,iCell))
                 endif
                KspsU(i2,k,iCell) = max(epsilon,KspsU(i1,k,iCell) + dt_small*KspsUtend(k,iCell))
                KspsD(i2,k,iCell) = max(epsilon,KspsD(i1,k,iCell) + dt_small*KspsDtend(k,iCell))

            enddo !nVertLevels
          enddo !nCells for second order moment tendencies

          ! update the third order w3 and mean fields
          do iCell = 1,nCells
            do k=1,nVertLevels
              w3check = (w2(i2,k,iCell) + w2(i2,k+1,iCell))**1.5
              w3(i2,k,iCell) = min((w3(i1,k,iCell) + dt_small*w3tend(k,iCell)) / &
                (1.0_RKIND + dt_small*tauw3(k,iCell)),w3check)

              if(w2(i2,k,iCell) == 1.0E-8) w3(i2,k,iCell) = 0.0_RKIND

              if(abs(w3(i2,k,iCell)) > 1.0_RKIND) then
                call mpas_log_write("ERROR: w3 out of range, w3 = $r, location = $i, $i", &
                  MPAS_LOG_CRIT,realArgs=(/w3(i2,k,iCell)/),intArgs=(/k,iCell/))
              endif
            enddo
          enddo
      iterCount = iterCount + 1
      call swap_time_levels !probably a better way to do this, with supercycling maybe
!  enddo !end iteration loop -- substepping is done.
    !now that substepping is done, apply computed fluxes to update mean fields.
    !you can collapse this loop too.
    do iCell = 1, nCells
      do k=1,nVertLevels
          utemp = uvel(k,iCell)
          vtemp = vvel(k,iCell)
          uvel(k,iCell) = uvel(k,iCell) - dt_small*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) /  &
                    (ze(k,iCell) - ze(k+1,iCell)) !+ dt_small*fCell(iCell)*vtemp

          vvel(k,iCell) = vvel(k,iCell) - dt_small*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) /  &
                    (ze(k,iCell) - ze(k+1,iCell)) !- dt_small*fCell(iCell)*utemp

          wtSumUp = wt(i1,k,iCell) + areaFraction(k,iCell)*wt_spsD(k,iCell) + &
            (1.0_RKIND - areaFraction(k,iCell))*wt_spsU(k,iCell)
          wtSumDn = wt(i1,k+1,iCell) + areaFraction(k+1,iCell)*wt_spsD(k+1,iCell) + &
            (1.0_RKIND - areaFraction(k+1,iCell))*wt_spsU(k+1,iCell)
          wsSumUp = ws(i1,k,iCell) + areaFraction(k,iCell)*ws_spsD(k,iCell) + &
            (1.0_RKIND - areaFraction(k,iCell))*ws_spsU(k,iCell)
          wsSumDn = ws(i1,k+1,iCell) + areaFraction(k+1,iCell)*ws_spsU(k+1,iCell) + &
            (1.0_RKIND - areaFraction(k+1,iCell))*ws_spsD(k+1,iCell)

          activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt_small*(wtSumUp - &
            wtSumDn) / (ze(k,iCell) - ze(k+1,iCell))
          activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt_small*(wsSumUp - &
            wsSumDn) / (ze(k,iCell) - ze(k+1,iCell))
      enddo
    enddo
  enddo
  end subroutine compute_ADC_tends
end module ocn_adc_mixing_fused
