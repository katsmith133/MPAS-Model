module ocn_adc_mixing_fused

  use mpas_kind_types
  use ocn_constants
  use ocn_turbulence
  use mpas_log
  use ocn_config
  use mpas_timer

  implicit none

  logical :: defineFirst, stopflag

  integer :: i1,i2

  contains

  subroutine swap_time_levels
    call mpas_timer_start('swap_time_levels_fused',.false.)
    i1 = mod(iterCount - 1, 2) + 1
    i2 = mod(iterCount, 2) + 1
    call mpas_timer_stop('swap_time_levels_fused')
  end subroutine swap_time_levels

  subroutine dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(2,nVertLevels,nCells),intent(in) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS
    integer :: i,k, ij

    real,dimension(nVertLevels) :: B, Bup, Bdo
    real,dimension(nVertLevels+1) :: tke, BupEdge, BdoEdge
    real :: sav, tudav, sudav, Tup, Tdo, Sup, Sdo
    real :: s1, z1, zV, sumv, minlen

    !NOTE: will need to convert to some form of displaced density in the mpas framework soon
    !possibly go back to a more traditional length scale formulation
    call mpas_timer_start('dissipation_lengths2_fused',.false.)
    do i=1,nCells
       tke(:) = 0.5*(u2(2,:,i) + v2(2,:,i) + w2(2,:,i))
       do k=1,nVertLevels
          B(k) = -9.806*(-alphaT(k,i)*(activeTracers(1,k,i) - 15.0_RKIND) + betaS(k,i)*   &
                    (activeTracers(2,k,i) - 35.0_RKIND))

          sav = 0.5_RKIND*(areaFraction(k,i) + areaFraction(k+1,i))
          tudav = 0.5_RKIND*(tumd(k,i) + tumd(k+1,i))
          sudav = 0.5*(sumd(k,i) + sumd(k+1,i))

          Tup = activeTracers(1,k,i) + (1.0_RKIND - sav)*tudav
          Tdo = activeTracers(1,k,i) - sav*tudav
          Sup = activeTracers(2,k,i) + (1.0_RKIND - sav)*sudav
          Sdo = activeTracers(2,k,i) - sav*sudav

          Bup(k) = -9.806*(-alphaT(k,i)*(Tup - 15.0) + betaS(k,i)*(Sup - 35.0))
          Bdo(k) = -9.806*(-alphaT(k,i)*(Tdo - 15.0) + betaS(k,i)*(Sdo - 35.0))

          if(k>1) THEN
             BupEdge(k) = 0.5*(Bup(k-1) + Bup(k))
             BdoEdge(k) = 0.5*(Bdo(k-1) + Bdo(k))
          endif
       enddo

       BdoEdge(nVertLevels+1) = BdoEdge(nVertLevels)
       BupEdge(nVertLevels+1) = BupEdge(nVertLevels)

       BdoEdge(1) = BdoEdge(2)
       BupEdge(1) = BupEdge(2)

       do k=2,nVertLevels

          sumv = 0
          ij=k
          lenup(k,i) = 0.0_RKIND
          do while(sumv <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv + (BupEdge(k) - Bup(ij))*(ze(ij-1,i)-ze(ij,i))
             lenup(k,i) =  lenup(k,i) + abs(ze(ij-1,i)-ze(ij,i))
             ij = ij + 1

             if(sumv > tke(k)) THEN
                ij = ij - 1
                s1 = sumv
                z1 = ze(ij,i)
                zV = ze(ij-1,i)

                sumv = sumv - (BupEdge(k) - Bup(ij))*(ze(ij-1,i)-ze(ij,i))
                lenup(k,i) = lenup(k,i) - abs(ze(ij,i)-ze(ij+1,i))
                lenup(k,i) = max(0.55,lenup(k,i) + abs((z1-zV)/(s1 - sumv)*(tke(k)-sumv)))
            !    if(Bup(k-1) - Bup(k) < 0) then
            !            minlen = abs(zm(k-1,i) - zm(k,i))
            !    else
            !            minlen = 0.5_RKIND
            !    endif
            !    if(k==5) then
            !      print *, 'len pre = ',lenup(k,i)
            !    endif
!
!                lenup(k,i) = max(minlen, lenup(k,i) + tke(k)/(1.0E-15 + sumv - (BupEdge(k) -         &
!                                Bup(ij))))
!                              if(k==5) then
!                                print *, lenup(k,i),minlen, tke(k)
!                              endif
                              exit
             endif

         end do

        !find lendown
        sumv = 0
        ij=k
        lendn(k,i) = 0.0_RKIND
        do while(sumv <= tke(k) .and. ij>1)
           sumv = sumv - (BdoEdge(k) - Bdo(ij-1))*(ze(ij-1,i)-ze(ij,i))
           lendn(k,i) = lendn(k,i) + abs(ze(ij-1,i)-ze(ij,i))
           ij = ij - 1

           if(sumv > tke(k)) THEN
              ij = ij + 1
              s1 = sumv
              z1 = ze(ij,i)
              zV = ze(ij-1,i)
              sumv = sumv + (BdoEdge(k) - Bdo(ij-1))*(ze(ij-1,i)-ze(ij,i))
              lendn(k,i) = lendn(k,i) - abs(ze(ij-1,i)-ze(ij,i))
              lendn(k,i) = max(0.55,lendn(k,i) + ((zV-z1)/(s1 - sumv)*(tke(k) - sumv)))

           exit
           endif
        enddo

        length(k,i) = (2.0_RKIND*lenup(k,i)*lendn(k,i))  &
          / (lenup(k,i) + lendn(k,i))
      enddo
   enddo

   length(1,:) = 0.55_RKIND
   length(nVertLevels+1,:) = 0.55_RKIND
   call mpas_timer_stop('dissipation_lengths2_fused')   

  end subroutine dissipation_lengths2

  subroutine build_diagnostic_arrays(nCells,nVertLevels,BVF,wtsfc,wssfc, &
          uwsfc, vwsfc, alphaT,betaS, boundaryLayerDepth)
    !construct dTdz, dSdz, dbdz
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, uwsfc, vwsfc
    real,dimension(nCells),intent(in) :: boundaryLayerDepth
    real,dimension(nVertLevels,nCells),intent(in) :: BVF, alphaT, betaS
    integer :: iCell, k, idx, idx2
    
    real,dimension(nCells) :: wstar
    logical :: first
    real :: coeff, maximum, Q, KE

    call mpas_timer_start('build_diagnostic_arrays_fused',.false.)
    first = .true.

    do iCell=1,nCells
      Q = grav*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))* &
        boundaryLayerDepth(iCell)
      if(Q > 0) then
        wstar(iCell) = abs(Q)**(1.0/3.0)
      else
        wstar(iCell) = 0.0
      endif

      u2(:,1,iCell) = 4.0*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      v2(:,1,iCell) =4.*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      uw(:,1,iCell) = -uwsfc(iCell)
      vw(:,1,iCell) = vwsfc(iCell)
      wt(:,1,iCell) = wtsfc(iCell)
      ws(:,1,iCell) = wssfc(iCell)
      KE = 0.5*(u2(2,1,iCell) + v2(2,1,iCell))
      eps(:,1,iCell) = KE**1.5/(0.5*(ze(1,iCell) - ze(2,iCell) + 1.0E-10_RKIND))
    enddo
    call mpas_timer_stop('build_diagnostic_arrays_fused')

  end subroutine build_diagnostic_arrays

  subroutine build_areaFraction_properties(nCells,nVertLevels,wtsfc, wssfc, alphaT,  &
      betaS, boundaryLayerDepth, tlev)

      integer,intent(in) :: nCells, nVertLevels, tlev
      real,dimension(nCells),intent(in) :: wtsfc, wssfc, boundaryLayerDepth
      real, dimension(nVertLevels,nCells) :: alphaT, betaS

      real :: wtav, McAv, sigav, tumdav, wumdav, sumdav, wb, bld, wstar
      real :: Sw, w3av, lareaFraction, wsav
      integer :: iCell,k

      call mpas_timer_start('build_areaFraction_properties_fused',.false.)
      do iCell = 1,nCells
        tumd(1,iCell) = 0.0_RKIND
        wumd(1,iCell) = 0.0_RKIND
        areaFraction(1,iCell) = 0.5_RKIND
        Mc(1,iCell) = 0.0_RKIND
        wb = grav*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))

        if(wb > 0.0) then
          wb = grav*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))
          wstar = (abs(0.4_RKIND*boundaryLayerDepth(iCell)*wb))**(1.0_RKIND/3.0_RKIND)
          w2t(1,iCell) = -0.3_RKIND*wstar * wtsfc(iCell)
          !Below FIXME!
          w2s(1,iCell) = 0.3_RKIND*wstar * wssfc(iCell)
        else
          w2t(1,iCell) = 0.0_RKIND
          w2s(1,iCell) = 0.0_RKIND
        endif

        !try new boundary condition derived from PDF
        sigav = 0.5_RKIND*(areaFraction(1,iCell) + areaFraction(2,iCell))
        wtav = 0.5_RKIND*(wt(tlev,1,iCell) + wt(tlev,2,iCell))
        McAv = 0.5_RKIND*(w2(tlev,1,iCell) + w2(tlev,2,iCell))
        w2t(1,iCell) = (1.0_RKIND - 2.0_RKIND*sigav)*wtav*sqrt(McAv) / (EPSILON + &
          sigav*(1.0_RKIND-sigav))

        do k=2,nVertLevels
          w3av = 0.5_RKIND*(w3(tlev,k-1,iCell) + w3(tlev,k,iCell))

          Sw = w3av / (max(w2(tlev,k,iCell)**1.5_RKIND,epsilon))
          lareaFraction = 0.5_RKIND - 0.5_RKIND*Sw / sqrt(4.0_RKIND + Sw**2.0_RKIND)

          if(lareaFraction < 0.01_RKIND) lareaFraction = 0.01_RKIND
          if(lareaFraction > 0.99_RKIND) lareaFraction = 0.99_RKIND

          areaFraction(k,iCell) = lareaFraction
          wumd(k,iCell) = sqrt(w2(tlev,k,iCell) / (areaFraction(k,iCell) * &
             (1.0_RKIND - areaFraction(k,iCell))))
          Mc(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - &
             areaFraction(k,iCell)) * wumd(k,iCell)

          tumd(k,iCell) = wt(tlev,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))
          sumd(k,iCell) = ws(tlev,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))
        enddo

        areaFraction(nVertLevels+1, iCell) = areaFraction(nVertLevels, iCell)
        do k=2,nVertLevels
          sigav = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
          tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
          sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
          wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
          w2t(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*tumdav
          w2s(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*sumdav
        enddo
      enddo !iCell loop end
      call mpas_timer_stop('build_areaFraction_properties_fused')
      
  end subroutine build_areaFraction_properties

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
      real :: B, Cval, diff, wtav, dzmid, Ksps, Sz, Tz
      real :: w3av, tempMoment

      real,dimension(nVertLevels,nCells) :: taupt, taups, taupv
      real,dimension(nVertLevels+1,nCells) :: wttend, wstend, uwtend, vwtend, u2tend, v2tend, &
        uttend, vttend, ustend, vstend, uvtend
      real,dimension(nVertLevels) :: epstend
      real,dimension(nVertLevels+1,nCells) :: KspsUtend, KspsDtend
    
!      call mpas_timer_start('compute_ADC_tends_fused',.false.)
      dt_small = config_adc_timestep
      niter = dt / dt_small

!      call swap_time_levels
! at start of iterations make sure arrays are synchronized. MAY NOT BE NEEDED

      do iCell=1,nCells
        do k=1,nVertLevels
          KspsU(2,k,iCell) = KspsU(1,k,iCell)
          KspsD(2,k,iCell) = KspsD(1,k,iCell)
          eps(2,k,iCell) = eps(1,k,iCell)
          u2(2,k,iCell) = u2(1,k,iCell)
          v2(2,k,iCell) = v2(1,k,iCell)
          w2(2,k,iCell) = w2(1,k,iCell)
          ut(2,k,iCell) = ut(1,k,iCell)
          vt(2,k,iCell) = vt(1,k,iCell)
          wt(2,k,iCell) = wt(1,k,iCell)
          us(2,k,iCell) = us(1,k,iCell)
          vs(2,k,iCell) = vs(1,k,iCell)
          ws(2,k,iCell) = ws(1,k,iCell)
          w3(2,k,iCell) = w3(1,k,iCell)
          uw(2,k,iCell) = uw(1,k,iCell)
          vw(2,k,iCell) = vw(1,k,iCell)
          t2(2,k,iCell) = t2(1,k,iCell)
          s2(2,k,iCell) = s2(1,k,iCell)
          uv(2,k,iCell) = uv(1,k,iCell)
          ts(2,k,iCell) = ts(1,k,iCell)
        enddo
      enddo

      call mpas_timer_start('compute_ADC_tends_fused',.false.)
      do iIter=1,niter

          !First do all the center based quantities (third order moments)
          !Call these few things now but inline later
          call build_diagnostic_arrays(nCells,nVertLevels,BVF,wtsfc,wssfc,  &
            uwsfc,vwsfc,alphaT,betaS,boundaryLayerDepth)

          call build_areaFraction_properties(nCells, nVertLevels,wtsfc, wssfc, &
                                    alphaT, betaS, boundaryLayerDepth, 2)

          !note to ashwath -- as you do the porting you need an update device here for the previous
          !two function calls since those calculations will be on host.
          !Now do TOM update
          do iCell = 1,nCells
            do k=1,nVertLevels
              Eav = 0.5*(Entrainment(k+1,iCell) + Entrainment(k,iCell))
              Dav = 0.5*(Detrainment(k+1,iCell) + Detrainment(k,iCell))
              u2av = 0.5*(u2(2,k,iCell) + u2(2,k+1,iCell))
              v2av = 0.5*(v2(2,k,iCell) + v2(2,k+1,iCell))
              w2av = 0.5*(w2(2,k,iCell) + w2(2,k+1,iCell))

              sigav = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
              wumdav = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
              tumdav = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
              sumdav = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
              KspsUav = 0.5*(KspsU(2,k,iCell) + KspsU(2,k+1,iCell))
              KspsDav = 0.5*(KspsD(2,k,iCell) + KspsD(2,k+1,iCell))
              Mcav = 0.5*(Mc(k,iCell) + Mc(k+1,iCell))
              lenav = 0.5*(length(k,iCell) + length(k+1,iCell))
              if(k==nVertLevels) then
                sigavp1 = 0.5*(areaFraction(k,iCell))
                wumdAvp1 = 0.5*(wumd(k,iCell))
                tumdAvp1 = 0.5*(tumd(k,iCell))
                sumdAvp1 = 0.5*(sumd(k,iCell))
              else
                sigavp1 = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
                wumdAvp1 = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
                tumdAvp1 = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
                sumdAvp1 = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
              endif

              KEsps = sigav*KspsUav+ (1.0 - sigav)*KspsDav
              KE = sqrt((u2av+v2av+w2av) + 0.0*KEsps)

              dz = ze(k,iCell) - ze(k+1,iCell)
              !KE = sqrt(areaFraction(k,iCell)*KspsUav + (1.0 - areaFraction(k,iCell))*KspsDav)
              Swup = - 2.0/3.0*(KspsU(2,k,iCell) - KspsU(2,k+1,iCell)) / dz &
                 - 2.0/3.0*KspsUav*(log(areaFraction(k,iCell)) -           &
                log(areaFraction(k+1,iCell))) / dz +         &
                2.0/3.0*(KspsD(2,k,iCell) - KspsD(2,k+1,iCell)) / dz        &
               + 2.0/3.0*KspsDav*(log(1.0-areaFraction(k,iCell)) -         &
                log(1.0-areaFraction(k+1,iCell))) /  dz

              w3tend1(k,iCell) = wumdav**3.0*(Eav*(3.0*sigav - 2.0) + Dav*(3.0*sigav - 1.0))
              w3tend2(k,iCell) =  wumdav**3.0*(6.0*sigav**2.0 - 6.0*sigav + 1)*             &
                    (areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))*  &
                    wumd(k,iCell) - areaFraction(k+1,iCell)*(1.0_RKIND - &
                    areaFraction(k+1,iCell))*wumd(k+1,iCell))/dz
              w3tend3(k,iCell) = -1.5_RKIND*sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)* &
                    wumdav**2.0_RKIND*((1.0_RKIND - 2.0_RKIND*areaFraction(k,iCell))*      &
                    wumd(k,iCell)**2.0_RKIND - (1.0_RKIND -                                &
                    2.0_RKIND*areaFraction(k+1,iCell))*wumd(k+1,iCell)**2.0_RKIND) / dz
              w3tend4(k,iCell) = 3.0_RKIND*(1.0_RKIND - 2.0_RKIND*sigav)*Mcav*wumdav*Swup -          &
                    C_mom_w3*KE / (1E-15_RKIND+ sqrt(2.0_RKIND)*lenAv)*w3(2,k,iCell)
              w3tend5(k,iCell) = 3.0_RKIND*grav*(alphaT(k,iCell)*w2t(k,iCell) - &
                    betaS(k,iCell)*w2s(k,iCell))*0.9_RKIND

              w3tend(k,iCell) = w3tend1(k,iCell) + w3tend2(k,iCell) + w3tend3(k,iCell) + &
                w3tend4(k,iCell) + w3tend5(k,iCell)

              if(k>1 .and. k < nVertLevels) then
                 w3tend(k,iCell) = w3tend(k,iCell) + kappa_w3*(w3(2,k-1,iCell) &
                    - w3(2,k+1,iCell)) / (zm(k-1,iCell) - zm(k+1,iCell))**2.0
              endif

            ! now get all the downgradient TOMs
              Ksps = 0.5_RKIND*((areaFraction(k,iCell)*KspsU(2,k,iCell) + (1.0_RKIND - &
                areaFraction(k,iCell))*KspsD(2,k,iCell)) + (areaFraction(k+1,iCell)* &
                KspsU(2,k+1,iCell) + (1.0_RKIND-areaFraction(k+1,iCell))* &
                KspsD(2,k+1,iCell)))
              KE = 0.5_RKIND*(u2(2,k,iCell) + v2(2,k,iCell) + w2(2,k,iCell))
              KEp1 = 0.5_RKIND*(u2(2,k+1,iCell) + v2(2,k+1,iCell) + w2(2,k+1,iCell))
              lenav = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))
              diff = C_mom * sqrt(0.5_RKIND*(KE + KEp1)) / (1.0E-10_RKIND + lenav)
              dz = ze(k,iCell) - ze(k+1,iCell)
              uw2(k,iCell) = -diff*(uw(2,k,iCell) - uw(2,k+1,iCell)) / dz
              vw2(k,iCell) = -diff*(vw(2,k,iCell) - vw(2,k+1,iCell)) / dz
              u2w(k,iCell) = -diff*(u2(2,k,iCell) - u2(2,k+1,iCell)) / dz
              v2w(k,iCell) = -diff*(v2(2,k,iCell) - v2(2,k+1,iCell)) / dz
              uvw(k,iCell) = -diff*(uv(2,k,iCell) - uv(2,k+1,iCell)) / dz

              diff = C_therm*sqrt(0.5*(KE + KEp1)) / (1.0E-10_RKIND + lenav)
              uwt(k,iCell) = -diff*(ut(2,k,iCell) - ut(2,k+1,iCell)) / dz
              vwt(k,iCell) = -diff*(vt(2,k,iCell) - vt(2,k+1,iCell)) / dz
              uws(k,iCell) = -diff*(us(2,k,iCell) - us(2,k+1,iCell)) / dz
              vws(k,iCell) = -diff*(vs(2,k,iCell) - vs(2,k+1,iCell)) / dz

            enddo !nVertLevels
          enddo !nCells

          if(.not. config_adc_use_old_length_scale)then
            call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
          endif

          !next all second order moment tendencies also apply tendencies
          do iCell=1,nCells
            epstend(1) = 0.0_RKIND
            u2cliptend(:,iCell) = 0.0_RKIND
            v2cliptend(:,iCell) = 0.0_RKIND

            do k=2,nVertLevels
              dzmid = (zm(k-1,iCell) - zm(k,iCell))
              dz = ze(k,iCell) - ze(k+1,iCell)
              B = grav*(alphaT(k,iCell)*wt(2,k,iCell) - betaS(k,iCell)*ws(2,k,iCell))

              sigav = 0.5*(areaFraction(k,iCell) + areaFraction(k-1,iCell))
              wumdav = 0.5*(wumd(k,iCell) + wumd(k-1,iCell))
              tumdav = 0.5*(tumd(k,iCell) + tumd(k-1,iCell))
              sumdav = 0.5*(sumd(k,iCell) + sumd(k-1,iCell))
              KspsUav = 0.5*(KspsU(2,k,iCell) + KspsU(2,k-1,iCell))
              KspsDav = 0.5*(KspsD(2,k,iCell) + KspsD(2,k-1,iCell))
              Mcav = 0.5*(Mc(k,iCell) + Mc(k-1,iCell))

              sigavp1 = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
              KspsUavp1 = 0.5*(KspsU(2,k,iCell) + KspsU(2,k+1,iCell))
              KspsDavp1 = 0.5*(KspsD(2,k,iCell) + KspsD(2,k+1,iCell))

              Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dzmid
              Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dzmid
              Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dzmid
              Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dzmid

              KEsps = areaFraction(k,iCell)*KspsU(2,k,iCell) +  &
                        (1.0_RKIND-areaFraction(k,iCell))*KspsD(2,k,iCell)
              KE = sqrt((u2(2,k,iCell) + v2(2,k,iCell) + w2(2,k,iCell)) + 0.0*KEsps)
              Swup = grav*alphaT(k,iCell)*tumd(k,iCell) - grav*        &
                betaS(k,iCell)*sumd(k,iCell) - 2.0_RKIND/3.0_RKIND*(1.0_RKIND/        &
                areaFraction(k,iCell)*(sigAv*KspsUav - sigavp1*KspsUavp1) / dzmid -     &
                1.0_RKIND/(1.0_RKIND - areaFraction(k,iCell))*((1.0_RKIND - sigav)* &
                KspsDav - (1.0 - sigavp1)*KspsDavp1) / dzmid)

              w2tend1(k,iCell) = -wumd(k,iCell)**2.0_RKIND*(         &
                Entrainment(k,iCell) + Detrainment(k,iCell))
              w2tend2(k,iCell) = -(w3(2,k-1,iCell) - w3(2,k,iCell)) / dzmid
              w2tend3(k,iCell) = -C_1*KE / (1.0E-15_RKIND + sqrt(2.0_RKIND)* &
                length(k,iCell))*(w2(2,k,iCell) - KE**2.0_RKIND/3.0_RKIND)
              w2tend4(k,iCell) = 2.0_RKIND*Mc(k,iCell)*Swup + 4.0_RKIND / &
                3.0_RKIND*C_2*areaFraction(k,iCell)*(1.0_RKIND -          &
                areaFraction(k,iCell))*wumd(k,iCell)*(grav*     &
                alphaT(k,iCell)*tumd(k,iCell) - grav*betaS(k,iCell)*      &
                sumd(k,iCell))
              w2tend5(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 -            &
                alpha2)*(uw(2,k,iCell)*Uz + vw(2,k,iCell)*Vz)

              w2tend(k,iCell) = w2tend1(k,iCell) + w2tend2(k,iCell) + &
                w2tend3(k,iCell) + w2tend4(k,iCell) + w2tend5(k,iCell)

              taupt(k,iCell) = KE / (sqrt(2.0)*c_pt*length(k,iCell))
              taups(k,iCell) = KE / (sqrt(2.0)*c_ps*length(k,iCell))
              taupv(k,iCell) = KE / (c_pv*length(k,iCell))

              wttend1(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / dzmid

              wttend2(k,iCell) = -w2(2,k,iCell)*Tz
              wttend3(k,iCell) = (1.0_RKIND - gamma1)*grav*(alphaT(k,iCell) &
                *t2(2,k,iCell) - betaS(k,iCell)*ts(2,k,iCell))
              wttend4(k,iCell) = -alpha3/4.0*(ut(2,k,iCell)*Uz + &
                vt(2,k,iCell)*Vz)
              wttend5(k,iCell) = kappa_FL*(wt(2,k-1,iCell) -       &
                    wt(2,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0

              wttend(k,iCell) = wttend1(k,iCell) + wttend2(k,iCell) +        &
                wttend3(k,iCell) + wttend4(k,iCell) + wttend5(k,iCell)

              wstend1(k,iCell) = -(w2s(k-1,iCell) - w2s(k,iCell)) / dzmid
              wstend2(k,iCell) = -w2(2,k,iCell)*Sz
              wstend3(k,iCell) = (1.0_RKIND - gamma1)*grav*(alphaT(k,iCell)* &
                ts(2,k,iCell) - betaS(k,iCell)*s2(2,k,iCell))
              wstend4(k,iCell) = -alpha3/4.0*(us(2,k,iCell)*Uz + &
                vs(2,k,iCell)*Vz)
              wsTend5(k,iCell) = kappa_FL*(ws(2,k-1,iCell) -     &
                ws(2,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

              wstend(k,iCell) = wstend1(k,iCell) + wstend2(k,iCell) + &
                wstend3(k,iCell) + wstend4(k,iCell) + wstend5(k,iCell)

              uwtend1(k,iCell) = -(uw2(k-1,iCell) - uw2(k,iCell)) / dzmid
              uwtend2(k,iCell) = 0.5*((0.8-4.0/3.0*alpha1)*0.5*KE**2.0 +  &
                (alpha1 - alpha2)*u2(2,k,iCell) + (alpha1 +  &
                alpha2 - 2.0_RKIND)*w2(2,k,iCell))*Uz
              uwtend3(k,iCell) = 0.5_RKIND*(alpha1 - alpha2)*    &
                uv(2,k,iCell)*Vz
              uwtend4(k,iCell) = beta5*grav*(alphaT(k,iCell)*   &
                ut(2,k,iCell) - betaS(k,iCell)*us(2,k,iCell))
              uwtend5(k,iCell) = - 2.0_RKIND*taupv(k,iCell)*uw(2,k,iCell) + &
                kappa_FL*(uw(2,k-1,iCell) - uw(2,k+1,iCell)) / &
                (ze(k-1,iCell) - ze(k+1,iCell))**2.0

              uwtend(k,iCell) = uwtend1(k,iCell) + uwtend2(k,iCell) + &
                uwtend3(k,iCell) + uwtend4(k,iCell) + uwtend5(k,iCell)

              vwtend(k,iCell) = (-(vw2(k-1,iCell) - vw2(k,iCell)) / dzmid +   &
                0.5_RKIND*((0.8_RKIND-4.0_RKIND/3.0_RKIND*alpha1)*0.5_RKIND*KE**2.0_RKIND +  &
                (alpha1 - alpha2)*v2(2,k,iCell) + (alpha1 -   &
                alpha2 - 2.0_RKIND)*w2(2,k,iCell))*Vz + 0.5_RKIND*(alpha1 &
                - alpha2)*uv(2,k,iCell)*Uz + beta5*grav*       &
                (alphaT(k,iCell)*vt(2,k,iCell) - betaS(k,iCell)*vs(2,k,iCell))) -            &
                taupv(k,iCell)*vw(2,k,iCell) + kappa_FL*(vw(2,k-1,iCell) &
                - vw(2,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0

              uvtend(k,iCell) = (-(uvw(k-1,iCell) - uvw(k,iCell)) / dz + &
                (1.0_RKIND - 0.5_RKIND*(alpha1+alpha2))*(uw(2,k,iCell)*Vz &
                + vw(2,k,iCell)*Uz)) - taupv(k,iCell)*uv(2,k,iCell) +           &
                kappa_VAR*(uv(2,k-1,iCell) - uv(2,k+1,iCell)) /       &
                (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

              u2tend1(k,iCell) = -(u2w(k-1,iCell) - u2w(k,iCell)) / dzmid
              u2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - &
                2.0_RKIND)*uw(2,k,iCell)*Uz
              u2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*vw(2,k,iCell)*Vz
              u2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1.0_RKIND-beta5)*B
              u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(2,k,iCell) + &
                taupv(k,iCell)*(KE**2.0_RKIND/3.0_RKIND - u2(2,k,iCell))

              u2tend(k,iCell) = u2tend1(k,iCell) + u2tend2(k,iCell) + &
                u2tend3(k,iCell) + u2tend4(k,iCell) + u2tend5(k,iCell)

              v2tend1(k,iCell) = -(v2w(k-1,iCell) - v2w(k,iCell)) / dzmid
              v2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 +alpha2 - &
                2.0_RKIND)*vw(2,k,iCell)*Vz
              v2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*uw(2,k,iCell)*Uz
              v2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1-beta5)*B
              v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(2,k,iCell) + taupv(k,iCell)* &
                (KE**2/3. - v2(2,k,iCell))

              v2tend(k,iCell) = v2tend1(k,iCell) + v2tend2(k,iCell) + &
                v2tend3(k,iCell) + v2tend4(k,iCell) + v2tend5(k,iCell)

              uttend(k,iCell) = (-(uwt(k-1,iCell) - uwt(k,iCell))/dz -  &
                uw(2,k,iCell)*Tz - (1.0_RKIND - alpha3)*wt(2,k,iCell) &
                *Uz) - ut(2,k,iCell)*taupt(k,iCell)

              vttend(k,iCell) = (-(vwt(k-1,iCell) - vwt(k,iCell))/dz -  &
                vw(2,k,iCell)*Tz - (1.0_RKIND - alpha3)*wt(2,k,iCell) &
                *Vz) - vt(2,k,iCell)*taupt(k,iCell)

              ustend(k,iCell) = (-(uws(k-1,iCell) - uws(k,iCell))/dz -  &
                uw(2,k,iCell)*Sz - (1.0_RKIND - alpha3)*ws(2,k,iCell) &
                *Uz) - us(2,k,iCell)*taups(k,ICell)

              vstend(k,iCell) = (-(vws(k-1,iCell) - vws(k,iCell))/dz -  &
                vw(2,k,iCell)*Sz - (1.0_RKIND - alpha3)*ws(2,k,iCell)* &
                Vz) - vs(2,k,iCell)*taups(k,iCell)

              t2(2,k,iCell) = tumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
                (1.0_RKIND-areaFraction(k,iCell))
              s2(2,k,iCell) = sumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
                (1.0_RKIND-areaFraction(k,iCell))
              ts(2,k,iCell) = tumd(k,iCell)*sumd(k,iCell)* &
                areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell))

              if(.not. config_adc_use_old_length_scale)then
              !epsilon tendency is next
              KEm1 = 0.5_RKIND*(u2(2,k-1,iCell) + v2(2,k-1,iCell) + w2(2,k-1,iCell))
              KE = 0.5_RKIND*(u2(2,k,iCell) + v2(2,k,iCell) + w2(2,k,iCell))
              KEp1 = 0.5_RKIND*(u2(2,k+1,iCell) + v2(2,k+1,iCell) + w2(2,k+1,iCell))
              tauM1 = 2.0_RKIND*KEm1 / (eps(2,k-1,iCell) + 1.0E-15_RKIND)
              tau = 2.0_RKIND*KE / (eps(2,k,iCell) + 1.0E-15_RKIND)
              tauP1 = 2.0_RKIND*KEp1 / (eps(2,k+1,iCell) + 1.0E-15_RKIND)
              tauAV = 0.5_RKIND*(tauUP + tauDN)

              tomUP = 0.5577/1.3_RKIND*(0.5_RKIND*(KEm1+KE))**2/(0.5_RKIND*(eps(2,k-1,iCell) + &
                eps(2,k,iCell)) + 1.0E-15_RKIND)*(eps(2,k-1,iCell) - &
                eps(2,k,iCell)) / (ze(k-1,iCell) - ze(k,iCell) + 1.0E-10_RKIND)
              tomDN = 0.5577/1.3_RKIND*(0.5_RKIND*(KE+KEp1))**2/(0.5_RKIND*(eps(2,k,iCell) + &
                eps(2,k+1,iCell)) + 1.0E-15_RKIND)*(eps(2,k,iCell) - &
                eps(2,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell) + 1.0E-10_RKIND)

              epstend(k) = (tomUP - tomDN) / (zm(k-1,iCell) - zm(k,iCell) + 1.0E-10_RKIND) - &
              !-0.8_RKIND*((tomUP / (0.5_RKIND*(tauM1 + tau) + 1.0E-10_RKIND) - tomDN / &
               ! (0.5_RKIND*(tau + tauP1) + 1.0E-10_RKIND)) / (zm(k-1,iCell) - zm(k,iCell))) - &
                2.88/(tau + 1.0E-10_RKIND)*uw(2,k,iCell)*(uvel(k-1,iCell) - uvel(k,iCell)) / &
                (zm(k-1,iCell) - zm(k,iCell)) - &
                2.88/(tau + 1.0E-10_RKIND)*vw(2,k,iCell)*(vvel(k-1,iCell) - vvel(k,iCell)) / &
                (zm(k-1,iCell) - zm(k,iCell)) + &
                2.88/(tau + 1.0E-10_RKIND)*grav*(alphaT(k,iCell)*wt(2,k,iCell) - betaS(k,iCell)* &
                ws(2,k,iCell)) - 3.84*eps(2,k,iCell)/(tau + 1.0E-10_RKIND)
              endif

              !finally update all subplume fluxes
              if(BVF(k,iCell) <= 0.0_RKIND) then
                lenspsU(k,iCell) = dzmid
                lenspsD(k,iCell) = dzmid
              else
                lenspsU(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsU(2,k,iCell)/BVF(k,iCell)))
                lenspsD(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsD(2,k,iCell)/BVF(k,iCell)))
              endif

              KmU(k,iCell) = 0.1_RKIND*lenspsU(k,iCell)*sqrt( KspsU(2,k,iCell) )
              KhU(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsU(k,iCell)/dz)*KmU(k,iCell)
              wt_spsU(k,iCell) =  -KhU(k,iCell)*Tz
              ws_spsU(k,iCell) =  -KhU(k,iCell)*Sz

              KmD(k,iCell) = 0.1_RKIND*lenspsD(k,iCell)*sqrt( KspsD(2,k,iCell) )
              KhD(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsD(k,iCell)/dz)*KmD(k,iCell)
              wt_spsD(k,iCell) = -KhD(k,iCell)*Tz
              ws_spsD(k,iCell) = -KhD(k,iCell)*Sz

              Entrainment(k,iCell) = Cww_E*areaFraction(k,iCell)*(1.0_RKIND- &
                areaFraction(k,iCell))*Mc(k,iCell) / ( lenspsD(k,iCell) + EPSILON )
              Detrainment(k,iCell) = Cww_D*areaFraction(k,iCell)*(1.0_RKIND- &
                areaFraction(k,iCell))*Mc(k,iCell) / ( lenspsU(k,iCell) + EPSILON )

              if(k==2) then
                Cval = 3.96
              else
                Cval = (0.19_RKIND+0.51_RKIND*lenspsU(k,iCell)/dzmid)
              endif

              KspsUtend(k,iCell) = grav*(alphaT(k,iCell)*wt_spsU(k,iCell) - &
                betaS(k,iCell)*ws_spsU(k,ICell)) + ((KmU(k-1,iCell) +                 &
                KmU(k,iCell))*(KspsU(2,k-1,iCell) - KspsU(2,k,iCell)) /         &
                (ze(k-1,iCell) - ze(k,iCell)) - (KmU(k,iCell) +             &
                KmU(k+1,iCell)) * (KspsU(2,k,iCell) - KspsU(2,k+1,iCell)) /     &
                (ze(k,iCell) - ze(k+1,iCell))) / dzmid -                                 &
                Cval*KspsU(2,k,iCell)**1.5_RKIND/(1.0E-15 + lenspsU(k,iCell)) +                     &
                eps(2,k,iCell) / (2.0_RKIND*areaFraction(k,iCell))

              if(k==2) then
                Cval = 3.96_RKIND
              else
                Cval = (0.19_RKIND+0.51_RKIND*lenspsD(k,iCell)/dz)
              endif

              KspsDtend(k,iCell) = grav*(alphaT(k,iCell)*wt_spsD(k,iCell) - &
                betaS(k,iCell)*ws_spsD(k,iCell)) + ((KmD(k-1,iCell) +                 &
                KmD(k,iCell))*(KspsD(2,k-1,iCell) - KspsD(2,k,iCell)) /         &
                (ze(k-1,iCell) - ze(k,iCell)) - (KmD(k,iCell) +             &
                KmD(k+1,iCell)) * (KspsD(2,k,iCell) - KspsD(2,k+1,iCell)) /     &
                (ze(k,iCell) - ze(k+1,iCell))) / dz -                                 &
                Cval*KspsD(2,k,iCell)**1.5_RKIND / lenspsD(k,iCell) +                   &
                eps(2,k,iCell) / (2.0_RKIND*(1.0_RKIND - areaFraction(k,iCell)))

                !update second order moment tendency here
                w2(2,k,iCell) = w2(1,k,iCell) + dt_small*w2tend(k,iCell)
                if(w2(2,k,iCell) < epsilon) then
                  w2cliptend(k,iCell) = epsilon-w2(2,k,iCell)
                  w2(2,k,iCell) = epsilon
                endif

                if(abs(w2(2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: w2 out of range, w2 = $r, location = $i, $i", &
                    MPAS_LOG_CRIT,realArgs=(/w2(2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                u2(2,k,iCell) = u2(1,k,iCell) + dt_small*u2tend(k,iCell)
                if(u2(2,k,iCell) < 0.0_RKIND) then
                   u2cliptend(k,iCell) = -u2(2,k,iCell)
                   u2(2,k,iCell) = 0.0_RKIND
                endif

                v2(2,k,iCell) = v2(1,k,iCell) + dt_small*v2tend(k,iCell)
                if(v2(2,k,iCell) < 0.0_RKIND) then
                   v2cliptend(k,iCell) = -v2(2,k,iCell)
                   v2(2,k,iCell) = 0.0_RKIND
                endif

                uw(2,k,iCell) = uw(1,k,iCell) + dt_small*uwtend(k,iCell)
                vw(2,k,iCell) = vw(1,k,iCell) + dt_small*vwtend(k,iCell)
                uv(2,k,iCell) = uv(1,k,iCell) + dt_small*uvtend(k,iCell)
                ut(2,k,iCell) = ut(1,k,iCell) + dt_small*uttend(k,iCell)
                wt(2,k,iCell) = (wt(1,k,iCell) + dt_small*wttend(k,iCell)) / &
                  (1.0_RKIND + dt_small*taupt(k,iCell))
                vt(2,k,iCell) = vt(1,k,iCell) + dt_small*vttend(k,iCell)
                us(2,k,iCell) = us(1,k,iCell) + dt_small*ustend(k,iCell)
                vs(2,k,iCell) = vs(1,k,iCell) + dt_small*vstend(k,iCell)
                ws(2,k,iCell) = (ws(1,k,iCell) + dt_small*wstend(k,iCell)) / &
                   (1.0_RKIND + dt_small*taups(k,iCell))
                if(abs(wt(2,k,iCell)) > 1.0_RKIND) then
                  call mpas_log_write("ERROR: wt out of range, wt = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/wt(2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                if(abs(ws(2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: ws out of range, ws = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/ws(2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                if(abs(u2(2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: u2 out of range, u2 = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/u2(2,k,iCell)/),intArgs=(/k,iCell/))
                endif

                if(abs(v2(2,k,iCell)) > 1.0_RKIND) then
                   call mpas_log_write("ERROR: v2 out of range, v2 = $r, location = $i, $i", &
                     MPAS_LOG_CRIT,realArgs=(/v2(2,k,iCell)/),intArgs=(/k,iCell/))
                endif

								if(.not. config_adc_use_old_length_scale)then
								   call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
                   eps(2,k,iCell) = eps(1,k,iCell) + dt_small*epstend(k)
                   if (eps(2,k,iCell) < 1.0E-10_RKIND) then
                      eps(2,k,iCell) = 1.0E-10_RKIND
                   endif
                   KE = 0.5_RKIND*(u2(2,k,iCell) + v2(2,k,iCell) + w2(2,k,iCell))
                   length(k,iCell) = KE**1.5 / eps(2,k,iCell)
								endif
                KspsU(2,k,iCell) = max(epsilon,KspsU(1,k,iCell) + dt_small*KspsUtend(k,iCell))
                KspsD(2,k,iCell) = max(epsilon,KspsD(1,k,iCell) + dt_small*KspsDtend(k,iCell))

            enddo !nVertLevels
            lenspsU(nVertLevels,iCell) = lenspsU(nVertLevels-1,iCell)
            lenspsD(nVertLevels,iCell) = lenspsD(nVertLevels-1,iCell)

            length(nVertLevels+1,iCell) = 1e-15_RKIND
            length(1,iCell) = length(2,iCell) !1e-1
          enddo !nCells for second order moment tendencies

          ! update the third order w3 and mean fields
          do iCell = 1,nCells
            do k=1,nVertLevels
              w3check = (w2(2,k,iCell) + w2(2,k+1,iCell))**1.5
              w3(2,k,iCell) = min(w3(1,k,iCell) + dt_small*w3tend(k,iCell),w3check)

              if(abs(w3(2,k,iCell)) > 1.0_RKIND) then
                call mpas_log_write("ERROR: w3 out of range, w3 = $r, location = $i, $i", &
                  MPAS_LOG_CRIT,realArgs=(/w3(2,k,iCell)/),intArgs=(/k,iCell/))
              endif
            enddo
          enddo
      iterCount = iterCount + 1
      !copy 2 into 1 next for the next sub iteration.
call mpas_timer_start('compute_ADC_tends_copy_nextLoop',.false.)
      do iCell=1,nCells
        do k=1,nVertLevels
          KspsU(1,k,iCell) = KspsU(2,k,iCell)
          KspsD(1,k,iCell) = KspsD(2,k,iCell)
          eps(1,k,iCell) = eps(2,k,iCell)
          u2(1,k,iCell) = u2(2,k,iCell)
          v2(1,k,iCell) = v2(2,k,iCell)
          w2(1,k,iCell) = w2(2,k,iCell)
          ut(1,k,iCell) = ut(2,k,iCell)
          vt(1,k,iCell) = vt(2,k,iCell)
          wt(1,k,iCell) = wt(2,k,iCell)
          us(1,k,iCell) = us(2,k,iCell)
          vs(1,k,iCell) = vs(2,k,iCell)
          ws(1,k,iCell) = ws(2,k,iCell)
          w3(1,k,iCell) = w3(2,k,iCell)
          uw(1,k,iCell) = uw(2,k,iCell)
          vw(1,k,iCell) = vw(2,k,iCell)
          t2(1,k,iCell) = t2(2,k,iCell)
          s2(1,k,iCell) = s2(2,k,iCell)
          uv(1,k,iCell) = uv(2,k,iCell)
          ts(1,k,iCell) = ts(2,k,iCell)
        enddo
      enddo
      call mpas_timer_stop('compute_ADC_tends_copy_nextLoop')
!      call swap_time_levels !probably a better way to do this, with supercycling maybe
    enddo !end iteration loop -- substepping is done.

    !now that substepping is done, apply computed fluxes to update mean fields.
    do iCell = 1, nCells
      do k=1,nVertLevels
          utemp = uvel(k,iCell)
          vtemp = vvel(k,iCell)
          uvel(k,iCell) = uvel(k,iCell) - dt*(uw(2,k,iCell) - uw(2,k+1,iCell)) /  &
                    (ze(k,iCell) - ze(k+1,iCell)) !+ dt_small*fCell(iCell)*vtemp

          vvel(k,iCell) = vvel(k,iCell) - dt*(vw(2,k,iCell) - vw(2,k+1,iCell)) /  &
                    (ze(k,iCell) - ze(k+1,iCell)) !- dt_small*fCell(iCell)*utemp

          activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt*(wt(2,k,iCell) - &
            wt(2,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell))
          activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt*(ws(2,k,iCell) - &
            ws(2,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell))
      enddo
    enddo

    call mpas_timer_stop('compute_ADC_tends_fused')

  end subroutine compute_ADC_tends
end module ocn_adc_mixing_fused
