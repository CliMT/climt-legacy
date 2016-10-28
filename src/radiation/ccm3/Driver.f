c
c This file is used by f2py to generate the Python C wrapper. 
c 
      subroutine driver(
c Input
     $     km,
     $     jm,
     $     im,
     $     idosw,
     $     idolw,
     $     pmidm1,
     $     ps,
     $     tm1,
     $     tg,
     $     qm1, 
     $     o3mmr,
     $     cldf,
     $     clwp,
     $     ciwp,
     $     aldif,
     $     aldir,
     $     asdif,
     $     asdir,
     $     zen,
     $     solin,
     $     r_liq,
     $     r_ice,
     $     co2vmr,
     $     n2ovmr,
     $     ch4vmr,
     $     f11vmr,
     $     f12vmr,
     $     tauvis,
     $     gravit,
     $     cpair,
     $     epsilo,
     $     stebol,
     $     dt,
c Output
     $     tinc,
     $     tdot,
     $     srfflx,
     $     qrs,
     $     qrl,
     $     swflx_out,
     $     lwflx_out,
     $     sw_cf_toa,
     $     sw_cf_srf,
     $     lw_cf_toa,
     $     lw_cf_srf,
     $     sw_toa,
     $     sw_srf,
     $     lw_toa,
     $     lw_srf)

c     Input
      integer km,jm,im
      integer idosw,idolw
      real*8 aldif(im,jm)
      real*8 aldir(im,jm)
      real*8 asdif(im,jm)
      real*8 asdir(im,jm)
      real*8 zen(im,jm)
      real*8 solin(im,jm)
      real*8 ps(im,jm)
      real*8 tg(im,jm)
      real*8 cldf(im,jm,km)
      real*8 clwp(im,jm,km)
      real*8 ciwp(im,jm,km)
      real*8 o3mmr(im,jm,km)
      real*8 r_liq(im,jm,km)
      real*8 r_ice(im,jm,km)
      real*8 pmidm1(im,jm,km)
      real*8 qm1(im,jm,km)
      real*8 tm1(im,jm,km)
      real*8 co2vmr
      real*8 n2ovmr
      real*8 ch4vmr
      real*8 f11vmr
      real*8 f12vmr
      real*8 tauvis
      real*8 gravit
      real*8 cpair
      real*8 epsilo
      real*8 stebol
      real*8 dt
cf2py intent(in,hide)  km,jm,im
cf2py intent(in) aldif, aldir, asdif, asdir, zen, solin, cldf, clwp, ciwp, o3mmr, r_liq, r_ice
cf2py intent(in) pmidm1, ps, qm1, tg, tm1, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr
cf2py intent(in) tauvis, idosw, idolw, gravit, cpair, epsilo, stebol, dt

c    Output
      real*8 tinc(im,jm,km)
      real*8 tdot(im,jm,km)
      real*8 qrs(im,jm,km)
      real*8 qrl(im,jm,km)
      real*8 swflx_out(im,jm,km)
      real*8 lwflx_out(im,jm,km)
      real*8 sw_toa(im,jm)
      real*8 sw_srf(im,jm)
      real*8 lw_toa(im,jm)
      real*8 lw_srf(im,jm)
      real*8 sw_cf_toa(im,jm)
      real*8 sw_cf_srf(im,jm)
      real*8 lw_cf_toa(im,jm)
      real*8 lw_cf_srf(im,jm)
      real*8 srfflx(im,jm)
cf2py intent(out) tinc, tdot, srfflx, qrs, qrl, swflx_out, lwflx_out, 
cf2py intent(out) sw_toa, sw_srf, lw_toa, lw_srf
cf2py intent(out) sw_cf_toa, sw_cf_srf, lw_cf_toa, lw_cf_srf

c Local 
      real*8 swflx(km+1),lwflx(km+1)
      real*8 cldf_col(km), clwp_col(km), ciwp_col(km), o3mmr_col(km)
      real*8 r_liq_col(km), r_ice_col(km), pmidm1_col(km), qm1_col(km)
      real*8 qrs_col(km), qrl_col(km), tm1_col(km)

      print *, 'lat, lon, lev: ' , im, jm, km
      do i=1,im
      do j=1,jm
      do k=1,km
      if (qm1(i,j,k).lt.0.) then
         print*,'qneg!',i,j,k
         qm1(i,j,k)=1.e-9
      endif
      enddo
      enddo
      enddo

      do i=1,im
      do j=1,jm

      cldf_col = cldf(i,j,km:1:-1)
      clwp_col = clwp(i,j,km:1:-1)
      ciwp_col = ciwp(i,j,km:1:-1)
      o3mmr_col = o3mmr(i,j,km:1:-1)
      r_liq_col = r_liq(i,j,km:1:-1)
      r_ice_col = r_ice(i,j,km:1:-1)
      pmidm1_col = pmidm1(i,j,km:1:-1)/100.
      qm1_col = qm1(i,j,km:1:-1)*1000.
      qm1_col(1:4) = 0.0
      qrs_col(:) = 0
      qrl_col(:) = 0
      tm1_col = tm1(i,j,km:1:-1)

      call crm(
     $     aldif(i,j),
     $     aldir(i,j),
     $     asdif(i,j),
     $     asdir(i,j),
     $     zen(i,j),
     $     solin(i,j),
     $     cldf_col,
     $     clwp_col,
     $     ciwp_col,
     $     o3mmr_col,
     $     r_liq_col,
     $     r_ice_col,
     $     pmidm1_col,
     $     ps(i,j)/100.,
     $     qm1_col,
     $     tg(i,j),
     $     tm1_col,
     $     co2vmr,
     $     n2ovmr,
     $     ch4vmr,
     $     f11vmr,
     $     f12vmr,
     $     tauvis,
     $     idosw,
     $     idolw,
     $     gravit,
     $     cpair,
     $     epsilo,
     $     stebol,
     $     qrs_col,
     $     qrl_col,
     $     swflx,
     $     lwflx,
     $     sw_cf_toa(i,j),
     $     sw_cf_srf(i,j),
     $     lw_cf_toa(i,j),
     $     lw_cf_srf(i,j),
     $     sw_toa(i,j),
     $     sw_srf(i,j),
     $     lw_toa(i,j),
     $     lw_srf(i,j) )

      if (idosw .ne. 1) then
         swflx  = solin(i,j)*(1.-asdir(i,j))
         sw_toa = solin(i,j)*(1.-asdir(i,j))
         sw_srf = solin(i,j)*(1.-asdir(i,j))
         sw_cf_toa = 0.
         sw_cf_srf = 0.
         qrs_col(:) = 0.
      endif
      swflx_out(i,j,:) = (swflx(1:km)+swflx(2:km+1))/2.
      swflx_out(i,j,:) = swflx_out(i,j,km:1:-1)
      lwflx_out(i,j,:) = (lwflx(1:km)+lwflx(2:km+1))/2.
      lwflx_out(i,j,:) = lwflx_out(i,j,km:1:-1)
      srfflx(i,j) = swflx(km+1) + lwflx(km+1)
      qrs(i,j,:) = qrs_col(km:1:-1)
      qrl(i,j,:) = qrl_col(km:1:-1)
!      do k=1,km
!        if (isnan(qrs_col(k))) then
!            print *, 'QRS'
!            print *, 'lat, lon, lev: ', i,j,k
!            print *, qrs_col
!            print *, 'Pressure: '
!            print *, pmidm1_col
!            print *, 'Temp: '
!            print *, tm1_col
!            print *, 'Vapour: '
!            print *, qm1_col
!            print *, 'Surf Pressure: '
!            print *, ps(i,j)
!        end if
!        if (isnan(qrl_col(k))) then
!            print *, 'QRL'
!            print *, 'lat, lon, lev: ', i,j,k
!            print *, qrl_col
!            print *, 'Pressure: '
!            print *, pmidm1_col
!            print *, 'Temp: '
!            print *, tm1_col
!            print *, 'Vapour: '
!            print *, qm1_col
!            print *, 'Surf Pressure: '
!            print *, ps(i,j)
!        end if
!
!      enddo
      enddo
      enddo

      tdot   = qrs + qrl
      print *, 'qrs:' 
      print *, qrs(1,1,:)
      print *, 'qrl: '
      print *, qrl(1,1,:)
      print *, 'pmid: '
      print *, pmidm1(1,1,:)
      tinc   = dt*tdot
      qrs = qrs * 86400.
      qrl = qrl * 86400.
      tdot = tdot * 86400.

      end
c-------------------------------------------------------------
      integer function get_nlev()

      integer get_km
      external get_km

      get_nlev = get_km()

      end 
