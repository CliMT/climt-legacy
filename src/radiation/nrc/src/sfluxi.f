      SUBROUTINE SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,RSFI,WNOI,DWNI,
     *                  COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,
     *                  fluxupi,fluxdni,FZEROI,NGWI)

C  GCM2.0  Feb 2003

      implicit none

      include "grid.h"
      include "radinc.h"

      real*8 planckir(L_NSPECTI,8501)
      common / plnr / planckir

      integer NLEVRAD, L, NW, NG, NTS, NTT

      real*8 TLEV(L_LEVELS), PLEV(L_LEVELS)
      real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8 FMNETI(L_NLEVRAD)
      real*8 WNOI(L_NSPECTI), DWNI(L_NSPECTI)
      real*8 DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 FMUPI(L_NLEVRAD), FMDI(L_NLEVRAD)
      real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 GWEIGHT(L_NGAUSS), NFLUXTOPI
      real*8 FTOPUP

      real*8 UBARI, RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
      real*8 PLANCK, PLTOP
      real*8 fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD)
      real*8 FZEROI(L_NSPECTI)
      real*8 fzero, tlimit
      integer NGWI(L_NSPECTI)
      common /tlim1 / tlimit

C======================================================================C
 
      NLEVRAD = L_NLEVRAD
 
C     ZERO THE NET FLUXES
    
      NFLUXTOPI = 0.0
 
      DO L=1,L_NLAYRAD
        FMNETI(L)  = 0.0
        FLUXUPI(L) = 0.0
        FLUXDNI(L) = 0.0
      END DO
 
C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      TTOP  = TLEV(2)
      TSURF = TLEV(L_LEVELS)

      NTS   = TSURF*10.0D0-499
      NTT   = TTOP *10.0D0-499

      DO 501 NW=1,L_NSPECTI

C       SURFACE EMISSIONS - INDEPENDENT OF GAUSS POINTS

        BSURF = (1.-RSFI)*PLANCKIR(NW,NTS)
        PLTOP = PLANCKIR(NW,NTT)

C  If FZEROI(NW) = 1, then the k-coefficients are zero - skip to the
C  special Gauss point at the end.
 
        FZERO = FZEROI(NW)
        IF(FZERO.ge.0.99) goto 40

        DO NG=1,NGWI(NW)-1
          fzero = fzero + (1.0-FZEROI(NW))*GWEIGHT(NG)
        END DO
 
        DO NG=NGWI(NW),L_NGAUSS-1
         
C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
C         CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
C         OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY
 
          TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
          BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP
 
C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
          
          CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),
     *                TAUCUMI(1,NW,NG),
     *                WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,
     *                BSURF,FTOPUP,FMUPI,FMDI)
 
C         NOW CALCULATE THE CUMULATIVE IR NET FLUX

          NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*GWEIGHT(NG)*
     *                           (1.0-FZEROI(NW))

          DO L=1,L_NLEVRAD-1

C           CORRECT FOR THE WAVENUMBER INTERVALS

            FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*
     *                              GWEIGHT(NG)*(1.0-FZEROI(NW))
            FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*
     *                                (1.0-FZEROI(NW))
            FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)*
     *                                (1.0-FZEROI(NW))
          END DO

       END DO       !End NGAUSS LOOP

   40  CONTINUE

C      SPECIAL 17th Gauss point

       NG     = L_NGAUSS

       TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
       BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP

C      WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C      CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C      WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 

       CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),
     *                TAUCUMI(1,NW,NG),
     *                WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,
     *                BSURF,FTOPUP,FMUPI,FMDI)
 
C      NOW CALCULATE THE CUMULATIVE IR NET FLUX

       NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*FZERO

       DO L=1,L_NLEVRAD-1

C        CORRECT FOR THE WAVENUMBER INTERVALS

         FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*FZERO
         FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*FZERO
         FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*FZERO
       END DO
  
  501 CONTINUE      !End Spectral Interval LOOP

C *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****

      RETURN
      END
