SUBROUTINE CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
    IMPLICIT NONE
    
    INTEGER:: INFO
    REAL*8::Pi,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF
    
    IF (INFO.EQ.1)THEN
        IF(TETA1.LT.TETA2)THEN
            TETAi=TETA2
            TETAF=TETA1+2.D0*PI
        ELSE
            TETAi=TETA2
            TETAF=TETA1
        ENDIF
    ELSE IF(INFO.EQ.2) THEN
        IF(TETA1.LT.TETA2)THEN
            TETAi=TETA1
            TETAF=TETA2
        ELSE
            TETAi=TETA1
            TETAF=TETA2+2.D0*Pi
        ENDIF                 
    ELSE IF (INFO.EQ.3)THEN
        IF(TETA1.LT.TETA3)THEN
            TETAi=TETA1
            TETAF=TETA3
        ELSE
            TETAi=TETA1
            TETAF=TETA3+2.D0*Pi
        ENDIF   
    ELSE IF(INFO.EQ.31)THEN
        IF(TETA3.LT.TETA2)THEN
            TETAi=TETA3
            TETAF=TETA2
        ELSE
            TETAi=TETA3
            TETAF=TETA2+2.D0*Pi
        ENDIF   
    ELSE IF(INFO.EQ.4)THEN
        IF(TETA4.LT.TETA3)THEN
            TETAi=TETA3
            TETAF=TETA4+2.D0*PI
        ELSE
            TETAi=TETA3
            TETAF=TETA4
        ENDIF                   
    ELSE IF (INFO.EQ.41)THEN
        IF(TETA1.LT.TETA3)THEN
            TETAi=TETA1
            TETAF=TETA3
        ELSE
            TETAi=TETA1
            TETAF=TETA3+2.D0*Pi
        ENDIF                   
    ELSE IF (INFO.EQ.42)THEN
        IF(TETA4.LT.TETA2)THEN
            TETAi=TETA4
            TETAF=TETA2
        ELSE
            TETAi=TETA4
            TETAF=TETA2+2.D0*Pi
        ENDIF                   
    ENDIF
END SUBROUTINE CORRIGE_TETA