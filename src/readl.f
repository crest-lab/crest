!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

C     *****************************************************************         
                                                            
      SUBROUTINE READL(A1,X,N)                              
      IMPLICIT REAL*8 (A-H,O-Z)                             
      CHARACTER*(*) A1                                      
      DIMENSION X(*)                                        
      I=0                                                   
      IS=1                                                  
  10  I=I+1                                                 
      X(I)=READAA(A1,IS,IB,IE)                              
      IF(IB.GT.0 .AND. IE.GT.0) THEN                        
                                IS=IE                       
                                GOTO 10                     
      ENDIF                                                 
      N=I-1                                                 
      RETURN                                                
      END                                                   
                                                            
                                                            
      FUNCTION READAA(A,ISTART,IEND,IEND2)                  
      IMPLICIT REAL*8 (A-H,O-Z)                             
      REAL*8 READAA                                         
      CHARACTER*(*) A                                       
      NINE=ICHAR('9')                                       
      IZERO=ICHAR('0')                                      
      MINUS=ICHAR('-')                                      
      IDOT=ICHAR('.')                                       
      ND=ICHAR('D')                                         
      NE=ICHAR('E')                                         
      IBL=ICHAR(' ')                                        
      IEND=0                                                
      IEND2=0                                               
      IDIG=0                                                
      C1=0                                                  
      C2=0                                                  
      ONE=1.D0                                              
      X = 1.D0                                              
      NL=LEN(A) 
      DO 10 J=ISTART,NL-1                                   
         N=ICHAR(A(J:J))                                    
         M=ICHAR(A(J+1:J+1)) 
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20  
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO        
     1 .OR. M.EQ.IDOT)) GOTO 20                             
   10 CONTINUE                                              
      READAA=0.D0                                           
      RETURN                                                
   20 CONTINUE                                              
      IEND=J                                                
      DO 30 I=J,NL                                          
         N=ICHAR(A(I:I))                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                  
            IDIG=IDIG+1                                     
            IF (IDIG.GT.10) GOTO 60                         
            C1=C1*10+N-IZERO                                
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                 
            ONE=-1.D0                                       
         ELSEIF(N.EQ.IDOT) THEN                             
            GOTO 40                                         
         ELSE                                               
            GOTO 60                                         
         ENDIF                                              
   30 CONTINUE                                              
   40 CONTINUE                                              
      IDIG=0                                                
      DO 50 II=I+1,NL                                       
         N=ICHAR(A(II:II))                                  
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                  
            IDIG=IDIG+1                                     
            IF (IDIG.GT.10) GOTO 60                         
            C2=C2*10+N-IZERO                                
            X = X /10                                       
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                
            X=-X                                            
         ELSE                                               
            GOTO 60                                         
         ENDIF                                              
   50 CONTINUE                                              
C                                                           
C PUT THE PIECES TOGETHER                                   
C                                                           
   60 CONTINUE                                              
      READAA= ONE * ( C1 + C2 * X)                          
      DO 55 J=IEND,NL                                       
         N=ICHAR(A(J:J))                                    
         IEND2=J                                            
         IF(N.EQ.IBL)RETURN                                 
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                       
      RETURN                                                
                                                            
   57 C1=0.0D0                                              
      ONE=1.0D0                                             
      DO 31 I=J+1,NL                                        
         N=ICHAR(A(I:I))                                    
         IEND2=I                                            
         IF(N.EQ.IBL)GOTO 70                                
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO  
         IF(N.EQ.MINUS)ONE=-1.0D0                           
   31 CONTINUE                                              
   61 CONTINUE                                              
   70 READAA=READAA*10**(ONE*C1)                            
      RETURN                                                
      END                                                   
