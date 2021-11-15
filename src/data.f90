!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme
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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Previously rcov.f
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine setrcov(rcov)
      real*8 :: rcov(94)

      rcov(  1 )= 0.80628308d0 
      rcov(  2 )= 1.15903197d0
      rcov(  3 )= 3.02356173d0
      rcov(  4 )= 2.36845659d0
      rcov(  5 )= 1.94011865d0
      rcov(  6 )= 1.88972601d0
      rcov(  7 )= 1.78894056d0
      rcov(  8 )= 1.58736983d0
      rcov(  9 )= 1.61256616d0
      rcov( 10 )= 1.68815527d0
      rcov( 11 )= 3.52748848d0
      rcov( 12 )= 3.14954334d0
      rcov( 13 )= 2.84718717d0
      rcov( 14 )= 2.62041997d0
      rcov( 15 )= 2.77159820d0
      rcov( 16 )= 2.57002732d0
      rcov( 17 )= 2.49443835d0
      rcov( 18 )= 2.41884923d0
      rcov( 19 )= 4.43455700d0
      rcov( 20 )= 3.88023730d0
      rcov( 21 )= 3.35111422d0
      rcov( 22 )= 3.07395437d0
      rcov( 23 )= 3.04875805d0
      rcov( 24 )= 2.77159820d0
      rcov( 25 )= 2.69600923d0
      rcov( 26 )= 2.62041997d0
      rcov( 27 )= 2.51963467d0
      rcov( 28 )= 2.49443835d0
      rcov( 29 )= 2.54483100d0
      rcov( 30 )= 2.74640188d0
      rcov( 31 )= 2.82199085d0
      rcov( 32 )= 2.74640188d0
      rcov( 33 )= 2.89757982d0
      rcov( 34 )= 2.77159820d0
      rcov( 35 )= 2.87238349d0
      rcov( 36 )= 2.94797246d0
      rcov( 37 )= 4.76210950d0
      rcov( 38 )= 4.20778980d0
      rcov( 39 )= 3.70386304d0
      rcov( 40 )= 3.50229216d0
      rcov( 41 )= 3.32591790d0
      rcov( 42 )= 3.12434702d0
      rcov( 43 )= 2.89757982d0
      rcov( 44 )= 2.84718717d0
      rcov( 45 )= 2.84718717d0
      rcov( 46 )= 2.72120556d0
      rcov( 47 )= 2.89757982d0
      rcov( 48 )= 3.09915070d0
      rcov( 49 )= 3.22513231d0
      rcov( 50 )= 3.17473967d0
      rcov( 51 )= 3.17473967d0
      rcov( 52 )= 3.09915070d0
      rcov( 53 )= 3.32591790d0
      rcov( 54 )= 3.30072128d0
      rcov( 55 )= 5.26603625d0
      rcov( 56 )= 4.43455700d0
      rcov( 57 )= 4.08180818d0
      rcov( 58 )= 3.70386304d0
      rcov( 59 )= 3.98102289d0
      rcov( 60 )= 3.95582657d0
      rcov( 61 )= 3.93062995d0
      rcov( 62 )= 3.90543362d0
      rcov( 63 )= 3.80464833d0
      rcov( 64 )= 3.82984466d0
      rcov( 65 )= 3.80464833d0
      rcov( 66 )= 3.77945201d0
      rcov( 67 )= 3.75425569d0
      rcov( 68 )= 3.75425569d0
      rcov( 69 )= 3.72905937d0
      rcov( 70 )= 3.85504098d0
      rcov( 71 )= 3.67866672d0
      rcov( 72 )= 3.45189952d0
      rcov( 73 )= 3.30072128d0
      rcov( 74 )= 3.09915070d0
      rcov( 75 )= 2.97316878d0
      rcov( 76 )= 2.92277614d0
      rcov( 77 )= 2.79679452d0
      rcov( 78 )= 2.82199085d0
      rcov( 79 )= 2.84718717d0
      rcov( 80 )= 3.32591790d0
      rcov( 81 )= 3.27552496d0
      rcov( 82 )= 3.27552496d0
      rcov( 83 )= 3.42670319d0
      rcov( 84 )= 3.30072128d0
      rcov( 85 )= 3.47709584d0
      rcov( 86 )= 3.57788113d0
      rcov( 87 )= 5.06446567d0
      rcov( 88 )= 4.56053862d0
      rcov( 89 )= 4.20778980d0
      rcov( 90 )= 3.98102289d0
      rcov( 91 )= 3.82984466d0
      rcov( 92 )= 3.85504098d0
      rcov( 93 )= 3.88023730d0
      rcov( 94 )= 3.90543362d0

      return 

end subroutine setrcov
 
subroutine setr2r4(r2r4)
      real*8 :: r2r4(94)

      r2r4( 1 )= 2.00734898
      r2r4( 2 )= 1.56637132
      r2r4( 3 )= 5.01986934
      r2r4( 4 )= 3.85379032
      r2r4( 5 )= 3.64446594
      r2r4( 6 )= 3.10492822
      r2r4( 7 )= 2.71175247
      r2r4( 8 )= 2.59361680
      r2r4( 9 )= 2.38825250
      r2r4( 10 )= 2.21522516
      r2r4( 11 )= 6.58585536
      r2r4( 12 )= 5.46295967
      r2r4( 13 )= 5.65216669
      r2r4( 14 )= 4.88284902
      r2r4( 15 )= 4.29727576
      r2r4( 16 )= 4.04108902
      r2r4( 17 )= 3.72932356
      r2r4( 18 )= 3.44677275
      r2r4( 19 )= 7.97762753
      r2r4( 20 )= 7.07623947
      r2r4( 21 )= 6.60844053
      r2r4( 22 )= 6.28791364
      r2r4( 23 )= 6.07728703
      r2r4( 24 )= 5.54643096
      r2r4( 25 )= 5.80491167
      r2r4( 26 )= 5.58415602
      r2r4( 27 )= 5.41374528
      r2r4( 28 )= 5.28497229
      r2r4( 29 )= 5.22592821
      r2r4( 30 )= 5.09817141
      r2r4( 31 )= 6.12149689
      r2r4( 32 )= 5.54083734
      r2r4( 33 )= 5.06696878
      r2r4( 34 )= 4.87005108
      r2r4( 35 )= 4.59089647
      r2r4( 36 )= 4.31176304
      r2r4( 37 )= 9.55461698
      r2r4( 38 )= 8.67396077
      r2r4( 39 )= 7.97210197
      r2r4( 40 )= 7.43439917
      r2r4( 41 )= 6.58711862
      r2r4( 42 )= 6.19536215
      r2r4( 43 )= 6.01517290
      r2r4( 44 )= 5.81623410
      r2r4( 45 )= 5.65710424
      r2r4( 46 )= 5.52640661
      r2r4( 47 )= 5.44263305
      r2r4( 48 )= 5.58285373
      r2r4( 49 )= 7.02081898
      r2r4( 50 )= 6.46815523
      r2r4( 51 )= 5.98089120
      r2r4( 52 )= 5.81686657
      r2r4( 53 )= 5.53321815
      r2r4( 54 )= 5.25477007
      r2r4( 55 )= 11.02204549
      r2r4( 56 )= 10.15679528
      r2r4( 57 )= 9.35167836
      r2r4( 58 )= 9.06926079
      r2r4( 59 )= 8.97241155
      r2r4( 60 )= 8.90092807
      r2r4( 61 )= 8.85984840
      r2r4( 62 )= 8.81736827
      r2r4( 63 )= 8.79317710
      r2r4( 64 )= 7.89969626
      r2r4( 65 )= 8.80588454
      r2r4( 66 )= 8.42439218
      r2r4( 67 )= 8.54289262
      r2r4( 68 )= 8.47583370
      r2r4( 69 )= 8.45090888
      r2r4( 70 )= 8.47339339
      r2r4( 71 )= 7.83525634
      r2r4( 72 )= 8.20702843
      r2r4( 73 )= 7.70559063
      r2r4( 74 )= 7.32755997
      r2r4( 75 )= 7.03887381
      r2r4( 76 )= 6.68978720
      r2r4( 77 )= 6.05450052
      r2r4( 78 )= 5.88752022
      r2r4( 79 )= 5.70661499
      r2r4( 80 )= 5.78450695
      r2r4( 81 )= 7.79780729
      r2r4( 82 )= 7.26443867
      r2r4( 83 )= 6.78151984
      r2r4( 84 )= 6.67883169
      r2r4( 85 )= 6.39024318
      r2r4( 86 )= 6.09527958
      r2r4( 87 )= 11.79156076
      r2r4( 88 )= 11.10997644
      r2r4( 89 )= 9.51377795
      r2r4( 90 )= 8.67197068
      r2r4( 91 )= 8.77140725
      r2r4( 92 )= 8.65402716
      r2r4( 93 )= 8.53923501
      r2r4( 94 )= 8.85024712

      return 

end subroutine setr2r4


