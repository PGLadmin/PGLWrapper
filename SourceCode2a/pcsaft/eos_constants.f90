!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module EOS_CONSTANTS
!> \brief  Module contains EOS-constants
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module eos_constants

  use PARAMETERS, only: dp
  implicit none
  save

  private
  public :: load_eos_constants, load_dd_constants, load_qq_constants, load_dq_constants

  real(dp), public, dimension( 0:6, 3 )             :: ap, bp

  real(dp), public, allocatable, dimension(:,:,:)   :: qqp2, qqp4, ddp2, ddp4, dqp2, dqp4
  real(dp), public, allocatable, dimension(:,:,:,:) :: qqp3, ddp3, dqp3

contains


!> \file eos_const.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine eos_const
!
!> \brief provides the constants of the PC-SAFT EOS
!
! This subroutine provides the constants of the PC-SAFT EOS.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_eos_constants !( ap, bp )

  use pcsaft_pure_and_binary_parameters, only: dipole, qudpole, dipole_quad
  !-----------------------------------------------------------------------------

  call load_pc_saft_constants

  if ( dipole ) call load_dd_constants !( ddp2, ddp3, ddp4 )

  if ( qudpole ) call load_qq_constants !( qqp2, qqp3, qqp4 )

  if ( dipole_quad ) call load_dq_constants !( dqp2, dqp3, dqp4 )
  
end subroutine load_eos_constants



!> \file eos_const.f90
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine eos_const
!
!> \brief provides the constants of the PC-SAFT EOS
!
! This subroutine provides the constants of the PC-SAFT EOS.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_pc_saft_constants

  !-----------------------------------------------------------------------------
  ! constants of the dispersion term
  !-----------------------------------------------------------------------------

  ap(0,1) =  0.91056314451539_dp
  ap(0,2) = -0.30840169182720_dp
  ap(0,3) = -0.09061483509767_dp
  ap(1,1) =  0.63612814494991_dp
  ap(1,2) =  0.18605311591713_dp
  ap(1,3) =  0.45278428063920_dp
  ap(2,1) =  2.68613478913903_dp
  ap(2,2) = -2.50300472586548_dp
  ap(2,3) =  0.59627007280101_dp
  ap(3,1) = -26.5473624914884_dp
  ap(3,2) =  21.4197936296668_dp
  ap(3,3) = -1.72418291311787_dp
  ap(4,1) =  97.7592087835073_dp
  ap(4,2) = -65.2558853303492_dp
  ap(4,3) = -4.13021125311661_dp
  ap(5,1) = -159.591540865600_dp
  ap(5,2) =  83.3186804808856_dp
  ap(5,3) =  13.7766318697211_dp
  ap(6,1) =  91.2977740839123_dp
  ap(6,2) = -33.7469229297323_dp
  ap(6,3) = -8.67284703679646_dp

  bp(0,1) =  0.72409469413165_dp
  bp(0,2) = -0.57554980753450_dp
  bp(0,3) =  0.09768831158356_dp
  bp(1,1) =  2.23827918609380_dp
  bp(1,2) =  0.69950955214436_dp
  bp(1,3) = -0.25575749816100_dp
  bp(2,1) = -4.00258494846342_dp
  bp(2,2) =  3.89256733895307_dp
  bp(2,3) = -9.15585615297321_dp
  bp(3,1) = -21.0035768148465_dp
  bp(3,2) = -17.2154716477721_dp
  bp(3,3) =  20.6420759743972_dp
  bp(4,1) =  26.8556413626615_dp
  bp(4,2) =  192.672264465249_dp
  bp(4,3) = -38.8044300520628_dp
  bp(5,1) =  206.551338406619_dp
  bp(5,2) = -161.826461648765_dp
  bp(5,3) =  93.6267740770146_dp
  bp(6,1) = -355.602356122079_dp
  bp(6,2) = -165.207693455561_dp
  bp(6,3) = -29.6669055851473_dp

end subroutine load_pc_saft_constants



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine dq_const
!
!> \brief provides the constants of the dipole-quadrupole term
!
! This subr. provides the constants of the dipole-quadrupole term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_dq_constants

  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mseg

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k
  real(dp)                                   :: mdq_i, mdq_j, mdq_k
  real(dp)                                   :: mf1, mf2, msegij
  !-----------------------------------------------------------------------------

  do i = 1, ncomp

     mdq_i = min( mseg(i), 2.0_dp )

     do j = 1, ncomp

        mdq_j = min( mseg(j), 2.0_dp )

        msegij = ( mdq_i * mdq_j )**0.5
        mf1 = ( msegij - 1.0_dp ) / msegij
        mf2 = mf1 * ( msegij - 2.0_dp ) / msegij

        dqp2(i,j,0) =  0.697094963_dp + mf1*(-0.673459279_dp) + mf2*0.670340770_dp
        dqp2(i,j,1) = -0.633554144_dp + mf1*(-1.425899106_dp) + mf2*(-4.338471826_dp)
        dqp2(i,j,2) =  2.945509028_dp + mf1 * 4.19441392_dp   + mf2*7.234168360_dp
        dqp2(i,j,3) = -1.467027314_dp + mf1 * 1.0266216_dp
        dqp2(i,j,4) = 0.0_dp

        dqp4(i,j,0) = -0.484038322_dp + mf1 * 0.67651011_dp   + mf2*(-1.167560146_dp)
        dqp4(i,j,1) =  1.970405465_dp + mf1*(-3.013867512_dp) + mf2*2.13488432_dp
        dqp4(i,j,2) = -2.118572671_dp + mf1 * 0.46742656_dp
        dqp4(i,j,3) = 0.0_dp
        dqp4(i,j,4) = 0.0_dp


        do k = 1, ncomp
           mdq_k = min( mseg(k), 2.0_dp )
           msegij = ( mdq_i * mdq_j * mdq_k )**(1.0_dp/3.0_dp)
           mf1 = ( msegij - 1.0_dp ) / msegij
           mf2 = ( msegij - 2.0_dp ) / msegij
           dqp3(i,j,k,0) = 0.795009692_dp + mf1*(-2.099579397_dp)
           dqp3(i,j,k,1) = 3.386863396_dp + mf1*(-5.941376392_dp)
           dqp3(i,j,k,2) = 0.475106328_dp + mf1*(-0.178820384_dp)
           dqp3(i,j,k,3) = 0.0_dp
           dqp3(i,j,k,4) = 0.0_dp
        end do

     end do

  end do

end subroutine load_dq_constants


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine dd_const
!
!> \brief provides the constants of the dipole-term
! This subroutine provides the constants of the dipole-term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_dd_constants

  use PARAMETERS, only: PI
  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mseg

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k
  real(dp)                                   :: mdd_i, mdd_j, mdd_k
  real(dp)                                   :: mf1, mf2, msegij, sin2t
  !-----------------------------------------------------------------------------

  sin2t = SIN( 0.0_dp * PI / 180.0_dp )
  sin2t = sin2t*sin2t

  do i = 1, ncomp

     mdd_i = min( mseg(i), 2.0_dp )

     do j = 1, ncomp

        mdd_j = min( mseg(j), 2.0_dp )

        msegij = ( mdd_i * mdd_j )**0.5_dp
        mf1 = ( msegij - 1.0_dp ) / msegij
        mf2 = mf1 * ( msegij - 2.0_dp ) / msegij

        ddp2(i,j,0) =  0.30435038064_dp + mf1*(0.95346405973_dp+0.201436_dp*sin2t)  &
                      + mf2*(-1.16100802773_dp-1.74114_dp*sin2t)
        ddp2(i,j,1) = -0.13585877707_dp + mf1*(-1.83963831920_dp+1.31649_dp*sin2t)  &
                      + mf2*4.52586067320_dp
        ddp2(i,j,2) =  1.44933285154_dp + mf1 * 2.01311801180_dp  + mf2*0.97512223853_dp
        ddp2(i,j,3) =  0.35569769252_dp + mf1*(-7.37249576667_dp) + mf2*(-12.2810377713_dp)
        ddp2(i,j,4) = -2.06533084541_dp + mf1 * 8.23741345333_dp  + mf2*5.93975747420_dp

        ddp4(i,j,0) =  0.21879385627_dp + mf1*(-0.58731641193_dp) + mf2*3.48695755800_dp
        ddp4(i,j,1) = -1.18964307357_dp + mf1 * 1.24891317047_dp  + mf2*(-14.9159739347_dp)
        ddp4(i,j,2) =  1.16268885692_dp + mf1*(-0.50852797392_dp) + mf2*15.3720218600_dp
        ddp4(i,j,3) =  0.0_dp
        ddp4(i,j,4) =  0.0_dp

        do k = 1, ncomp

           mdd_k = min( mseg(k), 2.0_dp )

           msegij = ( mdd_i * mdd_j * mdd_k )**(1.0_dp/3.0_dp)
           mf1 = ( msegij - 1.0_dp ) / msegij
           mf2 = mf1 * ( msegij - 2.0_dp ) / msegij
           ddp3(i,j,k,0) = -0.06467735252_dp + mf1*(-0.95208758351_dp+0.28503_dp*sin2t)  &
                + mf2*(-0.62609792333_dp+2.2195*sin2t)
           ddp3(i,j,k,1) =  0.19758818347_dp + mf1 * 2.99242575222_dp  + mf2*1.29246858189_dp
           ddp3(i,j,k,2) = -0.80875619458_dp + mf1*(-2.38026356489_dp) + mf2*1.65427830900_dp
           ddp3(i,j,k,3) =  0.69028490492_dp + mf1*(-0.27012609786_dp) + mf2*(-3.43967436378_dp)
           ddp3(i,j,k,4) =  0.0_dp

        end do

     end do

  end do

end subroutine load_dd_constants


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine qq_constants
!> \brief  provides the constants of the quadrupole-term
! This subroutine provides the constants of the quadrupole-term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine load_qq_constants

  use BASIC_VARIABLES, only: ncomp
  use pcsaft_pure_and_binary_parameters, only: mseg

  !-----------------------------------------------------------------------------
  integer                                    :: i, j, k
  real(dp)                                   :: mqq_i, mqq_j, mqq_k
  real(dp)                                   :: mf1, mf2, msegij
  !-----------------------------------------------------------------------------

  do i = 1, ncomp

     mqq_i = min( mseg(i), 2.0_dp )

     do j = 1, ncomp

        mqq_j = min( mseg(j), 2.0_dp )

        msegij = ( mqq_i * mqq_j )**0.5_dp
        mf1 = ( msegij - 1.0_dp ) / msegij
        mf2 = mf1 * ( msegij - 2.0_dp ) / msegij

        qqp2(i,j,0) =  1.237830788_dp + mf1 * 1.285410878_dp  + mf2*1.794295401_dp
        qqp2(i,j,1) =  2.435503144_dp + mf1*(-11.46561451_dp) + mf2*0.769510293_dp
        qqp2(i,j,2) =  1.633090469_dp + mf1 *22.08689285_dp   + mf2*7.264792255_dp
        qqp2(i,j,3) = -1.611815241_dp + mf1 * 7.46913832_dp   + mf2*94.48669892_dp
        qqp2(i,j,4) =  6.977118504_dp + mf1*(-17.19777208_dp) + mf2*(-77.1484579_dp)

        qqp4(i,j,0) =  0.454271755_dp + mf1*(-0.813734006_dp) + mf2*6.868267516_dp
        qqp4(i,j,1) = -4.501626435_dp + mf1 * 10.06402986_dp  + mf2*(-5.173223765_dp)
        qqp4(i,j,2) =  3.585886783_dp + mf1*(-10.87663092_dp) + mf2*(-17.2402066_dp)
        qqp4(i,j,3) =  0.0_dp
        qqp4(i,j,4) =  0.0_dp

        do k = 1, ncomp

           mqq_k = min( mseg(k), 2.0_dp )

           msegij = ( mqq_i * mqq_j * mqq_k )**(1.0_dp/3.0_dp)
           mf1 = ( msegij - 1.0_dp ) / msegij
           mf2 = mf1 * ( msegij - 2.0_dp ) / msegij
           qqp3(i,j,k,0) = -0.500043713_dp + mf1 * 2.000209381_dp + mf2*3.135827145_dp
           qqp3(i,j,k,1) =  6.531869153_dp + mf1*(-6.78386584_dp) + mf2*7.247588801_dp
           qqp3(i,j,k,2) = -16.01477983_dp + mf1 * 20.38324603_dp + mf2*3.075947834_dp
           qqp3(i,j,k,3) =  14.42597018_dp + mf1*(-10.89598394_dp)
           qqp3(i,j,k,4) =  0.0_dp

        end do

     end do

  end do

end subroutine load_qq_constants

end module eos_constants
