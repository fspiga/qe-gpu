!
! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE global_version
  !
  IMPLICIT NONE
  !
  SAVE
  !
  CHARACTER (LEN=6) :: version_number = '6.1'
  CHARACTER (LEN=12) :: svn_revision = '13369'
  !
  CHARACTER (LEN=6) :: gpu_version_number = '1.0'
  !
END MODULE global_version
