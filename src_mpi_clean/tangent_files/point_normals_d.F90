!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7259) - 18 Jan 2019 09:36
!
MODULE POINT_NORMALS_MOD_DIFF
  USE DATA_STRUCTURE_MOD_DIFF
  IMPLICIT NONE

CONTAINS
!  Differentiation of compute_normals in forward (tangent) mode (with options fixinterface):
!   variations   of useful results: *(point.nx) *(point.ny)
!   with respect to varying inputs: *(point.x) *(point.y)
!   Plus diff mem management of: point.x:in point.y:in point.nx:in
!                point.ny:in
  SUBROUTINE COMPUTE_NORMALS_D()
    IMPLICIT NONE
    DOUBLE PRECISION :: lx, ly, mx, my, rx, ry
    DOUBLE PRECISION :: lxd, lyd, mxd, myd, rxd, ryd
    DOUBLE PRECISION :: nx1, nx2, ny1, ny2, nx, ny
    DOUBLE PRECISION :: nx1d, nx2d, ny1d, ny2d, nxd, nyd
    DOUBLE PRECISION :: det
    DOUBLE PRECISION :: detd
    INTEGER :: i, j, k, l, m, r
    INTRINSIC DSQRT
    DOUBLE PRECISION :: arg1
    DOUBLE PRECISION :: arg1d
    pointd%nx = 0.0_8
    pointd%ny = 0.0_8
!Finding the normals for the points on the shapes ..   
    DO i=1,wall_points
      m = wall_points_index(i)
      l = point%left(m)
      r = point%right(m)
      lxd = pointd%x(l)
      lx = point%x(l)
      lyd = pointd%y(l)
      ly = point%y(l)
      mxd = pointd%x(m)
      mx = point%x(m)
      myd = pointd%y(m)
      my = point%y(m)
      rxd = pointd%x(r)
      rx = point%x(r)
      ryd = pointd%y(r)
      ry = point%y(r)
      nx1d = myd - lyd
      nx1 = my - ly
      nx2d = ryd - myd
      nx2 = ry - my
      ny1d = mxd - lxd
      ny1 = mx - lx
      ny2d = rxd - mxd
      ny2 = rx - mx
      nxd = 0.5d0*(nx1d+nx2d)
      nx = 0.5d0*(nx1+nx2)
      nyd = 0.5d0*(ny1d+ny2d)
      ny = 0.5d0*(ny1+ny2)
      arg1d = nxd*nx + nx*nxd + nyd*ny + ny*nyd
      arg1 = nx*nx + ny*ny
      IF (arg1 .EQ. 0.0) THEN
        detd = 0.D0
      ELSE
        detd = arg1d/(2.D0*DSQRT(arg1))
      END IF
      det = DSQRT(arg1)
      nxd = -((nxd*det-nx*detd)/det**2)
      nx = -(nx/det)
      nyd = (nyd*det-ny*detd)/det**2
      ny = ny/det
      pointd%nx(m) = nxd
      point%nx(m) = nx
      pointd%ny(m) = nyd
      point%ny(m) = ny
    END DO
!	Finding the normals for the outer boundary points ..
    DO i=1,outer_points
      m = outer_points_index(i)
      l = point%left(m)
      r = point%right(m)
      lxd = pointd%x(l)
      lx = point%x(l)
      lyd = pointd%y(l)
      ly = point%y(l)
      mxd = pointd%x(m)
      mx = point%x(m)
      myd = pointd%y(m)
      my = point%y(m)
      rxd = pointd%x(r)
      rx = point%x(r)
      ryd = pointd%y(r)
      ry = point%y(r)
      nx1d = myd - lyd
      nx1 = my - ly
      nx2d = ryd - myd
      nx2 = ry - my
      ny1d = mxd - lxd
      ny1 = mx - lx
      ny2d = rxd - mxd
      ny2 = rx - mx
      nxd = 0.5d0*(nx1d+nx2d)
      nx = 0.5d0*(nx1+nx2)
      nyd = 0.5d0*(ny1d+ny2d)
      ny = 0.5d0*(ny1+ny2)
      arg1d = nxd*nx + nx*nxd + nyd*ny + ny*nyd
      arg1 = nx*nx + ny*ny
      IF (arg1 .EQ. 0.0) THEN
        detd = 0.D0
      ELSE
        detd = arg1d/(2.D0*DSQRT(arg1))
      END IF
      det = DSQRT(arg1)
      nxd = -((nxd*det-nx*detd)/det**2)
      nx = -(nx/det)
      nyd = (nyd*det-ny*detd)/det**2
      ny = ny/det
      pointd%nx(m) = nxd
      point%nx(m) = nx
      pointd%ny(m) = nyd
      point%ny(m) = ny
    END DO
    IF (interior_points_normal_flag .EQ. 0 .AND. format .NE. 2) THEN
      DO i=1,interior_points
        k = interior_points_index(i)
        pointd%nx(k) = 0.0_8
        point%nx(k) = 0.d0
        pointd%ny(k) = 0.0_8
        point%ny(k) = 1.d0
      END DO
    ELSE IF (interior_points_normal_flag .EQ. 1 .AND. format .NE. 2) &
&   THEN
      DO i=1,interior_points
        k = interior_points_index(i)
        pointd%nx(k) = 0.0_8
        point%nx(k) = 1.d0
        pointd%ny(k) = 0.0_8
        point%ny(k) = 0.d0
      END DO
    END IF
  END SUBROUTINE COMPUTE_NORMALS_D

  SUBROUTINE COMPUTE_NORMALS()
    IMPLICIT NONE
    DOUBLE PRECISION :: lx, ly, mx, my, rx, ry
    DOUBLE PRECISION :: nx1, nx2, ny1, ny2, nx, ny
    DOUBLE PRECISION :: det
    INTEGER :: i, j, k, l, m, r
    INTRINSIC DSQRT
    DOUBLE PRECISION :: arg1
!Finding the normals for the points on the shapes ..   
    DO i=1,wall_points
      m = wall_points_index(i)
      l = point%left(m)
      r = point%right(m)
      lx = point%x(l)
      ly = point%y(l)
      mx = point%x(m)
      my = point%y(m)
      rx = point%x(r)
      ry = point%y(r)
      nx1 = my - ly
      nx2 = ry - my
      ny1 = mx - lx
      ny2 = rx - mx
      nx = 0.5d0*(nx1+nx2)
      ny = 0.5d0*(ny1+ny2)
      arg1 = nx*nx + ny*ny
      det = DSQRT(arg1)
      nx = -(nx/det)
      ny = ny/det
      point%nx(m) = nx
      point%ny(m) = ny
    END DO
!	Finding the normals for the outer boundary points ..
    DO i=1,outer_points
      m = outer_points_index(i)
      l = point%left(m)
      r = point%right(m)
      lx = point%x(l)
      ly = point%y(l)
      mx = point%x(m)
      my = point%y(m)
      rx = point%x(r)
      ry = point%y(r)
      nx1 = my - ly
      nx2 = ry - my
      ny1 = mx - lx
      ny2 = rx - mx
      nx = 0.5d0*(nx1+nx2)
      ny = 0.5d0*(ny1+ny2)
      arg1 = nx*nx + ny*ny
      det = DSQRT(arg1)
      nx = -(nx/det)
      ny = ny/det
      point%nx(m) = nx
      point%ny(m) = ny
    END DO
    IF (interior_points_normal_flag .EQ. 0 .AND. format .NE. 2) THEN
      DO i=1,interior_points
        k = interior_points_index(i)
        point%nx(k) = 0.d0
        point%ny(k) = 1.d0
      END DO
    ELSE IF (interior_points_normal_flag .EQ. 1 .AND. format .NE. 2) &
&   THEN
      DO i=1,interior_points
        k = interior_points_index(i)
        point%nx(k) = 1.d0
        point%ny(k) = 0.d0
      END DO
    END IF
  END SUBROUTINE COMPUTE_NORMALS

END MODULE POINT_NORMALS_MOD_DIFF

