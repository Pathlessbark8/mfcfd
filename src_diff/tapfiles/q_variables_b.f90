!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
!
MODULE Q_VARIABLES_MOD_DIFF
  USE DATA_STRUCTURE_MOD_DIFF
  IMPLICIT NONE

CONTAINS
!  Differentiation of eval_q_variables in reverse (adjoint) mode (with options fixinterface):
!   gradient     of useful results: point.prim point.q
!   with respect to varying inputs: point.prim point.q
  SUBROUTINE EVAL_Q_VARIABLES_B()
    IMPLICIT NONE
    INTEGER :: k
    REAL*8 :: rho, u1, u2, pr, beta
    REAL*8 :: rhob, u1b, u2b, prb, betab
    REAL*8 :: two_times_beta
    REAL*8 :: two_times_betab
    INTRINSIC DLOG
    REAL*8 :: tempb
    DO k=max_points,1,-1
      two_times_betab = -pointb%q(4, k)
      pointb%q(4, k) = 0.0_8
      pr = point%prim(4, k)
      rho = point%prim(1, k)
      beta = 0.5d0*rho/pr
      two_times_beta = 2.0d0*beta
      u2 = point%prim(3, k)
      two_times_betab = two_times_betab + u2*pointb%q(3, k)
      u2b = two_times_beta*pointb%q(3, k)
      pointb%q(3, k) = 0.0_8
      u1 = point%prim(2, k)
      two_times_betab = two_times_betab + u1*pointb%q(2, k)
      u1b = two_times_beta*pointb%q(2, k)
      pointb%q(2, k) = 0.0_8
      betab = (2.5d0/beta-u2**2-u1**2)*pointb%q(1, k) + 2.0d0*&
&       two_times_betab
      u1b = u1b - beta*2*u1*pointb%q(1, k)
      u2b = u2b - beta*2*u2*pointb%q(1, k)
      tempb = 0.5d0*betab/pr
      rhob = tempb + pointb%q(1, k)/rho
      pointb%q(1, k) = 0.0_8
      prb = -(rho*tempb/pr)
      pointb%prim(4, k) = pointb%prim(4, k) + prb
      pointb%prim(3, k) = pointb%prim(3, k) + u2b
      pointb%prim(2, k) = pointb%prim(2, k) + u1b
      pointb%prim(1, k) = pointb%prim(1, k) + rhob
    END DO
  END SUBROUTINE EVAL_Q_VARIABLES_B

  SUBROUTINE EVAL_Q_VARIABLES()
    IMPLICIT NONE
    INTEGER :: k
    REAL*8 :: rho, u1, u2, pr, beta
    REAL*8 :: two_times_beta
    INTRINSIC DLOG
!
    DO k=1,max_points
      rho = point%prim(1, k)
      u1 = point%prim(2, k)
      u2 = point%prim(3, k)
      pr = point%prim(4, k)
      beta = 0.5d0*rho/pr
      point%q(1, k) = DLOG(rho) + DLOG(beta)*2.5d0 - beta*(u1*u1+u2*u2)
      two_times_beta = 2.0d0*beta
      point%q(2, k) = two_times_beta*u1
      point%q(3, k) = two_times_beta*u2
      point%q(4, k) = -two_times_beta
    END DO
  END SUBROUTINE EVAL_Q_VARIABLES

!  Differentiation of eval_q_derivatives in reverse (adjoint) mode (with options fixinterface):
!   gradient     of useful results: point.x point.y point.q point.dq
!                point.qm
!   with respect to varying inputs: point.x point.y point.q point.dq
!                point.qm
  SUBROUTINE EVAL_Q_DERIVATIVES_B()
    IMPLICIT NONE
    INTEGER :: i, k, r, nbh
    REAL*8 :: x_i, y_i, x_k, y_k
    REAL*8 :: x_ib, y_ib, x_kb, y_kb
    REAL*8 :: delx, dely, dist, weights
    REAL*8 :: delxb, delyb, distb, weightsb
    REAL*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
    REAL*8 :: sum_delx_sqrb, sum_dely_sqrb, sum_delx_delyb
    REAL*8 :: sum_delx_delq(4), sum_dely_delq(4)
    REAL*8 :: sum_delx_delqb(4), sum_dely_delqb(4)
    REAL*8 :: det, delq, temp
    REAL*8 :: detb
    REAL*8 :: one_by_det
    REAL*8 :: one_by_detb
    INTRINSIC DSQRT
    REAL*8 :: tempb
    REAL*8, DIMENSION(4) :: tempb0
    REAL*8 :: tempb1
    REAL*8 :: tempb2
    REAL*8, DIMENSION(4) :: tempb3
    REAL*8, DIMENSION(4) :: tempb4
    REAL*8, DIMENSION(4) :: tempb5
    INTEGER :: branch
    INTEGER :: ad_to
    DO i=1,local_points
      x_i = point%x(i)
      y_i = point%y(i)
      CALL PUSHREAL8(sum_delx_sqr)
      sum_delx_sqr = 0.d0
      CALL PUSHREAL8(sum_dely_sqr)
      sum_dely_sqr = 0.d0
      CALL PUSHREAL8(sum_delx_dely)
      sum_delx_dely = 0.d0
      CALL PUSHREAL8ARRAY(sum_delx_delq, 4)
      sum_delx_delq = 0.d0
      CALL PUSHREAL8ARRAY(sum_dely_delq, 4)
      sum_dely_delq = 0.d0
      CALL PUSHREAL8ARRAY(point%qm(1, :, i), 4)
      point%qm(1, :, i) = point%q(:, i)
      CALL PUSHREAL8ARRAY(point%qm(2, :, i), 4)
      point%qm(2, :, i) = point%q(:, i)
      DO k=1,point%nbhs(i)
        nbh = point%conn(i, k)
        DO r=1,4
          IF (point%q(r, nbh) .GT. point%qm(1, r, i)) THEN
            CALL PUSHREAL8(point%qm(1, r, i))
            point%qm(1, r, i) = point%q(r, nbh)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (point%q(r, nbh) .LT. point%qm(2, r, i)) THEN
            CALL PUSHREAL8(point%qm(2, r, i))
            point%qm(2, r, i) = point%q(r, nbh)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        x_k = point%x(nbh)
        y_k = point%y(nbh)
        delx = x_k - x_i
        dely = y_k - y_i
        CALL PUSHREAL8(dist)
        dist = DSQRT(delx*delx + dely*dely)
        weights = dist**power
        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
        sum_dely_sqr = sum_dely_sqr + dely*dely*weights
        sum_delx_dely = sum_delx_dely + delx*dely*weights
        sum_delx_delq = sum_delx_delq + weights*delx*(point%q(:, nbh)-&
&         point%q(:, i))
        sum_dely_delq = sum_dely_delq + weights*dely*(point%q(:, nbh)-&
&         point%q(:, i))
      END DO
      CALL PUSHINTEGER4(k - 1)
    END DO
    DO i=local_points,1,-1
      det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
      one_by_det = 1.0d0/det
      sum_delx_delqb = 0.0_8
      sum_dely_delqb = 0.0_8
      tempb4 = one_by_det*pointb%dq(2, :, i)
      one_by_detb = SUM((sum_dely_delq*sum_delx_sqr-sum_delx_delq*&
&       sum_delx_dely)*pointb%dq(2, :, i))
      pointb%dq(2, :, i) = 0.0_8
      tempb5 = one_by_det*pointb%dq(1, :, i)
      sum_dely_delqb = sum_delx_sqr*tempb4 - sum_delx_dely*tempb5
      sum_delx_delqb = sum_dely_sqr*tempb5 - sum_delx_dely*tempb4
      one_by_detb = one_by_detb + SUM((sum_delx_delq*sum_dely_sqr-&
&       sum_dely_delq*sum_delx_dely)*pointb%dq(1, :, i))
      pointb%dq(1, :, i) = 0.0_8
      detb = -(one_by_detb/det**2)
      sum_delx_sqrb = sum_dely_sqr*detb + SUM(sum_dely_delq*tempb4)
      sum_delx_delyb = -SUM(sum_dely_delq*tempb5) - 2*sum_delx_dely*detb&
&       - SUM(sum_delx_delq*tempb4)
      sum_dely_sqrb = sum_delx_sqr*detb + SUM(sum_delx_delq*tempb5)
      y_i = point%y(i)
      x_i = point%x(i)
      y_ib = 0.0_8
      x_ib = 0.0_8
      CALL POPINTEGER4(ad_to)
      DO k=ad_to,1,-1
        weights = dist**power
        nbh = point%conn(i, k)
        y_k = point%y(nbh)
        dely = y_k - y_i
        tempb = SUM((point%q(:, nbh)-point%q(:, i))*sum_dely_delqb)
        tempb0 = weights*dely*sum_dely_delqb
        pointb%q(:, nbh) = pointb%q(:, nbh) + tempb0
        pointb%q(:, i) = pointb%q(:, i) - tempb0
        x_k = point%x(nbh)
        delx = x_k - x_i
        tempb1 = SUM((point%q(:, nbh)-point%q(:, i))*sum_delx_delqb)
        weightsb = delx*tempb1 + dely**2*sum_dely_sqrb + delx**2*&
&         sum_delx_sqrb + delx*dely*sum_delx_delyb + dely*tempb
        tempb3 = weights*delx*sum_delx_delqb
        pointb%q(:, nbh) = pointb%q(:, nbh) + tempb3
        pointb%q(:, i) = pointb%q(:, i) - tempb3
        IF (dist .LE. 0.0 .AND. (power .EQ. 0.0 .OR. power .NE. INT(&
&           power))) THEN
          distb = 0.0
        ELSE
          distb = power*dist**(power-1)*weightsb
        END IF
        CALL POPREAL8(dist)
        IF (delx**2 + dely**2 .EQ. 0.0) THEN
          tempb2 = 0.0
        ELSE
          tempb2 = distb/(2.D0*DSQRT(delx**2+dely**2))
        END IF
        delyb = weights*delx*sum_delx_delyb + 2*dely*tempb2 + weights*2*&
&         dely*sum_dely_sqrb + weights*tempb
        delxb = weights*dely*sum_delx_delyb + 2*delx*tempb2 + weights*2*&
&         delx*sum_delx_sqrb + weights*tempb1
        y_kb = delyb
        y_ib = y_ib - delyb
        x_kb = delxb
        x_ib = x_ib - delxb
        pointb%y(nbh) = pointb%y(nbh) + y_kb
        pointb%x(nbh) = pointb%x(nbh) + x_kb
        DO r=4,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            CALL POPREAL8(point%qm(2, r, i))
            pointb%q(r, nbh) = pointb%q(r, nbh) + pointb%qm(2, r, i)
            pointb%qm(2, r, i) = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(point%qm(1, r, i))
            pointb%q(r, nbh) = pointb%q(r, nbh) + pointb%qm(1, r, i)
            pointb%qm(1, r, i) = 0.0_8
          END IF
        END DO
      END DO
      CALL POPREAL8ARRAY(point%qm(2, :, i), 4)
      pointb%q(:, i) = pointb%q(:, i) + pointb%qm(2, :, i)
      pointb%qm(2, :, i) = 0.0_8
      CALL POPREAL8ARRAY(point%qm(1, :, i), 4)
      pointb%q(:, i) = pointb%q(:, i) + pointb%qm(1, :, i)
      pointb%qm(1, :, i) = 0.0_8
      CALL POPREAL8ARRAY(sum_dely_delq, 4)
      CALL POPREAL8ARRAY(sum_delx_delq, 4)
      CALL POPREAL8(sum_delx_dely)
      CALL POPREAL8(sum_dely_sqr)
      CALL POPREAL8(sum_delx_sqr)
      pointb%y(i) = pointb%y(i) + y_ib
      pointb%x(i) = pointb%x(i) + x_ib
    END DO
  END SUBROUTINE EVAL_Q_DERIVATIVES_B

  SUBROUTINE EVAL_Q_DERIVATIVES()
    IMPLICIT NONE
    INTEGER :: i, k, r, nbh
    REAL*8 :: x_i, y_i, x_k, y_k
    REAL*8 :: delx, dely, dist, weights
    REAL*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
    REAL*8 :: sum_delx_delq(4), sum_dely_delq(4)
    REAL*8 :: det, delq, temp
    REAL*8 :: one_by_det
    INTRINSIC DSQRT
    DO i=1,local_points
      x_i = point%x(i)
      y_i = point%y(i)
      sum_delx_sqr = 0.d0
      sum_dely_sqr = 0.d0
      sum_delx_dely = 0.d0
      sum_delx_delq = 0.d0
      sum_dely_delq = 0.d0
      point%qm(1, :, i) = point%q(:, i)
      point%qm(2, :, i) = point%q(:, i)
      DO k=1,point%nbhs(i)
        nbh = point%conn(i, k)
        DO r=1,4
          IF (point%q(r, nbh) .GT. point%qm(1, r, i)) point%qm(1, r, i)&
&            = point%q(r, nbh)
          IF (point%q(r, nbh) .LT. point%qm(2, r, i)) point%qm(2, r, i)&
&            = point%q(r, nbh)
        END DO
        x_k = point%x(nbh)
        y_k = point%y(nbh)
        delx = x_k - x_i
        dely = y_k - y_i
        dist = DSQRT(delx*delx + dely*dely)
        weights = dist**power
        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
        sum_dely_sqr = sum_dely_sqr + dely*dely*weights
        sum_delx_dely = sum_delx_dely + delx*dely*weights
        sum_delx_delq = sum_delx_delq + weights*delx*(point%q(:, nbh)-&
&         point%q(:, i))
        sum_dely_delq = sum_dely_delq + weights*dely*(point%q(:, nbh)-&
&         point%q(:, i))
      END DO
      det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
      one_by_det = 1.0d0/det
      point%dq(1, :, i) = (sum_delx_delq*sum_dely_sqr-sum_dely_delq*&
&       sum_delx_dely)*one_by_det
      point%dq(2, :, i) = (sum_dely_delq*sum_delx_sqr-sum_delx_delq*&
&       sum_delx_dely)*one_by_det
    END DO
  END SUBROUTINE EVAL_Q_DERIVATIVES

!  Differentiation of qtilde_to_primitive in reverse (adjoint) mode (with options fixinterface):
!   gradient     of useful results: u1 u2 pr rho
!   with respect to varying inputs: qtilde
  SUBROUTINE QTILDE_TO_PRIMITIVE_B(qtilde, qtildeb, u1, u1b, u2, u2b, &
&   rho, rhob, pr, prb)
    IMPLICIT NONE
    REAL*8 :: qtilde(4), u1, u2, rho, pr
    REAL*8 :: qtildeb(4), u1b, u2b, rhob, prb
    REAL*8 :: beta, temp, temp1, temp2
    REAL*8 :: betab, tempb, temp1b, temp2b
    REAL*8 :: q1, q2, q3, q4
    REAL*8 :: q1b, q2b, q3b, q4b
    INTRINSIC DLOG
    INTRINSIC DEXP
    q1 = qtilde(1)
    q2 = qtilde(2)
    q3 = qtilde(3)
    q4 = qtilde(4)
    beta = -(q4*0.5d0)
    temp = 0.5d0/beta
    CALL PUSHREAL8(u1)
    u1 = q2*temp
    CALL PUSHREAL8(u2)
    u2 = q3*temp
    temp1 = q1 + beta*(u1*u1+u2*u2)
    temp2 = temp1 - DLOG(beta)/(gamma-1)
    CALL PUSHREAL8(rho)
    rho = DEXP(temp2)
    rhob = rhob + temp*prb
    temp2b = DEXP(temp2)*rhob
    temp1b = temp2b
    q1b = temp1b
    u1b = u1b + beta*2*u1*temp1b
    u2b = u2b + beta*2*u2*temp1b
    tempb = q3*u2b + q2*u1b + rho*prb
    CALL POPREAL8(rho)
    betab = (u1**2+u2**2)*temp1b - 0.5d0*tempb/beta**2 - temp2b/((gamma-&
&     1)*beta)
    CALL POPREAL8(u2)
    q3b = temp*u2b
    CALL POPREAL8(u1)
    q2b = temp*u1b
    q4b = -(0.5d0*betab)
    qtildeb = 0.0_8
    qtildeb(4) = qtildeb(4) + q4b
    qtildeb(3) = qtildeb(3) + q3b
    qtildeb(2) = qtildeb(2) + q2b
    qtildeb(1) = qtildeb(1) + q1b
  END SUBROUTINE QTILDE_TO_PRIMITIVE_B

  SUBROUTINE QTILDE_TO_PRIMITIVE(qtilde, u1, u2, rho, pr)
    IMPLICIT NONE
    REAL*8 :: qtilde(4), u1, u2, rho, pr
    REAL*8 :: beta, temp, temp1, temp2
    REAL*8 :: q1, q2, q3, q4
    INTRINSIC DLOG
    INTRINSIC DEXP
    q1 = qtilde(1)
    q2 = qtilde(2)
    q3 = qtilde(3)
    q4 = qtilde(4)
    beta = -(q4*0.5d0)
    temp = 0.5d0/beta
    u1 = q2*temp
    u2 = q3*temp
    temp1 = q1 + beta*(u1*u1+u2*u2)
    temp2 = temp1 - DLOG(beta)/(gamma-1)
    rho = DEXP(temp2)
    pr = rho*temp
  END SUBROUTINE QTILDE_TO_PRIMITIVE

END MODULE Q_VARIABLES_MOD_DIFF

