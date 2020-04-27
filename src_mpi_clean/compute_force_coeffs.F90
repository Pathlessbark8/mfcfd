module compute_force_coeffs_mod
    !#include <petsc/finclude/petscsys.h>
    
    use data_structure_mod
    use petsc_data_structure_mod
    
    contains
    
    
    
    subroutine compute_cl_cd_cm()
        
        
        implicit none
        
        integer :: i, j, k
        integer :: l, m, r
        real*8 :: cp, temp
        real*8 :: lx, ly, mx, my, rx, ry
        real*8 :: ds1, ds2, ds
        
        real*8, dimension(shapes) :: H, V, pitch_mom
        real*8, dimension(shapes) :: lCd, lCm
        real*8 :: nx, ny
        character(len=64) :: cp_file
        character(len=10) :: itos
        ! PetscErrorCode :: ierr
        
        cp_file = 'cp/'//'cp-file'
        if (proc>1) cp_file = 'cp/'//'cp-file'//trim(itos(4,rank))
        
        OPEN(UNIT=201,FILE=trim(cp_file),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
        !if(rank==0) OPEN(UNIT=202,FILE='clcdcm')
        
        temp = 0.5d0*rho_inf*Mach*Mach
        
        H = 0.d0
        V = 0.d0
        pitch_mom = 0.0d0
        
        do j = 1, shape_points
            
            m = shape_points_index(j)
            r = point%right(m) 
            l = point%left(m) 
            
            
            lx = point%x(l)
            ly = point%y(l)
            
            mx = point%x(m)
            my = point%y(m)
            
            rx = point%x(r)
            ry = point%y(r)
            
            ds1 = (mx - lx)**2 + (my - ly)**2
            ds1 = dsqrt(ds1)
            
            ds2 = (rx - mx)**2 + (ry - my)**2
            ds2 = dsqrt(ds2)
            
            ds = 0.5d0*(ds1 + ds2)
            
            nx = point%nx(m)
            ny = point%ny(m)
            
            cp = point%prim(4,m) - pr_inf
            cp = -cp/temp
            
            write(201, *) point%flag_2(m), point%x(m), cp
            
            H(point%flag_2(m)) = H(point%flag_2(m)) + cp*nx*ds
            V(point%flag_2(m)) = V(point%flag_2(m)) + cp*ny*ds
            
            pitch_mom(point%flag_2(m)) = pitch_mom(point%flag_2(m))&
            + (-cp*ny*ds*(mx - 0.25d0) + cp*nx*ds*(my))
            
        enddo
        
        lCl = V*dcos(theta) - H*dsin(theta)
        lCd = H*dcos(theta) + V*dsin(theta)
        lCm = pitch_mom
        
        ! call MPI_Reduce(lCl, Cl , shapes, MPI_DOUBLE, MPI_SUM, 0, &
        !         & PETSC_COMM_WORLD, ierr)
        ! call MPI_Reduce(lCd, Cd , shapes, MPI_DOUBLE, MPI_SUM, 0, &
        !         & PETSC_COMM_WORLD, ierr)
        ! call MPI_Reduce(lCm, Cm , shapes, MPI_DOUBLE, MPI_SUM, 0, &
        !         & PETSC_COMM_WORLD, ierr)
        
        ! if(rank == 0) then
        !         do j = 1, shapes
        !                 write(*,'(i4,3e30.20)') j, Cl, Cd, Cm
        !         end do
        ! end if
        
        CLOSE(UNIT=201)
        
    end subroutine 
    
end module compute_force_coeffs_mod
