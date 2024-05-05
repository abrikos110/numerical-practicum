function f(x)
    implicit none

    real, intent(in) :: x
    real :: f

    f = (x - 0.5) * (x - 0.7) * (x - 0.9)
    f = (3 * x**2 - 1) * exp(-x**2)
end function f


function calc_poly(a, n, x)
    implicit none

    integer, intent(in) :: n
    real, intent(in) :: x, a(n)
    real :: calc_poly
    integer :: i

    calc_poly = 0.0
    do i = n, 1, -1
        calc_poly = calc_poly * x
        calc_poly = calc_poly + a(i)
    end do
end function calc_poly


function mse(a, n_a, x1, x2, n_x)
    implicit none

    integer, intent(in) :: n_a, n_x
    real, intent(in) :: a(n_a), x1, x2

    real :: mse, calc_poly, f, q
    integer :: i

    mse = 0.0
    do i = 0, n_x
        q = x1 + (x2 - x1) / n_x * i
        mse = mse + (calc_poly(a, n_a, q) - f(q)) ** 2
    end do
end function mse


subroutine mse_grad(a, n_a, x1, x2, n_x, ans)
    implicit none

    integer, intent(in) :: n_a, n_x
    real, intent(in) :: a(n_a), x1, x2

    real, intent(out) :: ans(n_a)
    real :: calc_poly, f, q
    integer :: i, v

    ans = 0.0
    do v = 1, n_a
        do i = 0, n_x
            q = x1 + (x2 - x1) / n_x * i
            ans(v) = ans(v) + 2 * (calc_poly(a, n_a, q) - f(q)) * q**(v-1)
        end do
    end do
end subroutine mse_grad


function extremum_of_quadratic(fm1, f0, f1)
    implicit none

    real, intent(in) :: fm1, f0, f1
    real :: extremum_of_quadratic

    extremum_of_quadratic = (fm1 - f1) / 2 / (fm1 + f1 - 2 * f0)
end function extremum_of_quadratic


function dot(a, b, n)
    implicit none

    integer :: n
    real, intent(in) :: a(n), b(n)
    real :: dot

    dot = sum(a*b)
end function dot


subroutine approx(a, n_a, c, d, n_x, eps)
    implicit none

    integer, intent(in) :: n_a, n_x
    real, intent(in) :: c, d, eps
    real :: x, b, a(n_a), v(n_a), g(n_a), dot, mse, extremum_of_quadratic

    call mse_grad(a, n_a, c, d, n_x, g)
    v = 0

    write (*, *) "polynomial coeffs", a
    do while (1==1)
        b = dot(g, g, n_a)
        call mse_grad(a, n_a, c, d, n_x, g)
        b = dot(g, g, n_a) / b
        v = -g + v * b

        x = extremum_of_quadratic(mse(a-v, n_a, c, d, n_x), mse(a, n_a, c, d, n_x), mse(a+v, n_a, c, d, n_x))
        a = a + x * v

        !write (*, *) "polynomial coeffs", a
        !write (*, *) "|grad|", sqrt(dot(g, g, n_a))
        !write (*, *) "mse", mse(a, n_a, c, d, n_x)
        !write (*, *)
        write (*, *) mse(a, n_a, c, d, n_x), sqrt(dot(g, g, n_a)), maxval(abs(x*v)), maxval(abs(x*v)) / maxval(abs(a))

        if (dot(v, v, n_a) < eps*eps) then
            exit
        end if
    end do
    write (*, *) "polynomial coeffs", a

end subroutine approx


program prog
    implicit none

    real, allocatable :: a(:), g(:), v(:)
    real :: f, calc_poly, mse, b, dot, c, d, x, extremum_of_quadratic, eps

    integer :: N, n0, m, n_a, i  ! n1 = n0 + m, deg(P) in [n0 : n1)

    c = 0.1
    d = 1.1
    N = 25
    n0 = 3
    m = 2

    write (*, '(A)', advance='no') "Enter c: "
    read (*,*) c

    write (*, '(A)', advance='no') "Enter d: "
    read (*,*) d

    write (*, '(A)', advance='no') "Enter N: "
    read (*,*) N

    write (*, '(A)', advance='no') "Enter n0: "
    read (*,*) n0

    write (*, '(A)', advance='no') "Enter m (n1 = n0 + m): "
    read (*,*) m

    allocate(a(n0 + m))
    allocate(g(n0 + m))
    allocate(v(n0 + m))

    a = 0

    do n_a = n0+1, n0+m
        call mse_grad(a, n_a, c, d, N, g) ! DO NOT REMOVE THIS LINE ===============

        write (*, *)
        write (*, '(a, f5.2, a, f5.2, a)') "[", c, ":", d, "]"
        write (*, '(a, i2)') "N = ", N
        write (*, '(a, i2, a, i2, a)') "deg(P) = ", n_a
        write (*, *) mse(a, n_a, c, d, N)
        write (*, *)


        v = 0
        eps = 1e-3
        call approx(a, n_a, c, d, N, eps)

        write (*, *) "================="
        write (*, *) "x_i            ", "f(x_i)           ", "P(x_i)          ", " |f(x_i) - P(x_i)|"
        do i = 0, N
            x = i * (d - c) / N + c
            b = calc_poly(a, n_a, x)
            write (*, *) x, f(x), b, abs(b - f(x))
        end do
    end do

    write (*, *) 0.1

    deallocate(a)
    deallocate(g)
    deallocate(v)
end program prog
