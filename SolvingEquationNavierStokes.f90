!******************** 
! Common program    *
! 7/10/16           *
! Annotated by I.V. *
!********************

Program m2
 implicit none
 integer, parameter::n=20 , m=20 ! сетка 20х20
 integer i ,j, k
 real*8 F(n,m) ,eps, dt, h,  h2, x(n), y(m) ,Omega(n,m), Omega_n(n,m) ,Py,Ox,Px,Oy,L,h22
 open(1,file='1.txt')

! 5 массивов: 
! Три двумерных - F(n,m), Omega(n,m), Omega_n(n,m)
! Два одномерных x(n), y(m)
 
 F=0d0 ! занулили массив F(20,20)

 h=1d0/(n-1) ! шаг сетки
 h2=h*h      ! квадрат шага
 eps=1d-6    ! точность
 dt=0.2d0*h2 ! время
 h22=h*2d0   

 x(1)=0d0    ! в массиве х(20) - первый элемент нулевой
 y(1)=0d0    ! в массиве у(20) - первый элемент нулевой

   ! х(20) = {0, h,h,h,h,...h}
   do i=2,n
    x(i)=x(i-1)+h
   enddo 
   
   ! y(20) = {0, h,h,h,h,...h}
   do j=2,m	 
    y(j)=y(j-1)+h
   enddo 
    
! обнулили массивы Omega(n,m) и F(n,m)
 Omega=0d0
 F=0d0

 do k=1,1000 ! внешний цикл по времени
  
  do j=2,m-1 ! шаги по игрек
        
       do i=2,n-1 ! шаги по икс
            !*** свертка производных************ 
	    Py=(F(i,j+1)-F(i,j-1))/h22
            Ox=(Omega(i+1,j)-Omega(i-1,j))/h22
            Px=(F(i+1,j)-F(i-1,j))/h22
            Oy=(Omega(i,j+1)-Omega(i,j-1))/h22
            !***********************************  
            L=(Omega(i+1,j)-4d0*Omega(i,j)+Omega(i-1,j)+Omega(i,j+1)+Omega(i,j-1))/h2
	    Omega_n(i,j)=Omega(i,j)+dt*(L-Py*Ox+Px*Oy)
       enddo

  enddo  
	 
  call Phi(Omega_n,F,n,m,eps) ! вызов процедуры /* ?решаем здесь уравнение для Omega_n? */

  ! граничные условия для Omega_n
  Omega_n(1,:)=-2*F(2,:)/h2         ! первую строку заполняем тем, что справа
  Omega_n(n,:)=-2*F(n-1,:)/h2       ! n-ую строку заполняем тем, что справа
  Omega_n(:,1)=-2*F(:,2)/h2         ! первый столбец тем, что справа
  Omega_n(:,m)=-2*(F(:,m-1)+h)/h2   ! столбец m тем, что справа
  Omega=Omega_n
 enddo ! конец цикла по времени (k=1..1000)  
	
  ! запись массивов x(n), y(m), F(n,m) в файл // warum denn *x, *y,**F   
  do j=1,m
     do i=1,n
	write(1,*) x(i),y(j),F(i,j)
     enddo
  enddo

  close(1)
  pause
  end


!   сетка и --------> г.у.
! 
!         v=1 
!    *************
!    *           *
!v=0 *           * v=0
!    *           *
!    *           *
!    *************
!         v=0
! 
! --------------> должно обеспечиваться равенство нулю касательных и тангенциальных компонент

! программа на входе Omega, F(которая на самом деле Psi )
subroutine Phi(Omega,F,n,m,eps)
implicit none
integer i, j, n, m
real*8 F(n,m), F_n(n,m) ,eps ,dt, h, Fs, h2, Fs_n, Omega(n,m)

h=1d0/(n-1)
h2=h*h
dt=0.2d0*h2

Fs=1d0
Fs_n=0d0

do while (abs(Fs_n-Fs)>abs(eps*Fs))
 
 do j=2,m-1
    do i=2,n-1
       F_n(i,j)= F(i,j)+(dt/h2)*(F(i+1,j)+F(i,j+1)+F(i-1,j)+F(i,j-1)-4d0*F(i,j)) + dt*Omega(i,j)
    enddo
  enddo

    ! граничные условия на Psi. Скорость? везде нуль
    F_n(:,1)=0d0
    F_n(1,:)=0d0
    F_n(n,:)=0d0
    F_n(:,m)=0d0
    Fs = Fs_n  	
    Fs_n = 0d0

	do j=1,m
           do i=1,n
	   Fs_n=Fs_n+F_n(i,j)
	   enddo
	enddo
	Fs_n = Fs_n/n/m
  	F=F_n
	
enddo
print*, Fs_n
end subroutine Phi
