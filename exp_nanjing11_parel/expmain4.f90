
    !***************************************************
    ! �����ݴŵ���
    !***************************************************
    Module ran_mod
    Implicit None
    ! ran return a uniform random number between 0-1
    ! normal return a normal distribution
    contains
    function ran()   !returns random number between 0 - 1
    implicit none
    integer , save :: flag = 0
    double precision :: ran
    if(flag==0) then
        call random_seed()
        flag = 1
    endif
    call random_number(ran)     ! built in fortran 90 random number function
    end function ran

    function normal(mean,sigma)
    implicit none
    integer :: flag
    double precision, parameter :: pi = 3.141592653589793239
    double precision :: u1, u2, y1, y2, normal, mean, sigma
    save flag
    data flag /0/
    u1 = ran(); u2 = ran()
    if (flag.eq.0) then
        y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
        normal = mean + sigma*y1
        flag = 1
    else
        y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
        normal = mean + sigma*y2
        flag = 0
    endif
    end function normal
    !The above codes are made in Fortran 90 language, if you have any question, you may write to sealin2008@hotmail.com
    End Module ran_mod


    program main
    use ran_mod

    use mpi
    integer:: ierr, rank, nproc_current

    integer:: ntc,nolayer,npara,ns,ns_real    ! number of time channel / layers / parameter per point / points
    integer:: counter,i,j,k  ! index
    integer:: tot_para,tot_ntc
    integer:: size0,size1,size2
    integer:: l,ka,ns_index,nturn1,nturn11

    REAL :: cpu_start, cpu_finish, elapsed_time
    real *8,allocatable:: rho_iter(:), hh_iter(:)

    real *8,allocatable:: rho_pro(:), hh_pro(:)

    real *8 temp0(1,1),pi,htt,t_st,t_ed, xr1, hr1, rt1, rr1

    real *8 para1,para2,para3,para4,para5,para6,para7

    real *8 rerror,rerror_pre,lam,epsi,eps

    real *8,allocatable::hz1_iter(:),hz1_iter_tmp(:)      ! dBz/dt HCP

    real *8,allocatable::time(:)
    real *8,allocatable::deltadobs(:,:),deltam(:,:)
    real *8,allocatable::jacobi(:,:),jacobi1(:,:)
    real *8,allocatable::EYE(:,:),b0(:,:)
    real *8,allocatable::Btemp(:,:), invBtemp(:,:)
    real *8,allocatable::Atemp(:,:),height(:,:),depth(:,:)
    real *8,allocatable::m_pre(:,:),m_next(:,:),Vobs(:,:),Vobs_tot(:,:)
    real *8,allocatable::deltaRp(:,:),deltaRd(:,:)
    real *8,allocatable::Wp(:,:), Wd(:,:),Rp(:,:),Rd(:,:)
    real *8,allocatable::s(:),e(:),work(:),u(:,:),v(:,:)

    real *8,allocatable:: point1set(:), point2set(:), point3set(:), point4set(:),ht(:)
    integer :: iTimes1, iTimes2, iTimes3, iTimes4, rate
    
  CALL system_clock(count_rate=rate)
  call SYSTEM_CLOCK(iTimes1)

    Data pi/3.1415926D0/

    open(35, file='parameter_settings.txt', status='old')

    read(35,*)para1,para2,para3,para4,para5,t_st,t_ed, xr1, hr1, rt1, rr1,para6,para7

    ntc = floor(para1)
    nolayer = floor(para2)
    ns_real = floor(para3)
    ns = floor(para4)
    ns_index = floor(para5)
    nturn1 = floor(para6)
    nturn11 = floor(para7)


    npara = 2*nolayer-1


    epsi=1.d-30
    tot_para = npara*ns   ! 450
    tot_ntc = ntc*ns      ! 140
    size0 = (ns-1)*(nolayer-1)  ! 796
    size1 = (ns-1)*npara

    eps = 1.d-10

    size2 = tot_ntc                ! ��Լ��

    ka=max(size2,tot_para)+1


    allocate(hz1_iter(ntc),hz1_iter_tmp(ntc),time(ntc))
    allocate(deltadobs(tot_ntc,1),deltam(tot_para,1))
    allocate(jacobi(tot_ntc,tot_para),jacobi1(ntc,npara))
    allocate(EYE(tot_para,tot_para),b0(size2,1))
    allocate(Btemp(tot_para,tot_para), invBtemp(tot_para,tot_para))
    allocate(Atemp(size2,tot_para),height(ns,nolayer-1),depth(ns,nolayer-1))
    allocate(m_pre(tot_para,1),m_next(tot_para,1),Vobs(tot_ntc,1),Vobs_tot(ns_real*ntc,1))
    allocate(deltaRp(size1,1),deltaRd(size0,1),Wp(size1,size1), Wd(size0,size0))
    allocate(Rp(size1,tot_para),Rd(size0,tot_para))
    allocate(s(ka))
    allocate(e(ka))
    allocate(work(ka))
    allocate(u(size2,size2))
    allocate(v(tot_para,tot_para))

    allocate(rho_iter(nolayer), hh_iter(nolayer))

    allocate(rho_pro(nolayer*ns_real), hh_pro(nolayer*ns_real))


    allocate(point1set(ns_real*2), point2set(ns_real*2), point3set(ns_real*2), point4set(ns_real*2),ht(ns_real))



    ! print*,"size0 = ",size0
    ! print*,"size1 = ",size1
    ! print*,"size2 = ",size2
    ! print*,"tot_ntc = ",tot_ntc
    ! print*,"tot_para = ",tot_para

    !******************************************************************
    Open (11, File='res2d.dat', Status='unknown')
    Open (12, File='err.dat', Status='unknown')
    Open (13, File='log.dat', Status='unknown')

    open (14, file='vobs_20ms.txt', status='old')

    Open (15, File='point1set.txt', Status='old')
    Open (16, File='point2set.txt', Status='old')
    Open (17, File='point3set.txt', Status='old')
    Open (18, File='point4set.txt', Status='old')

    Open (19, File='rho_pro_tunnel_20ms.txt', Status='old')
    Open (20, File='dep_pro_tunnel_20ms.txt', Status='old')


    ! Open (16, File='ht.txt', Status='old')

    Open (21, File='flag.dat', Status='unknown')
    Open (22, File='jac_rank0.dat', Status='unknown')
    Open (23, File='jac_rank1.dat', Status='unknown')

    do i=1,ns_real*nolayer
        read(19,*)rho_pro(i)
        read(20,*)hh_pro(i)
    end do


    ! do i=1,ns_real
    !     read(16,*)ht(i)
    ! end do




    do i=1,ns_real*2
        read(15,*)point1set(i)
        read(16,*)point2set(i)
        read(17,*)point3set(i)
        read(18,*)point4set(i)
    end do


    do i=1,ns_real*ntc
        read(14,*)Vobs_tot(i,1)
    end do

    Vobs_tot = Vobs_tot   ! ��mVת��ΪV


    Vobs(:,1) = Vobs_tot(1+(ns_index-1)*ntc:ns_index*ntc,1)


    !************************** initialize ***************************
    !��ʼ����ѡ���������ͬ�ľ��Ȱ�ռ�

    rho_iter = 10.d0

    hh_iter = 2.d0

    !************************** iteration begin ************************
    counter = 0
    lam = 10.0
    rerror_pre = 1.d8
    rerror = rerror_pre

    CALL CPU_TIME(cpu_start)

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc_current, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

101 Continue
    ! ��λ����
    !************************ ������������ ********************************
    EYE = 0.d0
    do i=1,tot_para
        EYE(i,i)=lam
    end do

    if(rerror .gt. rerror_pre)then
        lam = 2.0*lam
    else
        if(rerror .eq. rerror_pre)then
            lam = lam
        else
            lam = lam/2.0

        end if
    end if

    if(lam.le.1.d-30)then
        lam = 1.d-30
    end if

    rerror_pre = rerror
    counter = counter+1
    jacobi = 0.d0



    !*************************** �������ݽ�� *******************************
!    Write(,*)'calculate dB/dt in time domain ...'
    do k=1,ns

        htt = 0.01

        if(counter .eq. 1)then
            do i=1,nolayer
                rho_iter(i) = rho_pro(i+(ns_index-1)*nolayer)
                if(i .lt. nolayer)then
                    hh_iter(i) = hh_pro(i+(ns_index-1)*nolayer)
                end if
            end do

        end if


        if(counter .ge. 2)then
            !Write(*,*)'update parameter ...'
            do i=1,npara
                if(i .le. nolayer)then
                    rho_iter(i) = dexp(m_next(i+(k-1)*npara,1))
                else
                    hh_iter(i-nolayer) = dexp(m_next(i+(k-1)*npara,1))
                end if
            end do
        end if

        call forwardprocess(rho_iter,hh_iter,hz1_iter,nolayer,ntc,time,ns_index,&
        & point1set,point2set,point3set,point4set,htt,ns_real,t_st,t_ed,xr1, hr1,rt1, rr1,nturn1,nturn11)


        ! ��ǰʱ�̵Ĳ������絼��
        do i=1,npara
            if(i .le. nolayer)then
                m_pre(i+(k-1)*npara,1) = dlog(rho_iter(i))
            else

                m_pre(i+(k-1)*npara,1) = dlog(hh_iter(i-nolayer))

            end if
        end do


        ! ���㵱ǰ�������
        do i=1,ntc
            deltadobs(i+(k-1)*ntc,1) = Vobs(i+(k-1)*ntc,1)-hz1_iter(i)
        end do

        ! Write(13,*)hz1_iter
        ! Write(13,*)deltadobs


        
        ! �����ſɱȾ���
!        Write(6,*)'calculate Jacobi matrix ...'
        ! print*, 'this is it'
  call SYSTEM_CLOCK(iTimes2)
         call cal_jacobi(rank, nproc_current, 3,2 , rho_iter, hh_iter,jacobi1, nolayer, ntc, ns_index,&
         & point1set, point2set, point3set, point4set,htt,ns_real,t_st,t_ed,xr1, hr1,rt1, rr1,nturn1,nturn11)
  call SYSTEM_CLOCK(iTimes3)

        
        ! print*,'cur rank = ',rank
        ! print*,'cur jac size = ',size(jacobi1)
        
        ! if(rank==0)then
        !     write(23,*) jacobi1
        !     print*,'rank=',rank,', cur jac size = ',size(jacobi1)
        ! end if

        ! ! �ܵ��ſɱȾ���
        do i=1,ntc
            do j=1,npara
                jacobi(i+(k-1)*ntc,j+(k-1)*npara) = jacobi1(i,j)
            end do
        end do


        do i=1,nolayer-1
            height(k,i)=hh_iter(i)   ! �߶Ⱦ���
        end do

    end do

    if(rank==0)then
        Write(11,20)dexp(m_pre)
    end if

    do i=1,size2
        if(i .le. tot_ntc)then
            b0(i,1)=deltadobs(i,1)
            do j=1,tot_para
                Atemp(i,j)=jacobi(i,j)
            end do
        else
            if(i .le. tot_ntc+size1)then
                b0(i,1) = deltaRp(i-tot_ntc,1)
                do j=1,tot_para
                    Atemp(i,j)=Rp(i-tot_ntc,j)
                end do
            else
                b0(i,1) = deltaRd(i-tot_ntc-size1,1)
                do j=1,tot_para
                    Atemp(i,j)=Rd(i-tot_ntc-size1,j)
                end do
            end if
        endif
    end do


    call bmuav(Atemp,size2,tot_para,u,v,l,eps,ka,s,e,work)   !    

    Btemp = matmul(transpose(Atemp),Atemp)+EYE

    !****************** �����������Է��� *****************************
    call get_invMatrix(Btemp,tot_para,invBtemp)


    ! lambda ̫С��ʱ���ʹ������������
    deltam = matmul(matmul(matmul(matmul(transpose(v),invBtemp),&
        transpose(Atemp)),transpose(u)),b0)


    m_next = m_pre+deltam

    !************ calculate phid *****************

    temp0 = matmul(transpose(deltadobs),deltadobs)
    rerror = temp0(1,1)/tot_ntc
    if(rank==0)then
        Write(12, 100)rerror
    end if
    

    If(rerror>epsi .and. counter<3) Goto 101

    call MPI_Finalize(ierr)
    CALL CPU_TIME(cpu_finish)
    elapsed_time = cpu_finish - cpu_start
    WRITE(13,*) 'Elapsed CPU time:', elapsed_time, 'seconds'

    Write(21, *)1.0

  call SYSTEM_CLOCK(iTimes4)
    open (33, File='expmain4_record.txt', status='unknown')
   if(rank==0)then
    ! write(33,*) 'Wall Clock TIME T1andT4:', real(iTimes1)/real(rate), real(iTimes4)/real(rate)
    write(33,*) 'Wall Clock TIME T4-T1:', real(iTimes4-iTimes1)/real(rate)!80/1000!
    write(*,*) 'Wall Clock TIME T4-T1:', real(iTimes4-iTimes1)/real(rate)!80/1000!
    write(33,*) 'Wall Clock TIME T3-T2:', real(iTimes3-iTimes2)/real(rate)!80/1000!
    write(*,*) 'Wall Clock TIME T3-T2:', real(iTimes3-iTimes2)/real(rate)!80/1000!
  end if

10  Format(56E14.6)

20  Format(13E14.6)   ! ÿ�����Ĳ������� 2*nolayers-1

100 Format(E12.6)    ! Em.n��ʾ�ø�������������������ܹ�mλ��С�����nλ
    end program


    !!!!!!!!!!!!!!!!!!!!!----����ֵ�ֽ⡢��������Ӻ���----!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !																				  !
    !																				  !
    !   �Ӻ��������ԡ�Fortran�����㷨���򼯡����ڶ��棩��ʿ���������廪��ѧ������     !
    !																				  !
    !																				  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                          subroutine program 3 ����ֵ�ֽ�
    !!                           Householder�任��QR�㷨
    !// ForQuill v1.01 Beta www.fcode.cn
    Subroutine bmuav(a, m, n, u, v, l, eps, ka, s, e, work)
    Dimension a(m, n), u(m, m), v(n, n), s(ka), e(ka), work(ka)
    Real *8 a, u, v, s, e, d, work, dd, f, g, cs, sn, shh, sk, ek, b, c, sm, sm1, em1, eps
    it = 60
    k = n
    If (m-1<n) k = m - 1
    l = m
    If (n-2<m) l = n - 2
    If (l<0) l = 0
    ll = k
    If (l>k) ll = l
    If (ll>1) Then
        Do kk = 1, ll
            If (kk<=k) Then
                d = 0.0
                Do i = kk, m
                    d = d + a(i, kk)*a(i, kk)
                End Do
                s(kk) = sqrt(d)
                If (s(kk)/=0.0) Then
                    If (a(kk,kk)/=0.0) s(kk) = sign(s(kk), a(kk,kk))
                    Do i = kk, m
                        a(i, kk) = a(i, kk)/s(kk)
                    End Do
                    a(kk, kk) = 1.0 + a(kk, kk)
                End If
                s(kk) = -s(kk)
            End If
            If (n>=kk+1) Then
                Do j = kk + 1, n
                    If ((kk<=k) .And. (s(kk)/=0.0)) Then
                        d = 0.0
                        Do i = kk, m
                            d = d + a(i, kk)*a(i, j)
                        End Do
                        d = -d/a(kk, kk)
                        Do i = kk, m
                            a(i, j) = a(i, j) + d*a(i, kk)
                        End Do
                    End If
                    e(j) = a(kk, j)
                End Do
            End If
            If (kk<=k) Then
                Do i = kk, m
                    u(i, kk) = a(i, kk)
                End Do
            End If
            If (kk<=l) Then
                d = 0.0
                Do i = kk + 1, n
                    d = d + e(i)*e(i)
                End Do
                e(kk) = sqrt(d)
                If (e(kk)/=0.0) Then
                    If (e(kk+1)/=0.0) e(kk) = sign(e(kk), e(kk+1))
                    Do i = kk + 1, n
                        e(i) = e(i)/e(kk)
                    End Do
                    e(kk+1) = 1.0 + e(kk+1)
                End If
                e(kk) = -e(kk)
                If ((kk+1<=m) .And. (e(kk)/=0.0)) Then
                    Do i = kk + 1, m
                        work(i) = 0.0
                    End Do
                    Do j = kk + 1, n
                        Do i = kk + 1, m
                            work(i) = work(i) + e(j)*a(i, j)
                        End Do
                    End Do
                    Do j = kk + 1, n
                        Do i = kk + 1, m
                            a(i, j) = a(i, j) - work(i)*e(j)/e(kk+1)
                        End Do
                    End Do
                End If
                Do i = kk + 1, n
                    v(i, kk) = e(i)
                End Do
            End If
        End Do
    End If
    mm = n
    If (m+1<n) mm = m + 1
    If (k<n) s(k+1) = a(k+1, k+1)
    If (m<mm) s(mm) = 0.0
    If (l+1<mm) e(l+1) = a(l+1, mm)
    e(mm) = 0.0
    nn = m
    If (m>n) nn = n
    If (nn>=k+1) Then
        Do j = k + 1, nn
            Do i = 1, m
                u(i, j) = 0.0
            End Do
            u(j, j) = 1.0
        End Do
    End If
    If (k>=1) Then
        Do ll = 1, k
            kk = k - ll + 1
            If (s(kk)/=0.0) Then
                If (nn>=kk+1) Then
                    Do j = kk + 1, nn
                        d = 0.0
                        Do i = kk, m
                            d = d + u(i, kk)*u(i, j)/u(kk, kk)
                        End Do
                        d = -d
                        Do i = kk, m
                            u(i, j) = u(i, j) + d*u(i, kk)
                        End Do
                    End Do
                End If
                Do i = kk, m
                    u(i, kk) = -u(i, kk)
                End Do
                u(kk, kk) = 1.0 + u(kk, kk)
                If (kk-1>=1) Then
                    Do i = 1, kk - 1
                        u(i, kk) = 0.0
                    End Do
                End If
            Else
                Do i = 1, m
                    u(i, kk) = 0.0
                End Do
                u(kk, kk) = 1.0
            End If
        End Do
    End If
    Do ll = 1, n
        kk = n - ll + 1
        If ((kk<=l) .And. (e(kk)/=0.0)) Then
            Do j = kk + 1, n
                d = 0.0
                Do i = kk + 1, n
                    d = d + v(i, kk)*v(i, j)/v(kk+1, kk)
                End Do
                d = -d
                Do i = kk + 1, n
                    v(i, j) = v(i, j) + d*v(i, kk)
                End Do
            End Do
        End If
        Do i = 1, n
            v(i, kk) = 0.0
        End Do
        v(kk, kk) = 1.0
    End Do
    Do i = 1, m
        Do j = 1, n
            a(i, j) = 0.0
        End Do
    End Do
    m1 = mm
    it = 60
310 If (mm==0) Then
        l = 0
        If ((m>=n)) Then
            i = n
        Else
            i = m
        End If
        Do j = 1, i - 1
            a(j, j) = s(j)
            a(j, j+1) = e(j)
        End Do
        a(i, i) = s(i)
        If (m<n) a(i, i+1) = e(i)
        Do i = 1, n - 1
            Do j = i + 1, n
                Call real_exchange(v(i,j), v(j,i))
            End Do
        End Do
        Return
    End If
    !
    If (it==0) Then
        l = mm
        If (m>=n) Then
            i = n
        Else
            i = m
        End If
        Do j = 1, i - 1
            a(j, j) = s(j)
            a(j, j+1) = e(j)
        End Do
        a(i, i) = s(i)
        If (m<n) a(i, i+1) = e(i)
        Do i = 1, n - 1
            Do j = i + 1, n
                Call real_exchange(v(i,j), v(j,i))
            End Do
        End Do
        Return
    End If
    !
    kk = mm
320 kk = kk - 1
    If (kk/=0) Then
        d = abs(s(kk)) + abs(s(kk+1))
        dd = abs(e(kk))
        If (dd>eps*d) Goto 320
        e(kk) = 0.0
    End If
    If (kk==mm-1) Then
        kk = kk + 1
        If (s(kk)<0.0) Then
            s(kk) = -s(kk)
            Do i = 1, n
                v(i, kk) = -v(i, kk)
            End Do
        End If
335     If (kk/=m1) Then
            If (s(kk)<s(kk+1)) Then
                Call real_exchange(s(kk), s(kk+1))
                If (kk<n) Then
                    Do i = 1, n
                        Call real_exchange(v(i,kk), v(i,kk+1))
                    End Do
                End If
                If (kk<m) Then
                    Do i = 1, m
                        Call real_exchange(u(i,kk), u(i,kk+1))
                    End Do
                End If
                kk = kk + 1
                Goto 335
            End If
        End If
        it = 60
        mm = mm - 1
        Goto 310
    End If
    ks = mm + 1
360 ks = ks - 1
    If (ks>kk) Then
        d = 0.0
        If (ks/=mm) d = d + abs(e(ks))
        If (ks/=kk+1) d = d + abs(e(ks-1))
        dd = abs(s(ks))
        If (dd>eps*d) Goto 360
        s(ks) = 0.0
    End If
    If (ks==kk) Then
        kk = kk + 1
        d = abs(s(mm))
        If (abs(s(mm-1))>d) d = abs(s(mm-1))
        If (abs(e(mm-1))>d) d = abs(e(mm-1))
        If (abs(s(kk))>d) d = abs(s(kk))
        If (abs(e(kk))>d) d = abs(e(kk))
        sm = s(mm)/d
        sm1 = s(mm-1)/d
        em1 = e(mm-1)/d
        sk = s(kk)/d
        ek = e(kk)/d
        b = ((sm1+sm)*(sm1-sm)+em1*em1)/2.0
        c = sm*em1
        c = c*c
        shh = 0.0
        If ((b/=0.0) .Or. (c/=0.0)) Then
            shh = sqrt(b*b+c)
            If (b<0.0) shh = -shh
            shh = c/(b+shh)
        End If
        f = (sk+sm)*(sk-sm) - shh
        g = sk*ek
        Do i = kk, mm - 1
            Call sss(f, g, cs, sn)
            If (i/=kk) e(i-1) = f
            f = cs*s(i) + sn*e(i)
            e(i) = cs*e(i) - sn*s(i)
            g = sn*s(i+1)
            s(i+1) = cs*s(i+1)
            If ((cs/=1.0) .Or. (sn/=0.0)) Then
                Do j = 1, n
                    d = cs*v(j, i) + sn*v(j, i+1)
                    v(j, i+1) = -sn*v(j, i) + cs*v(j, i+1)
                    v(j, i) = d
                End Do
            End If
            Call sss(f, g, cs, sn)
            s(i) = f
            f = cs*e(i) + sn*s(i+1)
            s(i+1) = -sn*e(i) + cs*s(i+1)
            g = sn*e(i+1)
            e(i+1) = cs*e(i+1)
            If (i<m) Then
                If ((cs/=1.0) .Or. (sn/=0.0)) Then
                    Do j = 1, m
                        d = cs*u(j, i) + sn*u(j, i+1)
                        u(j, i+1) = -sn*u(j, i) + cs*u(j, i+1)
                        u(j, i) = d
                    End Do
                End If
            End If
        End Do
        e(mm-1) = f
        it = it - 1
        Goto 310
    End If
    If (ks==mm) Then
        kk = kk + 1
        f = e(mm-1)
        e(mm-1) = 0.0
        Do ll = kk, mm - 1
            i = mm + kk - ll - 1
            g = s(i)
            Call sss(g, f, cs, sn)
            s(i) = g
            If (i/=kk) Then
                f = -sn*e(i-1)
                e(i-1) = cs*e(i-1)
            End If
            If ((cs/=1.0) .Or. (sn/=0.0)) Then
                Do j = 1, n
                    d = cs*v(j, i) + sn*v(j, mm)
                    v(j, mm) = -sn*v(j, i) + cs*v(j, mm)
                    v(j, i) = d
                End Do
            End If
        End Do
        Goto 310
    End If
    kk = ks + 1
    f = e(kk-1)
    e(kk-1) = 0.0
    Do i = kk, mm
        g = s(i)
        Call sss(g, f, cs, sn)
        s(i) = g
        f = -sn*e(i)
        e(i) = cs*e(i)
        If ((cs/=1.0) .Or. (sn/=0.0)) Then
            Do j = 1, m
                d = cs*u(j, i) + sn*u(j, kk-1)
                u(j, kk-1) = -sn*u(j, i) + cs*u(j, kk-1)
                u(j, i) = d
            End Do
        End If
    End Do
    Goto 310
    End Subroutine bmuav

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!        subroutine program 4 subroutine program for ����ֵ�ֽ�
    Subroutine sss(f, g, cs, sn)
    Real *8 f, g, cs, sn, d, r
    If ((abs(f)+abs(g))==0.0) Then
        cs = 1.0
        sn = 0.0
        d = 0.0
    Else
        d = sqrt(f*f+g*g)
        If (abs(f)>abs(g)) d = sign(d, f)
        If (abs(g)>=abs(f)) d = sign(d, g)
        cs = f/d
        sn = g/d
    End If
    r = 1.0
    If (abs(f)>abs(g)) Then
        r = sn
    Else
        If (cs/=0.0) r = 1.0/cs
    End If
    f = d
    g = r
    Return
    End Subroutine sss

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  subroutine program 6 double precision exchange
    Subroutine real_exchange(t1, t2)
    Implicit None
    Double Precision t1, t2, cup
    cup = t1
    t1 = t2
    t2 = cup

    end
    subroutine get_invMatrix(A,n,invA)
    implicit none
    integer :: n
    real*8  :: A(n,n), invA(n,n)

    call solve(A, invA, n)

    end subroutine get_invMatrix


    subroutine solve(A, invA, n)
    implicit none
    integer :: i, n
    real*8  :: A(n,n), invA(n,n), E(n,n)

    E = 0.d0
    !����EΪ��λ����
    do i = 1, n
        E(i,i) = 1.d0
    end do

    call mateq(A, E, invA, n, n)

    end subroutine solve


    subroutine mateq(A, B, X, n, M)
    implicit none
    integer :: i, k, m, n, id_max
    real*8  :: A(n,n), B(n,M), X(n,M)

    real*8  :: Aup(n,n), elmax, temp
    !AbΪ������� [AB]
    real*8  :: AB(n,n+M), vtemp1(n+M), vtemp2(n+M), vtmp(n), xtmp(n)

    AB(1:n,1:n) = A

    AB(:,n+1:n+M) = B

    !##########################################################
    !  ���������Ԫ��ȥ���ĺ��Ĳ���
    do k = 1, n-1

        elmax  = abs(Ab(k,k))
        id_max = k

        !���Ϊ������Ԫ��
        !��γ������ҪĿ�Ĳ���Ϊ�˸�ֵ���Ԫ�ظ�elmax������Ϊ���ҳ����Ԫ�ض�Ӧ�ı��


        do i=k+1,n
            if (abs(Ab(i,k))>elmax) then
                elmax = Ab(i,k)

                id_max = i
            end if
        end do

        !���ˣ��Ѿ���ɲ������Ԫ�أ���������Ժ��� ��k�н���
        !��������Ԫ�أ���������
        vtemp1 = Ab(k,:)
        vtemp2 = Ab(id_max,:)


        Ab(k,:)      = vtemp2
        Ab(id_max,:) = vtemp1
        !
        !����һ�����Ϊ��������Ԫ�أ���������Ժ󼴰�����Ԫ������
        !#########################################################

        do i = k+1, n

            temp    = Ab(i,k)/Ab(k,k)

            Ab(i,:) = Ab(i,:)-temp*Ab(k,:)

        end do

    end do
    !-----------------------------
    ! ������һ����Ab�Ѿ���Ϊ������ʽ�ľ���
    !            | *  *  *  *  #  #|
    !     [A b]= | 0  *  *  *  #  #|
    !            | 0  0  *  *  #  #|
    !            | 0  0  0  *  #  #|
    !
    Aup(:,:) = AB(1:n,1:n)

    do i=1,m
        !�����������Ƿ�����Ļش�����
        vtmp = AB(:,n+i)
        call uptri(Aup,vtmp,xtmp,n)
        !�Ѽ�������ֵ��X
        X(:,i) = xtmp
    end do

    end subroutine mateq


    subroutine uptri(A,b,x,n)
    implicit none
    integer :: i, j, n
    real*8  :: A(n,n), b(n), x(n)

    x(n) = b(n)/A(n,n)
    !�ش�����
    do i = n-1, 1, -1

        x(i) = b(i)
        do j = i + 1,n
            x(i) = x(i)-a(i,j)*x(j)
        end do
        x(i) = x(i)/A(i,i)

    end do

    end subroutine uptri



    subroutine cal_jacobi(rank, nproc_current, Proc1, Proc2, m1,m2,jac,nlayer,ntc0,ns_id1,p11,p21,p31,&
    p41,htt1,ns_real1,t_st1,t_ed1,xr1,hr1,rt1,rr1,nturn1,nturn11)
    use mpi
    integer:: i,j,k,no,ntc0,ns_id1,ns_real1,nturn1,nturn11
    
    real *8 htt1
    Real *8 m1(nlayer), m2(nlayer), delta,jac_perRank(ntc0)
    real *8 jac(ntc0,2*nlayer-1) ! ע���ǰ����ŵ��ʵ��ܲ�����Ŀ
    Real *8 t_st1,t_ed1,xr1, hr1,rt1, rr1

    real *8,allocatable::time1(:)

    real *8 p11(ns_real1*2), p21(ns_real1*2), p31(ns_real1*2), p41(ns_real1*2)

    integer:: ierr, rank, nproc_current
    integer:: Proc1, Proc2 ! procedure num for work1/2
    integer:: workPProc1, workPProc2 ! work per procedure
    integer:: workrank1, workrank2
    Real *8,allocatable:: zeta(:)
    real :: start_time, end_time, total_time, compute_time, comm_time, sum_compute, sum_comm

    integer,allocatable:: npara_array(:)
    integer:: para_rank
    real*8, allocatable:: res1(:, :) ! 
    real*8, allocatable:: res2(:, :) ! 
    integer, dimension(2) :: starts, sizes, subsizes
    integer (kind=MPI_Address_kind) :: start, extent
    integer :: blocktype, resizedtype, doublePsize
    integer :: gridIdx, blockIdx, blockstride
    workPProc1 = (2*nlayer-1)/Proc1
    workPProc2 = 20/Proc2

    allocate(npara_array(workPProc1))
    allocate(zeta(workPProc2))
    allocate(time1(ntc0))
    allocate(res1(ntc0, workPProc1)) ! 
    allocate(res2(ntc0, workPProc1 * Proc1 * Proc2)) ! 

    delta = 2.d-06  ! 
    
    workrank1 = rank / Proc2
    workrank2 = mod(rank , Proc2)

    ! describe what these subblocks look like inside the full concatenated array
    ! sizes    = [ ntc0, (2*nlayer-1)*workPProc2 ]
    ! subsizes = [ ntc0, (2*nlayer-1)*workPProc2 / workPProc1 ]
    sizes    = [ ntc0, workPProc1 * Proc1 * Proc2 ]
    subsizes = [ ntc0, workPProc1 ]

    starts   = [ 0, 0 ]
    call MPI_Type_create_subarray( 2, sizes, subsizes, starts,     &
                                   MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                   blocktype, ierr)

    start = 0
    call MPI_Type_size(MPI_DOUBLE_PRECISION, doublePsize, ierr)
    extent = doublePsize * ntc0 * workPProc1
    ! print*, extent
    call MPI_Type_create_resized(blocktype, start, extent, resizedtype, ierr)
    call MPI_Type_commit(resizedtype, ierr)


    do i = 1,workPProc1
        npara_array(i) = workrank1 * workPProc1 + i
    end do

    do i = 1,workPProc2
        zeta(i) = workrank2 * workPProc2 + i
    end do

    do j=1,size(npara_array)
        para_rank = npara_array(j)
        jac_perRank = 0.d0

        call cal_diff(m1,m2,zeta,workPProc2,delta,para_rank,jac_perRank,nlayer,ntc0,ns_id1,p11,p21,p31,p41,htt1,ns_real1,&
        & t_st1,t_ed1,xr1,hr1,rt1,rr1,nturn1,nturn11)
        
        call MPI_barrier(MPI_COMM_WORLD, ierr)


        do i=1,ntc0
            res1(i,j)=jac_perRank(i)
        end do
    end do

    call MPI_Gather( res1, ntc0*workPProc1, MPI_DOUBLE_PRECISION, &  ! everyone send 3*2 ints
                    res2, 1, resizedtype,                   &  ! root gets 1 resized type from everyone
                    0, MPI_COMM_WORLD, ierr)

    ! call MPI_barrier(MPI_COMM_WORLD, ierr)
    if(rank==0)then
        jac = 0.d0
        do i = 1,ntc0
            do j = 1, Proc1*workPProc1
                blockIdx = (j-1) / workPProc1
                ! blockstride = workPProc1
                blockstride = workPProc1 * Proc2
                do k = 1,Proc2 !1,2
                    gridIdx = (k-1)*workPProc1
                    jac(i,j) = jac(i,j) + res2(i, mod(j-1,workPProc1)+1 + gridIdx + blockIdx * blockstride)
                end do
            end do
        end do
        write(22,*) "jac from rank0:"
        write(22,*) jac
    end if

    call MPI_Bcast(jac, ntc0*(Proc1*workPProc1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if(rank==1)then
        write(23,*) "jac from rank1:"
        write(23,*) jac
    end if
    end subroutine cal_jacobi

    subroutine cal_diff(rhoo, hhh, zeta1,workPProc2, del, jth, temp, nlayer,ntc1,ns_id,p1,p2,p3,p4,htt2,ns_real2,&
        & t_st2,t_ed2,xr1, hr1,rt1, rr1,nturn1,nturn11)

    integer:: workPProc2
    integer:: i,k,jth,ntc1, ns_id,ns_real2,nturn1,nturn11
    Real *8 zeta1(workPProc2), del, rhoo(nlayer), hhh(nlayer), temp(ntc1),htt2
    Real *8 rhonew1(nlayer), hhnew1(nlayer), rhonew2(nlayer), hhnew2(nlayer)

    Real *8,allocatable:: hzz1(:),hzz11(:), time0(:)

    real *8 p1(ns_real2*2), p2(ns_real2*2), p3(ns_real2*2), p4(ns_real2*2)
    Real *8 t_st2,t_ed2,xr1, hr1,rt1, rr1

    allocate(hzz1(ntc1))
    allocate(hzz11(ntc1), time0(ntc1))

    do i=1,nlayer
        rhonew1(i) = rhoo(i)
        hhnew1(i) = hhh(i)

        rhonew2(i) = rhoo(i)
        hhnew2(i) = hhh(i)

    end do

    temp=0.d0
    do i=1,size(zeta1)  ! different zeta
        if(jth .le. nlayer)then
            rhonew1(jth) = dexp(dlog(rhoo(jth))+zeta1(i)*del)
            rhonew2(jth) = dexp(dlog(rhoo(jth))-zeta1(i)*del)
        else
            hhnew1(jth-nlayer) = dexp(dlog(hhh(jth-nlayer))+zeta1(i)*del)
            hhnew2(jth-nlayer) = dexp(dlog(hhh(jth-nlayer))-zeta1(i)*del)
        endif

        call forwardprocess(rhonew1, hhnew1, hzz1, nlayer,ntc1,time0,ns_id,p1,p2,p3,p4,htt2,&
        &ns_real2,t_st2,t_ed2,xr1, hr1,rt1, rr1,nturn1,nturn11)
        call forwardprocess(rhonew2, hhnew2, hzz11, nlayer,ntc1,time0,ns_id,p1,p2,p3,p4,htt2,&
        & ns_real2,t_st2,t_ed2,xr1, hr1,rt1, rr1,nturn1,nturn11)

        do k=1,ntc1
            temp(k)=temp(k)+(hzz1(k)-hzz11(k))*1.d0/2/20/zeta1(i)/del
        end do

    end do

    end subroutine cal_diff



    Subroutine forwardprocess(rho, hh, hz1, nlayer,nt,t,ns_id,p1,p2,p3,p4,ht,ns1,t_st,t_ed,xr, hr,rt, rr,nturn,nturn1)
    integer::nlayer, nt, npls, ns_id, ns1
    Real *8 t_st,t_ed, delta_t
    Real *8 t(nt), hz1(nt),tlog(nt)

    Real *8 s2, ss2
    Real *8 r, rplus, tt1, tt2, tt3, tt4, pi
    Real *8 rho(nlayer), hh(nlayer), frq(67)

    Real *8 xr, yr, ht, hr, zplus, zminus

    real *8 slope1, slope2, slope3

    real *8 p1(ns1*2), p2(ns1*2), p3(ns1*2), p4(ns1*2)

    real *8  rt, rr
    integer::nturn, nturn1


    Complex *16 func(5, 67)

    Data pi/3.1415926D0/
    Data ngau/10000/

    Common /para/r

    Common /funn/frq, func

    !**************************************************************
    Call filter
    !**************************************************************

    !hr = 0.01

    !xr = 0.58
    yr = 0.

    zplus = ht - hr
    zminus = ht + hr
    r = dsqrt(xr*xr+yr*yr)
    rplus = dsqrt(r*r+zplus*zplus)

    !**************************************************************
    ! ������Ȧ
    !rt = 0.5             ! ��Ȧ�İ뾶
    !nturn = 3.           ! ��Ȧ������

    !rr = 0.25            ! �뾶
    !nturn1 = 20          ! ������Ȧ������

!    t_st = 2.024e-3
!    t_ed = 20e-3

    npls = 1

    delta_t = (dlog10(t_ed)-dlog10(t_st))/(nt-1)


    do i=1,nt
        tlog(i)=dlog10(t_st)+(i-1)*delta_t
        t(i)=10**tlog(i)
    end do


    ic = 3

    If (ic==2) Then   ! �����Է���
        tt1 = 2e-3
        tt2 = t_ed-tt1
    End If

    If (ic==3) Then   ! �����Է���

        tt1 = p2(1+(ns_id-1)*2)-p1(1+(ns_id-1)*2)
        tt2 = p3(1+(ns_id-1)*2)-p2(1+(ns_id-1)*2)
        tt3 = p4(1+(ns_id-1)*2)-p3(1+(ns_id-1)*2)
        tt4 = t_ed-tt1-tt2-tt3

        slope1 = (p2(2+(ns_id-1)*2)-p1(2+(ns_id-1)*2))/tt1
        slope2 = (p2(2+(ns_id-1)*2)-p3(2+(ns_id-1)*2))/tt2
        slope3 = (p3(2+(ns_id-1)*2)-p4(2+(ns_id-1)*2))/tt3

    End If


    !********************* impulse and step wave ************************
    If (ic==0 .Or. ic==1) Then
        ik = 0
        Do i = 1,nt

            Call frt(rho, hh, t(i), hz1(i), 2, zplus, zminus, ic, ik, nlayer)

            ik = ik + 1

            ! Transformation of cylindrical coordinate system
            hz1(i) = (4*rt*rt*nturn)*hz1(i)*(pi*rr**2*nturn1)

        End Do

    Else If (ic==3) Then
        Do i = 1, nt

            ss2 = 0.D0

            ik = 0

            Do ip = 1, npls
                If (t(i)>(ip-1)*(tt1+tt2+tt3+tt4) .And. t(i)<ip*(tt1+tt2+tt3+tt4)) kpls = ip
            End Do

            Do ip = 1, kpls - 1

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4), s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + s2*slope1

                ik = ik + 1

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4)-tt1, s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + (-1)*s2*(slope1+slope2)

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4)-tt1-tt2, s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + (-1)*s2*(slope3-slope2)

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4)-tt1-tt2-tt3, s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + s2*slope3
            End Do

            !****************** contribution from resting pulses**************************
            If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4))Then

                Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4), s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + s2*slope1

                ik = ik + 1

                If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4)+tt1)Then

                    Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4)-tt1, s2, 2, zplus, zminus, 1, ik, nlayer)

                    ss2 = ss2 + (-1)*s2*(slope1+slope2)


                    If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4)+tt1+tt2)Then

                        Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4)-tt1-tt2, s2, 2, zplus, zminus, 1, ik, nlayer)

                        ss2 = ss2 + (-1)*s2*(slope3-slope2)

                        If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4)+tt1+tt2+tt3)Then

                            Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4)-tt1-tt2-tt3, s2, 2, zplus, zminus, 1, ik, nlayer)

                            ss2 = ss2 + s2*slope3

                        End If
                    End If
                End If
            End If

            ! dB/dt for trapezoid
            hz1(i) = 4*rt*rt*nturn*ss2*(pi*rr**2*nturn1)
        End Do

    End If


    End Subroutine forwardprocess





    Subroutine gauleg(x1, x2, x, w, n)
    Implicit Real *8(A-H, O-Z)
    Real *8 x1, x2, xm, xl, x(n), w(n)
    Parameter (eps=3.D-14)
    m = (n+1)/2
    xm = 0.5D0*(x2+x1)
    xl = 0.5D0*(x2-x1)
    Do i = 1, m
        z = cos(3.141592654D0*(i-0.25D0)/(n+0.5D0))
1       Continue
        p1 = 1.D0
        p2 = 0.D0
        Do j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
        End Do
        pp = n*(z*p1-p2)/(z*z-1.D0)
        z1 = z
        z = z1 - p1/pp
        If (abs(z-z1)>eps) Goto 1
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = 2.D0*xl/((1.D0-z*z)*pp*pp)
        w(n+1-i) = w(i)
    End Do
    Return
    End Subroutine gauleg

    !        spl(nfrq, nc, frq, funr0, f, funr1)
    Subroutine spl(nx, n2, x, fx, x2, fx2)
    Real *8 x(nx), fx(nx), c(3, nx), x2(n2), fx2(n2), xint, xl1, xl2
    Call splin1(nx, fx, c)
    xl1 = dlog10(x(1))
    xl2 = dlog10(x(nx))
    Do ix = 1, n2
        xint = dlog10(x2(ix))
        Call splin2(nx, xint, xl1, xl2, c, fx2(ix))
    End Do
    End Subroutine spl

    Subroutine splin1(n, y, c)
    Real *8 y(n), c(3, n), p
    n1 = n - 1
    Do i = 2, n1
        c(1, i) = y(i+1) - 2.*y(i) + y(i-1)
    End Do
    c(2, 1) = 0.D0
    c(3, 1) = 0.D0
    Do i = 2, n1
        p = 4. + c(2, i-1)
        c(2, i) = -1.D0/p
        c(3, i) = (c(1,i)-c(3,i-1))/p
    End Do
    c(1, n) = 0.D0
    Do ii = 2, n1
        i = n + 1 - ii
        c(1, i) = c(2, i)*c(1, i+1) + c(3, i)
    End Do
    c(1, 1) = 0.
    Do i = 1, n1
        c(2, i) = y(i+1) - y(i) - c(1, i+1) + c(1, i)
        c(3, i) = y(i) - c(1, i)
    End Do
    c(3, n) = y(n)
    Return
    End Subroutine splin1

    Subroutine splin2(n, xint, x1, x2, c, yint)
    Real *8 c(3, n), xint, x1, x2, yint, h, u, p, q
    h = (x2-x1)/dble(float(n-1))
    If (xint<x1) Goto 10
    If (xint>=x2) Goto 20
    u = (xint-x1)/h
    i = 1 + int(u)
    p = u - i + 1
    q = 1.D0 - p
    yint = c(1, i)*q**3 + c(1, i+1)*p**3 + c(2, i)*p + c(3, i)
    Return
10  p = (xint-x1)/h
    yint = c(2, 1)*p + c(3, 1)
    Return
20  p = (xint-x2)/h
    yint = c(2, n-1)*p + c(3, n)
    Return
    End Subroutine splin2



    Subroutine frt(rho1, hh1, t, ft, item, zplus, zminus, ic, ik, nlayer)
    Complex *16 fun, iomega, func(5, 67)
    Real *8 t, ft, zplus, zminus, pi, q, rho1(nlayer), hh1(nlayer)
    Real *8 frq(67), funr0(67), funi0(67)
    Real *8 f(160), omega(160), funr1(160), funi1(160), h(200)
    Common /funn/frq, func
    Data pi, q/3.141592654D0, 1.258925412D0/
    Data ncnull, nc, ndec, (h(i), i=1, 160)/80,160,10,&
        & 2.59511139938829d-13,3.66568771323555d-13,5.17792876616242d-13,&
        & 7.31400730405791d-13,1.03313281156235d-12,1.45933600088387d-12,&
        & 2.06137146234699d-12,2.91175733962418d-12,4.11297804457870d-12,&
        & 5.80971771117984d-12,8.20647323099742d-12,1.15919058389365d-11,&
        & 1.63740746547780d-11,2.31288803930431d-11,3.26705938902288d-11,&
        & 4.61481520721098d-11,6.51864545047052d-11,9.20775899532545d-11,&
        & 1.30064200980219d-10,1.83718747396255d-10,2.59512512377884d-10,&
        & 3.66566596154242d-10,5.17796324027279d-10,7.31395266627501d-10,&
        & 1.03314147106736d-09,1.45932227649333d-09,2.06139321404013d-09,&
        & 2.91172286551380d-09,4.11303268236158d-09,5.80963111612975d-09,&
        & 8.20661047490285d-09,1.15916883220051d-08,1.63744193958818d-08,&
        & 2.31283340152144d-08,3.26714598407299d-08,4.61467796330556d-08,&
        & 6.84744728867720d-08,5.46574677490374d-08,1.13319898777493d-07,&
        & 2.16529974157527d-07,2.88629942214140d-07,3.42872728051125d-07,&
        & 4.79119488706262d-07,7.42089418889752d-07,1.07736520535271d-06,&
        & 1.46383231306575d-06,2.01727682134668d-06,2.89058197617431d-06,&
        & 4.15237808867022d-06,5.84448989361742d-06,8.18029430348419d-06,&
        & 1.15420854481494d-05,1.63897017145322d-05,2.31769096113890d-05,&
        & 3.26872676331330d-05,4.60786866701851d-05,6.51827321351636d-05,&
        & 9.20862589540037d-05,1.30169142615951d-04,1.83587481111627d-04,&
        & 2.59595544393723d-04,3.66324383719323d-04,5.18210697462501d-04,&
        & 7.30729969562531d-04,1.03385239132389d-03,1.45738764044730d-03,&
        & 2.06298256402732d-03,2.90606401578959d-03,4.11467957883740d-03,&
        & 5.79034253321120d-03,8.20005721235220d-03,1.15193892333104d-02,&
        & 1.63039398900789d-02,2.28256810984487d-02,3.22248555163692d-02,&
        & 4.47865101670011d-02,6.27330674874545d-02,8.57058672847471d-02,&
        & 1.17418179407605d-01,1.53632645832305d-01,1.97718111895102d-01,&
        & 2.28849924263247d-01,2.40310905012422d-01,1.65409071929404d-01,&
        & 2.84709685167114d-03,-2.88015846269687d-01,-3.69097391853225d-01,&
        & -2.50109865922601d-02,5.71811109500426d-01,-3.92261390212769d-01,&
        & 7.63282774297327d-02,5.16233692927851d-02,-6.48015160576432d-02,&
        & 4.89045522502552d-02,-3.26934307794750d-02,2.10542570949745d-02,&
        & -1.33862848934736d-02,8.47098801479259d-03,-5.35134515919751d-03,&
        & 3.37814023806349d-03,-2.13157364002470d-03,1.34506352474558d-03,&
        & -8.48929743771803d-04,5.35521822356713d-04,-3.37744799986382d-04,&
        & 2.13268792633204d-04,-1.34629969723156d-04,8.47737416679279d-05,&
        & -5.34940635827096d-05,3.3904416298191d-05,-2.13315638358794d-05,&
        & 1.33440911625019d-05,-8.51629073825634d-06,5.44362672273211d-06,&
        & -3.32112278417896d-06,2.07147190852386d-06,-1.42009412555511d-06,&
        & 8.78247754998004d-07,-4.5566280473703d-07,3.38598103040009d-07,&
        & -2.87407830772251d-07,1.07866150545699d-07,-2.4724024185358d-08,&
        & 5.35535110396030d-08,-3.3789981131378d-08,2.13200367531820d-08,&
        & -1.34520337740075d-08,8.48765950790546d-09,-5.35535110396018d-09,&
        & 3.37899811131383d-09,-2.13200367531819d-09,1.34520337740075d-09,&
        & -8.48765950790576d-10,5.35535110396015d-10,-3.37899811131382d-10,&
        & 2.13200367531811d-10,-1.34520337740079d-10,8.48765950790572d-11,&
        & -5.35535110396034d-11,3.37899811131381d-11,-2.13200367531818d-11,&
        & 1.34520337740074d-11,-8.48765950790571d-12,5.35535110396031d-12,&
        & -3.37899811131379d-12,2.13200367531817d-12,-1.34520337740073d-12,&
        & 8.48765950790567d-13,-5.35535110396029d-13,3.37899811131377d-13,&
        & -2.13200367531816d-13,1.34520337740078d-13,-8.48765950790596d-14,&
        & 5.35535110396007d-14,-3.37899811131377d-14,2.13200367531816d-14,&
        & -1.34520337740083d-14,8.4876550790558d-15,-5.35535110396025d-15,&
        & 3.37899811131389d-15/
    Data nfrq, (frq(i), i=1, 67)/67,&
        & 0.10000000d-02,0.14677993d-02,0.21544347d-02,0.31622777d-02,&
        & 0.46415888d-02,0.68129207d-02,0.10000000d-01,0.14677993d-01,&
        & 0.21544347d-01,0.31622777d-01,0.46415888d-01,0.68129207d-01,&
        & 0.10000000d+00,0.14677993d+00,0.21544347d+00,0.31622777d+00,&
        & 0.46415888d+00,0.68129207d+00,0.10000000d+01,0.14677993d+01,&
        & 0.21544347d+01,0.31622777d+01,0.46415888d+01,0.68129207d+01,&
        & 0.10000000d+02,0.14677993d+02,0.21544347d+02,0.31622777d+02,&
        & 0.46415888d+02,0.68129207d+02,0.10000000d+03,0.14677993d+03,&
        & 0.21544347d+03,0.31622777d+03,0.46415888d+03,0.68129207d+03,&
        & 0.10000000d+04,0.14677993d+04,0.21544347d+04,0.31622777d+04,&
        & 0.46415888d+04,0.68129207d+04,0.10000000d+05,0.14677993d+05,&
        & 0.21544347d+05,0.31622777d+05,0.46415888d+05,0.68129207d+05,&
        & 0.10000000d+06,0.14677993d+06,0.21544347d+06,0.31622777d+06,&
        & 0.46415888d+06,0.68129207d+06,0.10000000d+07,0.14677993d+07,&
        & 0.21544347d+07,0.31622777d+07,0.46415888d+07,0.68129207d+07,&
        & 0.10000000d+08,0.14677993d+08,0.21544347d+08,0.31622777d+08,&
        & 0.46415888d+08,0.68129207d+08,0.1000000d+09/
    If (ik==0) Then
        Do i = 1, nfrq
            Call forward(rho1, hh1, frq(i), func(item,i), item, zplus, zminus, nlayer) ! mu0*H = B
        End Do
    End If
    Do nn = 1, nc
        n = -nc + ncnull + nn
        omega(nn) = q**(-(n-1))/t
        f(nn) = omega(nn)/(2.D0*pi)
    End Do

    Do i = 1, nfrq
        funr0(i) = dreal(func(item,i))
        funi0(i) = dimag(func(item,i))
    End Do
    Call spl(nfrq, nc, frq, funr0, f, funr1)  !interpolation
    Call spl(nfrq, nc, frq, funi0, f, funi1)
    ft = 0.D0
    Do nn = 1, nc

        If (ic==0) Then
            iomega = (1.D0, 0.D0)
        Else If (ic==1) Then
            iomega = 1.D0/((0.,-1.D0)*omega(nn)) ! -1/iw
        End If

        fun = dcmplx(funr1(nn), funi1(nn))*iomega
        ita = max0(1, nn-nc+1)
        ite = min0(1, nn)

        Do it = ita, ite
            itn = nc - nn + it
            ft = ft + dimag(fun)*dsqrt(omega(nn))*h(itn) ! primary field st
        End Do

    End Do
    ft = -ft*dsqrt(2.D0/pi/t)
    Return
    End Subroutine frt

    Subroutine forward(rho2, hh2,  f, fun, item, zplus, zminus, nlayer)
    Complex *16 t3, t5, t6, hf, fun
    Real *8 pi, f, zplus, zminus
    Real *8 r, rho2(nlayer), hh2(nlayer)
    Common /para/r
    pi = 3.1415926D0
    If (item==1)Then
        hf = t6(rho2, hh2,f, zminus, nlayer)/(4.D0*pi)
    Else If (item==2)Then
        hf = -t3(rho2, hh2,f, zminus, nlayer)/(4.D0*pi)
    Else If (item==3)Then
        hf = (-t3(rho2, hh2,f,zminus, nlayer)+t5(rho2, hh2,f,zminus, nlayer)/r)/(4.D0*pi)
        If (zplus<0.D0) hf = -hf
    Else If (item==4)Then
        hf = t5(rho2, hh2,f, zminus, nlayer)/(4.D0*pi*r)
        If (zplus<0.D0) hf = -hf
    Else If (item==5) Then
        hf = -t6(rho2, hh2,f, zminus, nlayer)/(4.D0*pi)
        If (zplus<0.D0) hf = -hf
    Else
        Print *, 'item must be between 1 and 5'
    End If
    fun = hf*4.D-7*pi
    Return
    End Subroutine forward


    Complex *16 Function t3(rho3, hh3, f, z, nlayer)
    Complex *16 s, s1, b
    Real *8 r, rho3(nlayer), hh3(nlayer)
    Real *8 h0, h1, f, z, u, fac, expc
    Common /para/r
    Common /hankel/nc, ncnull, h0(100), h1(100)
    fac = 0.1D0*dlog(10.D0)
    s = (0.D0, 0.D0)
    Do nn = 1, nc
        nu = nn
        mn = nc - nn + 1
        nnn = ncnull - nc + nu
        u = expc(-(nnn-1)*fac)/r
        s1 = (b(rho3, hh3,f,u, nlayer)-u)/(b(rho3, hh3,f,u, nlayer)+u)*expc(-u*z)*u*u
        s = s + s1*h0(mn)
    End Do
    t3 = s/r
    Return
    End Function t3

    Complex *16 Function t5(rho5, hh5, f, z, nlayer)
    Complex *16 s, s1, b
    Real *8 r, rho5(nlayer), hh5(nlayer)
    Real *8 h0, h1, f, z, u, fac, expc
    Common /para/r
    Common /hankel/nc, ncnull, h0(100), h1(100)
    fac = 0.1D0*dlog(10.D0)
    s = (0.D0, 0.D0)
    Do nn = 1, nc
        nu = nn
        mn = nc - nn + 1
        nnn = ncnull - nc + nu
        u = expc(-(nnn-1)*fac)/r
        s1 = (b(rho5, hh5,f,u, nlayer)-u)/(b(rho5, hh5,f,u, nlayer)+u)*expc(-u*z)*u
        s = s + s1*h1(mn)
    End Do
    t5 = s/r
    Return
    End Function t5

    Complex *16 Function t6(rho6, hh6, f, z, nlayer)
    Complex *16 s, s1, b
    Real *8 r, rho6(nlayer), hh6(nlayer)
    Real *8 h0, h1, f, z, u, fac, expc
    Common /para/r
    Common /hankel/nc, ncnull, h0(100), h1(100)
    fac = 0.1D0*dlog(10.D0)
    s = (0.D0, 0.D0)
    Do nn = 1, nc
        nu = nn
        mn = nc - nn + 1
        nnn = ncnull - nc + nu
        u = expc(-(nnn-1)*fac)/r
        s1 = (b(rho6, hh6, f,u, nlayer)-u)/(b(rho6, hh6, f,u, nlayer)+u)*expc(-u*z)*u*u
        s = s + s1*h1(mn)
    End Do
    t6 = s/r
    Return
    End Function t6

    ! output: BE
    ! input: f - frequency
    ! input: u - k
    ! input: mu - �ŵ���
    ! nlayer -  no. of layers
    Complex *16 Function b(rho4, hh4, f, u, nlayer)
    Complex *16 alpha, s1, s2
    Real *8 f, u, pi
    Real *8 r, rho4(nlayer), hh4(nlayer)
    real *8 mu0
    Common /para/r

    pi = 3.1415926D0
    mu0 = pi*4d-07

    b = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/rho4(nlayer))  ! alpha nlayer omega = 2*pi*f
    If (nlayer==1) Return
    Do i = nlayer-1, 1, -1

        alpha = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/rho4(i))
        s1 = (0.D0, 0.D0)
        If (dreal(2.D0*alpha*hh4(i))<400.D0) s1 = cdexp(-2.D0*alpha*hh4(i))
        s2 = (1.D0-s1)/(1.D0+s1)  ! tanh
        b = alpha*(b+alpha*s2)/(alpha+b*s2)
    End Do

    End Function b

    Real *8 Function expc(x)
    Real *8 x, x1
    x1 = x
    If (dabs(x1)>650.D0) x1 = dsign(650.D0, x1)
    expc = dexp(x1)
    Return
    End Function expc

    Subroutine filter
    Common /hankel/nc, ncnull, h0(100), h1(100)
    Real *8 h0, h1

    Data(h0(i),i=1,48)/&
        & 2.89878288d-07,3.64935144d-07,4.59426126d-07,5.78383226d-07,&
        & 7.28141338d-07,9.16675639d-07,1.15402625d-06,1.45283298d-06,&
        & 1.82900834d-06,2.30258511d-06,2.89878286d-06,3.64935148d-06,&
        & 4.59426119d-06,5.78383236d-06,7.28141322d-06,9.16675664d-06,&
        & 1.15402621d-05,1.45283305d-05,1.82900824d-05,2.30258527d-05,&
        & 2.89878259d-05,3.64935186d-05,4.59426051d-05,5.78383329d-05,&
        & 7.28141144d-05,9.16675882d-05,1.15402573d-04,1.45283354d-04,&
        & 1.82900694d-04,2.30258630d-04,2.89877891d-04,3.64935362d-04,&
        & 4.59424960d-04,5.78383437d-04,7.28137738d-04,9.16674828d-04,&
        & 1.15401453d-03,1.45282561d-03,1.82896826d-03,2.30254535d-03,&
        & 2.89863979d-03,3.64916703d-03,4.59373308d-03,5.78303238d-03,&
        & 7.27941497d-03,9.16340705d-03,1.15325691d-02,1.45145832d-02/

    Data(h0(i),i=49,100)/&
        & 1.82601199d-02,2.29701042d-02,2.88702619d-02,3.62691810d-02,&
        & 4.54794031d-02,5.69408192d-02,7.09873072d-02,8.80995426d-02,&
        & 1.08223889d-01,1.31250483d-01,1.55055715d-01,1.76371506d-01,&
        & 1.85627738d-01,1.69778044d-01,1.03405245d-01,-3.02583233d-02,&
        & -2.27574393d-01,-3.62173217d-01,-2.05500446d-01,3.37394873d-01,&
        & 3.17689897d-01,-5.13762160d-01,3.09130264d-01,-1.26757592d-01,&
        & 4.61967890d-02,-1.80968674d-02,8.35426050d-03,-4.47368304d-03,&
        & 2.61974783d-03,-1.60171357d-03,9.97717882d-04,-6.26275815d-04,&
        & 3.94338818d-04,-2.48606354d-04,1.56808604d-04,-9.89266288d-05,&
        & 6.24152398d-05,-3.93805393d-05,2.48472358d-05,-1.56774945d-05,&
        & 9.89181741d-06,-6.24131160d-06,3.93800058d-06,-2.48471018d-06,&
        & 1.56774609d-06,-9.89180896d-07,6.24130948d-07,-3.93800005d-07,&
        & 2.48471005d-07,-1.56774605d-07,9.89180888d-08,-6.24130946d-08/

    Data(h1(i),i=1,48)/&
        & 1.84909557d-13,2.85321327d-13,4.64471808d-13,7.16694771d-13,&
        & 1.16670043d-12,1.80025587d-12,2.93061898d-12,4.52203829d-12,&
        & 7.36138206d-12,1.13588466d-11,1.84909557d-11,2.85321327d-11,&
        & 4.64471808d-11,7.16694771d-11,1.16670043d-10,1.80025587d-10,&
        & 2.93061898d-10,4.52203829d-10,7.36138206d-10,1.13588466d-09,&
        & 1.84909557d-09,2.85321326d-09,4.64471806d-09,7.16694765d-09,&
        & 1.16670042d-08,1.80025583d-08,2.93061889d-08,4.52203807d-08,&
        & 7.36138149d-08,1.13588452d-07,1.84909521d-07,2.85321237d-07,&
        & 4.64471580d-07,7.16694198d-07,1.16669899d-06,1.80025226d-06,&
        & 2.93060990d-06,4.52201549d-06,7.36132477d-06,1.13587027d-05,&
        & 1.84905942d-05,2.85312247d-05,4.64449000d-05,7.16637480d-05,&
        & 1.16655653d-04,1.79989440d-04,2.92971106d-04,4.51975783d-04/

    Data(h1(i),i=49,100)/&
        & 7.35565435d-04,1.13444615d-03,1.84548306d-03,2.84414257d-03,&
        & 4.62194743d-03,7.10980590d-03,1.15236911d-02,1.76434485d-02,&
        & 2.84076233d-02,4.29770596d-02,6.80332569d-02,9.97845929d-02,&
        & 1.51070544d-01,2.03540581d-01,2.71235377d-01,2.76073871d-01,&
        & 2.16691977d-01,-7.83723737d-02,-3.40675627d-01,-3.60693673d-01,&
        & 5.13024526d-01,-5.94724729d-02,-1.95117123d-01,1.99235600d-01,&
        & -1.38521553d-01,8.79320859d-02,-5.50697146d-02,3.45637848d-02,&
        & -2.17527180d-02,1.37100291d-02,-8.64656417d-03,5.45462758d-03,&
        & -3.44138864d-03,2.17130686d-03,-1.36998628d-03,8.64398952d-04,&
        & -5.45397874d-04,3.44122545d-04,-2.17126585d-04,1.36997597d-04,&
        & -8.64396364d-05,5.45397224d-05,-3.44122382d-05,2.17126544d-05,&
        & -1.36997587d-05,8.64396338d-06,-5.45397218d-06,3.44122380d-06,&
        & -2.17126543d-06,1.36997587d-06,-8.64396337d-07,5.45397218d-07/

    nc = 100
    ncnull = 60
    Return
    End Subroutine filter
