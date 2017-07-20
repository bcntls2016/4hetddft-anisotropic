!------------------------------------------------------------------
!---                    Subroutine printoutc                    ---
!------------------------------------------------------------------

Subroutine Printoutc(pr)

Use Impressio

Type(info_printout) pr

  open(10,file=pr%namefile,form='FORMATTED')
    write(10,'("#  ")')
    write(10,'("#  Density after ",I15," iterations")')pr%it
    write(10,'("#  Actual deltatps............:",1p,E20.10)')pr%dtps
    Write(10,'("#  Total evolution time(ps)...:",1p,E20.12)')pr%time
!
    Write(10,'("#  Sigma potential.:  ",A)')pr%selec_sigma
    Write(10,'("#  r_cutoff_sigma..:",1p,E16.8)')pr%r_cutoff_sigma
    Write(10,'("#  umax_sigma......:",1p,E16.8)')pr%umax_sigma
!
    Write(10,'("#  Pi potential....:  ",A)')pr%selec_pi
    Write(10,'("#  r_cutoff_pi.....:",1p,E16.8)')pr%r_cutoff_pi
    Write(10,'("#  umax_pi.........:",1p,E16.8)')pr%umax_pi
!
    If(pr%Lstate.Eq.'D')Then
      Write(10,'("#  Delta potential.:  ",A)')pr%selec_delta
      Write(10,'("#  r_cutoff_delta..:",1p,E16.8)')pr%r_cutoff_delta
      Write(10,'("#  umax_delta......:",1p,E16.8)')pr%umax_delta
    Endif        
    Write(10,'("#  xcm,ycm,zcm.:",1p,3E16.8)')pr%cm
    Write(10,'("#  Ekin........:",1p,E16.8)')pr%ekin
    Write(10,'("#  Elj.........:",1p,E16.8)')pr%elj
    Write(10,'("#  Ealphas.....:",1p,E16.8)')pr%ealphas
    Write(10,'("#  Ecorr.......:",1p,E16.8)')pr%ecor
    Write(10,'("#  Esolid......:",1p,E16.8)')pr%esolid
    Write(10,'("#  Ekinx.......:",1p,E16.8)')pr%ekinx
    Write(10,'("#  Eso.........:",1p,E16.8)')pr%eso
    Write(10,'("#  Eimpu.......:",1p,E16.8)')pr%evx
    Write(10,'("#  Etot........:",1p,E16.8)')pr%etot
    Write(10,'("#  <Lx,y,z>....:",1p,3E16.8)')pr%ang
    Write(10,'("#  <x2,y2,z2>..:",1p,3E16.8)')pr%r2
    write(10,'("#  ")')
    write(10,*) pr%xmax,pr%ymax,pr%zmax,pr%hx,pr%hy,pr%hz,pr%nx,pr%ny,pr%nz
    write(10,*) pr%rimp
    write(10,*) pr%vimp
    write(10,*) pr%ninvar
    write(10,*) pr%invar
    write(10,*) pr%psi
  close(10)

end
