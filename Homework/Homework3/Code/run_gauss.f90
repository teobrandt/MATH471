program run_gauss

  implicit none
  call lglnodes(x,w,n)
  f = exp(cos(pi*pi*x))
  Integral_value = sum(f*w)
  write *, Integral_value

end program run_gauss
