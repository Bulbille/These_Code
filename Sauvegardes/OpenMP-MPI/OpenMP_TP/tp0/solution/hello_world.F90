program hello_world
!$ use OMP_LIB
  implicit none
  
  integer                      :: rang
!$ integer                     :: nb_taches

  rang = 0

  !$OMP PARALLEL PRIVATE(rang)
  !$OMP SINGLE
  !$ nb_taches = OMP_GET_NUM_THREADS()
  !$ print '(//,3X,"Execution de hello_world en parallele avec ",i2," threads")', nb_taches
  !$OMP END SINGLE

  !$ rang = OMP_GET_THREAD_NUM()
  print '(3X,"Bonjour depuis le thread de rang ",i2)', rang
  !$OMP END PARALLEL

end program hello_world
