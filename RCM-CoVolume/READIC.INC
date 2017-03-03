      OPEN(99, file='init.dat')
      READ(99,*) TUBLEN
      READ(99,*) M
      READ(99,*) NOTIST
      READ(99,*) NOPROF
      READ(99,*) TOL
      READ(99,*) CFLCOE
      READ(99,*) GAMMA
      READ(99,*) B
      READ(99,*) IC_DL0
      READ(99,*) IC_UL0
      READ(99,*) IC_PL0
      READ(99,*) IC_DR0
      READ(99,*) IC_UR0
      READ(99,*) IC_PR0
      READ(99,*) IC_X0
      CLOSE(99)
      