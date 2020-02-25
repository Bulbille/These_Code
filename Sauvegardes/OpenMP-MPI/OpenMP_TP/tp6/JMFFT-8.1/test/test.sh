#                               -*- Mode: Sh -*- 
# test.sh --- Passage systematique de tous les tests Fourier
# 
# Auteur          : Jean-Marie Teuler (CNRS) <Jean-Marie.Teulet@lcp.u-psud.fr>
# Créé le         : Thu May 13 14:51:24 2004
# Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
# Dern. mod. le   : Thu May 13 14:52:39 2004
# 
#!/bin/ksh

# Quelques macros
F90=$(grep "FC" ../arch/Make.inc | awk -F"=" '{print $2}')
F90FLAGS=$(grep "FFLAGS" ../arch/Make.inc | awk -F "=" '{print $2}')
LD=$(grep "LD" ../arch/Make.inc | awk -F "=" '{print $2}')
LDFLAGS=$(grep "LDFLAGS" ../arch/Make.inc | awk -F "=" '{print $2}')

function testf {

  test -f stop.dat && exit 1

  /bin/rm -f temp1 temp2 temp.f90 temp.log a.out

  file=t$1.f90
  n=$2
  m=$3
  l=$4
  sign=$5
  scale=$6
  echo $*
  sed \
    -e "s/integer, parameter :: n =.*$/integer, parameter :: n = $n/" \
    -e "s/integer, parameter :: m =.*$/integer, parameter :: m = $m/" \
    -e "s/integer, parameter :: l =.*$/integer, parameter :: l = $l/" \
    -e "s/isign = 1/isign = $sign/" \
    -e "s/scale = .*$/scale = $scale/" \
    < $file > temp.f90
  $F90 -c $F90FLAGS temp.f90 > temp.log 2>&1
  if (( $? != 0 )) ; then
    echo Problemes de compilation
    exit 1
  fi
  $LD $LDFLAGS temp.o -L../lib -ljmfft >> temp.log 2>&1
  if (( $? != 0 )) ; then
    echo Problemes editions de liens
    exit 1
  fi
  ./a.out
  ./check > check.log
  grep -q OK check.log && echo Test OK ||
    { echo Problemes d execution; touch stop.dat; }
}

rm -f stop.dat
rm -f check check.log
$F90 $F90FLAGS check.f90 -o check $LDFLAGS > check.log 2>&1 ||
  { echo Problemes compilation check; exit 1; }

# 1. Test de jmccfft
test=jmccfft
testf $test   2 0 0 +1 2
testf $test   3 0 0 -1 2
testf $test   4 0 0 -1 2
testf $test   5 0 0 -1 2
testf $test   6 0 0 +1 2
testf $test   8 0 0 -1 2
testf $test   9 0 0 +1 4
testf $test  10 0 0 +1 4
testf $test  12 0 0 -1 2
testf $test  15 0 0 -1 2
testf $test  16 0 0 +1 2
testf $test  18 0 0 +1 2
testf $test  20 0 0 +1 2
testf $test  25 0 0 -1 5
testf $test  27 0 0 -1 5
testf $test  30 0 0 +1 3
testf $test  32 0 0 +1 3
testf $test  36 0 0 -1 6
testf $test  60 0 0 -1 6
testf $test  72 0 0 +1 7
testf $test  90 0 0 +1 7
testf $test 108 0 0 -1 8
testf $test    7 0 0 -1 8
testf $test   11 0 0 +1 2
testf $test   13 0 0 -1 2
testf $test   49 0 0 +1 8
testf $test  343 0 0 -1 8
testf $test   14 0 0 +1 8
testf $test   42 0 0 -1 8
testf $test  210 0 0 +1 8
testf $test   98 0 0 -1 8
testf $test  294 0 0 +1 8
testf $test 1470 0 0 -1 8

# 2. Test de jmccfftm
test=jmccfftm
testf $test   16   8 0 +1 4
testf $test    8  16 0 -1 3
testf $test   16   8 0 +1 4
testf $test   16  16 0 -1 5
testf $test   30   8 0 +1 2
testf $test   30  16 0 -1 3
testf $test    9   9 0 +1 6
testf $test   27   9 0 -1 7
testf $test    9  27 0 +1 8
testf $test   27  27 0 -1 9
testf $test   36  36 0 +1 8
testf $test   36  72 0 -1 7
testf $test   72  36 0 +1 6
testf $test   72  72 0 -1 5
testf $test   75  36 0 +1 6
testf $test   75  72 0 -1 5
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 3. Test de jmccfft2d
test=jmccfft2d
testf $test    5   8 0 +1 2
testf $test    5  16 0 -1 3
testf $test    8   8 0 +1 2
testf $test    8  16 0 -1 3
testf $test   16   8 0 +1 4
testf $test   16  16 0 -1 5
testf $test    9   9 0 +1 6
testf $test   27   9 0 -1 7
testf $test   30   9 0 -1 7
testf $test    9  30 0 -1 7
testf $test    9  27 0 +1 8
testf $test   27  27 0 -1 9
testf $test   36  36 0 +1 8
testf $test   36  72 0 -1 7
testf $test   72  36 0 +1 6
testf $test   72  72 0 -1 5
testf $test   60  50 0 +1 5
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 4. Test de jmccfft3d
test=jmccfft3d
testf $test     4  1  1 +1 4
testf $test     1  4  1 +1 4
testf $test     1  1  4 +1 4
testf $test     4  4  1 +1 4
testf $test     4  1  4 +1 4
testf $test     1  4  4 +1 4
testf $test     4  4  4 +1 4
testf $test     8  4  4 +1 2
testf $test     4  8  4 +1 2
testf $test     4  4  8 +1 2
testf $test     4  8  8 +1 2
testf $test     8  4  8 +1 2
testf $test     8  8  4 +1 2
testf $test     8  8  8 +1 2
testf $test    16  8  4 +1 2
testf $test     8  9 36 -1 3
testf $test    12 36  8 +1 4
testf $test    24 12 12 -1 5
testf $test    30 10  6 -1 5
testf $test    10  9 15 -1 5
testf $test     9 10 15 -1 5
testf $test     9 15 10 -1 5
testf $test    30  7  7 -1 5
testf $test     7 30  7 +1 5
testf $test    30 49  7 -1 5
testf $test    49 30  7 +1 5
testf $test    42 49  7 -1 5
testf $test    49 42  7 +1 5
testf $test    42 98  7 -1 5
testf $test    98 42  7 +1 5
testf $test    14 42 49 +1 5

# 5. Test de jmscfft
test=jmscfft
testf $test  2  0  0  +1  3
testf $test  4  0  0  +1  3
testf $test  6  0  0  +1  3
testf $test  8  0  0  +1  3
testf $test 12  0  0  +1  3
testf $test 16  0  0  +1  3
testf $test 18  0  0  +1  3
testf $test 24  0  0  +1  3
testf $test 30  0  0  +1  3
testf $test   14 0 0 -1 8
testf $test   22 0 0 +1 2
testf $test   26 0 0 -1 2
testf $test   98 0 0 +1 8
testf $test  686 0 0 -1 8
testf $test   42 0 0 -1 8
testf $test  210 0 0 +1 8
testf $test   98 0 0 -1 8
testf $test  294 0 0 +1 8
testf $test 1470 0 0 -1 8

# 6. Test de jmscfftm
test=jmscfftm
testf $test     2  1 0 +1 3
testf $test     2  3 0 +1 3
testf $test     2  2 0 +1 3
testf $test     8  8 0 +1 3
testf $test     8 16 0 -1 3
testf $test    16  8 0 +1 3
testf $test     8  6 0 +1 3
testf $test     8 12 0 -1 3
testf $test    16 18 0 +1 3
testf $test     3 16 0 +1 3
testf $test     9 16 0 +1 3
testf $test    27 16 0 +1 3
testf $test     4  1 0 +1 3
testf $test     8  3 0 -1 2
testf $test    16  5 0 +1 1
testf $test    30 30 0 +1 1
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 7. Test de jmscfft2d
test=jmscfft2d
testf $test     2  1 0 +1 3
testf $test     1  2 0 -1 3
testf $test     2  2 0 +1 3
testf $test     3  2 0 +1 3
testf $test     3  4 0 +1 3
testf $test    30 20 0 +1 3
testf $test     2  3 0 +1 3
testf $test     3  4 0 +1 3
testf $test     4  3 0 +1 3
testf $test     4  4 0 -1 3
testf $test     4  8 0 +1 3
testf $test     8  4 0 -1 3
testf $test     8  8 0 +1 3
testf $test     8 16 0 -1 3
testf $test    16  8 0 +1 3
testf $test    30 30 0 -1 3
testf $test    27 20 0 -1 3
testf $test    20 27 0 +1 3
testf $test    27 20 0 +1 3
testf $test    20 27 0 -1 3
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 8. Test de jmscfft3d
test=jmscfft3d
testf $test     8  8  8 +1 3
testf $test     9  8  8 -1 3
testf $test     8  9  8 +1 3
testf $test     8  8  9 -1 3
testf $test     8  9  9 +1 3
testf $test     9  8  9 -1 3
testf $test     9  9  8 +1 3
testf $test    16  8  8 -1 3
testf $test     8 16  8 +1 3
testf $test     8  8 16 -1 3
testf $test    16 16  8 +1 3
testf $test    16  8 16 -1 3
testf $test     8 16 16 +1 3
testf $test    30 10  6 +1 3
testf $test    30  7  7 -1 5
testf $test     7 30  7 +1 5
testf $test    30 49  7 -1 5
testf $test    49 30  7 +1 5
testf $test    42 49  7 -1 5
testf $test    49 42  7 +1 5
testf $test    42 98  7 -1 5
testf $test    98 42  7 +1 5
testf $test    14 42 49 +1 5

# 9. Test de jmcsfft
test=jmcsfft
testf $test     2  0 0 -1 3
testf $test     4  0 0 +1 2
testf $test     8  0 0 -1 2
testf $test    16  0 0 -1 2
testf $test    32  0 0 -1 2
testf $test     6  0 0 -1 2
testf $test    12  0 0 +1 2
testf $test    18  0 0 +1 2
testf $test    30  0 0 +1 2
testf $test   14 0 0 +1 8
testf $test   42 0 0 -1 8
testf $test  210 0 0 +1 8
testf $test   98 0 0 -1 8
testf $test  294 0 0 +1 8
testf $test 1470 0 0 -1 8

# 9bis. Test de jmcsfft-safe
test=jmcsfft-safe
testf $test     2  0 0 -1 3
testf $test     4  0 0 +1 2
testf $test     8  0 0 -1 2
testf $test    16  0 0 -1 2
testf $test    32  0 0 -1 2
testf $test     6  0 0 -1 2
testf $test    12  0 0 +1 2
testf $test    18  0 0 +1 2
testf $test    30  0 0 +1 2
testf $test   14 0 0 +1 8
testf $test   42 0 0 -1 8
testf $test  210 0 0 +1 8
testf $test   98 0 0 -1 8
testf $test  294 0 0 +1 8
testf $test 1470 0 0 -1 8

# 10. Test de jmcsfftm
test=jmcsfftm
testf $test    48  1 0 +1 3
testf $test     2  1 0 +1 3
testf $test     4  1 0 -1 5
testf $test     6  1 0 +1 3
testf $test     8  1 0 +1 3
testf $test    18  1 0 +1 3
testf $test    16  1 0 +1 3
testf $test    32  1 0 +1 3
testf $test     2  2 0 -1 3
testf $test     3  4 0 +1 2
testf $test     4  8 0 -1 2
testf $test     6 16 0 +1 2
testf $test    12  8 0 -1 2
testf $test    18 10 0 +1 2
testf $test    16  8 0 -1 2
testf $test    27  8 0 +1 2
testf $test    30 30 0 +1 2
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 10bis. Test de jmcsfftm-safe
test=jmcsfftm-safe
testf $test    48  1 0 +1 3
testf $test     2  1 0 +1 3
testf $test     4  1 0 -1 5
testf $test     6  1 0 +1 3
testf $test     8  1 0 +1 3
testf $test    18  1 0 +1 3
testf $test    16  1 0 +1 3
testf $test    32  1 0 +1 3
testf $test     2  2 0 -1 3
testf $test     3  4 0 +1 2
testf $test     4  8 0 -1 2
testf $test     6 16 0 +1 2
testf $test    12  8 0 -1 2
testf $test    18 10 0 +1 2
testf $test    16  8 0 -1 2
testf $test    27  8 0 +1 2
testf $test    30 30 0 +1 2
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 11. Test de jmcsfft2d
test=jmcsfft2d
testf $test   2  2 0 +1 3
testf $test   2  3 0 -1 3
testf $test   3  2 0 +1 3
testf $test   4  4 0 -1 3
testf $test   6  4 0 +1 3
testf $test   6  6 0 -1 3
testf $test   8  8 0 +1 3
testf $test  12  8 0 -1 3
testf $test   8 12 0 +1 3
testf $test  12 12 0 -1 3
testf $test  30 30 0 -1 3
testf $test  20 27 0 -1 3
testf $test  27 20 0 +1 3
testf $test  20 27 0 +1 3
testf $test  27 20 0 -1 3
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 11bis. Test de jmcsfft2d-safe
test=jmcsfft2d-safe
testf $test   2  2 0 +1 3
testf $test   2  3 0 -1 3
testf $test   3  2 0 +1 3
testf $test   4  4 0 -1 3
testf $test   6  4 0 +1 3
testf $test   6  6 0 -1 3
testf $test   8  8 0 +1 3
testf $test  12  8 0 -1 3
testf $test   8 12 0 +1 3
testf $test  12 12 0 -1 3
testf $test  30 30 0 -1 3
testf $test  20 27 0 -1 3
testf $test  27 20 0 +1 3
testf $test  20 27 0 +1 3
testf $test  27 20 0 -1 3
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 12. Test de jmcsfft3d
test=jmcsfft3d
testf $test   2  2  2 +1 3
testf $test   2  3  2 -1 2
testf $test   3  2  2 -1 1
testf $test   2  2  3 +1 1
testf $test   4  4  4 -1 2
testf $test   6  4  4 +1 3
testf $test   4  6  4 -1 4
testf $test   4  4  6 +1 5
testf $test   6  6  4 -1 6
testf $test   6  4  6 +1 7
testf $test   4  6  6 -1 8
testf $test   6  9  9 +1 3
testf $test   9  6  9 -1 9
testf $test   9  9  6 +1 5
testf $test   6  6  9 -1 6
testf $test   6  9  6 +1 7
testf $test   9  6  6 -1 8
testf $test   6  6  6 +1 9
testf $test  12 12 18 -1 8
testf $test  12 18 12 +1 7
testf $test  18 12 12 -1 6
testf $test  12 18 18 +1 5
testf $test  18 18 12 -1 4
testf $test  18 18 12 +1 3
testf $test  18 18 18 -1 2
testf $test  30 10  6 -1 2
testf $test    30  7  7 -1 5
testf $test     7 30  7 +1 5
testf $test    30 49  7 -1 5
testf $test    49 30  7 +1 5
testf $test    42 49  7 -1 5
testf $test    49 42  7 +1 5
testf $test    42 98  7 -1 5
testf $test    98 42  7 +1 5
testf $test    14 42 49 +1 5

# 12bis. Test de jmcsfft3d-safe
test=jmcsfft3d-safe
testf $test   2  2  2 +1 3
testf $test   2  3  2 -1 2
testf $test   3  2  2 -1 1
testf $test   2  2  3 +1 1
testf $test   4  4  4 -1 2
testf $test   6  4  4 +1 3
testf $test   4  6  4 -1 4
testf $test   4  4  6 +1 5
testf $test   6  6  4 -1 6
testf $test   6  4  6 +1 7
testf $test   4  6  6 -1 8
testf $test   6  9  9 +1 3
testf $test   9  6  9 -1 9
testf $test   9  9  6 +1 5
testf $test   6  6  9 -1 6
testf $test   6  9  6 +1 7
testf $test   9  6  6 -1 8
testf $test   6  6  6 +1 9
testf $test  12 12 18 -1 8
testf $test  12 18 12 +1 7
testf $test  18 12 12 -1 6
testf $test  12 18 18 +1 5
testf $test  18 18 12 -1 4
testf $test  18 18 12 +1 3
testf $test  18 18 18 -1 2
testf $test  30 10  6 -1 2
testf $test    30  7  7 -1 5
testf $test     7 30  7 +1 5
testf $test    30 49  7 -1 5
testf $test    49 30  7 +1 5
testf $test    42 49  7 -1 5
testf $test    49 42  7 +1 5
testf $test    42 98  7 -1 5
testf $test    98 42  7 +1 5
testf $test    14 42 49 +1 5

# 13. Test de jmrfftmlt, sens C -> R
test=jmrfftmlt1
testf $test    4  2 0 +1 1
testf $test   10 20 0 +1 1
testf $test   20 20 0 +1 1
testf $test   30   7 0 -1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 13bis. Test de jmrfftmlt-safe, sens C -> R
test=jmrfftmlt1-safe
testf $test    4  2 0 +1 1
testf $test   10 20 0 +1 1
testf $test   20 20 0 +1 1
testf $test   30   7 0 -1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 14. Test de jmrfftmlt, sens R -> C
test=jmrfftmlt2
testf $test    4  2 0 -1 1
testf $test   10 20 0 -1 1
testf $test   20 10 0 -1 1
testf $test   30 10 0 -1 1
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 15. Test de jmcfftmlt
test=jmcfftmlt
testf $test   10 20 0 +1 1
testf $test   10 40 0 -1 1
testf $test   10 30 0 +1 1
testf $test   20 20 0 -1 1
testf $test   20 40 0 +1 1
testf $test   20 30 0 -1 1
testf $test    4  2 0 +1 1
testf $test    3  2 0 -1 1
testf $test   15 40 0 +1 1
testf $test   15 30 0 -1 1
testf $test   25 20 0 +1 1
testf $test   25 40 0 -1 1
testf $test   25 30 0 +1 1
testf $test   20 15 0 -1 1
testf $test   20 45 0 +1 1
testf $test   20 25 0 -1 1
testf $test   20 75 0 +1 1
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 16. Test de jmsinfftmlt
test=jmsinfftmlt
testf $test   10 20 0 +1 1
testf $test   10 40 0 -1 1
testf $test   10 30 0 +1 1
testf $test   20 20 0 -1 1
testf $test   20 40 0 +1 1
testf $test   20 30 0 -1 1
testf $test    4  2 0 +1 1
testf $test    3  2 0 -1 1
testf $test   15 40 0 +1 1
testf $test   15 30 0 -1 1
testf $test   25 20 0 +1 1
testf $test   25 40 0 -1 1
testf $test   25 30 0 +1 1
testf $test   20 15 0 -1 1
testf $test   20 45 0 +1 1
testf $test   20 25 0 -1 1
testf $test   20 75 0 +1 1
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5

# 17. Test de jmcosfftmlt
test=jmcosfftmlt
testf $test   10 20 0 +1 1
testf $test   10 40 0 -1 1
testf $test   10 30 0 +1 1
testf $test   20 20 0 -1 1
testf $test   20 40 0 +1 1
testf $test   20 30 0 -1 1
testf $test    4  2 0 +1 1
testf $test    3  2 0 -1 1
testf $test   15 40 0 +1 1
testf $test   15 30 0 -1 1
testf $test   25 20 0 +1 1
testf $test   25 40 0 -1 1
testf $test   25 30 0 +1 1
testf $test   20 15 0 -1 1
testf $test   20 45 0 +1 1
testf $test   20 25 0 -1 1
testf $test   20 75 0 +1 1
testf $test   30   7 0 -1 5
testf $test    7  30 0 +1 5
testf $test   30  49 0 -1 5
testf $test   49  30 0 +1 5
testf $test   42  49 0 -1 5
testf $test   49  42 0 +1 5
testf $test   42  98 0 -1 5
testf $test   98  42 0 +1 5
