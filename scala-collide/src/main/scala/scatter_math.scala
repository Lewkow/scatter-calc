package scatter_math

class ScatterMath() {

  def legendre(x: Double, n: Int): Array[Double] = {
    val L: Int = n+1
    val leg = new Array[Double](L)
    leg(0) = 1.0d
    leg(1) = x
    var i: Int = 2
    for (i <- 1 to (L-2)) {
      var nn: Double = (i-1).toDouble
      leg(i+1) = ( (2.0d*nn+1.0d) / (nn+1.0d) )*x*leg(i) - ( nn / (nn+1.0d) )*leg(i-1)
    } 
    leg
  }
 
  def spherical_bessel_1st(x: Double, j: Int) = {
    val IAccuracy: Int = 20 
    val big: Double = 1.0e20
    val jbig: Int = j + (math.sqrt((j*IAccuracy).toDouble)).toInt
    var bes: Double = 0.0d
    var temp: Double = 0.0d
    if (j == 0) {
      bes = math.sin(x)
    } else if (j == 1) {
      bes = math.sin(x)/x - math.cos(x)
    } else {
      if (x > (j.toDouble+0.5d)) {
        // use upward recurrance
        var fm1: Double = math.sin(x)/x - math.cos(x)
        var fm2: Double = math.sin(x)
        var N: Int = 2
        var f: Double = 0.0d
        for (N <- 2 to j) {
          f = (2.0d*N.toDouble-1.0d)*fm1/x - fm2
          fm2 = fm1
          fm1 = f
        }
        bes = f
      } else {
        // use downward recurrance
        var fp1: Double = 1.0e-20
        var fp2: Double = 1.0e-20
        var f: Double = 0.0d
        var N: Int = jbig
        while (N >= jbig) {
          f = (2.0d*N.toDouble+3.0d)*fp1/x - fp2
          if (N == j) {
            temp = f 
          }
          if (math.abs(f) > big) {
            fp2 = fp2/big
            fp1 = fp1/big
            f   = f/big
            if (N < j) {
              temp = temp/big
            }
          }
          N = N - 1
          fp2 = fp1
          fp1 = f
        }
        bes = temp*math.sin(x)/f 
      }
    }
    bes
  }

  def dx_spherical_bessel_1st(x: Double, j: Int) = {
    var bes: Double = 0.0d 
    if (j == 0) {
      bes = math.cos(x)
    } else {
      val pres: Double = spherical_bessel_1st(x,j)
      val futr: Double = spherical_bessel_1st(x,j+1)
      val past: Double = spherical_bessel_1st(x,j-1)
      val bes: Double = ( pres/x - (futr - past) ) / 2.0d
    }
    bes
  }

  def spherical_bessel_2nd(x: Double, j: Int): Double = {
    var bes: Double = 0.0d
    if (j == 0) {
      bes = -math.cos(x)
    } else if (j == 1) {
      bes = -math.cos(x)/x - math.sin(x)
    } else {
      var fm2: Double = -math.cos(x)
      var fm1: Double = -math.cos(x)/x - math.sin(x)
      var f: Double = 0.0d
      var N: Int = 2
      for (N <- 2 to j) {
        f = (2.0d*N.toDouble-1.0d)*fm1/x - fm2
        fm2 = fm1
        fm1 = f
      } 
      bes = f 
    }
    bes
  }

  def dx_spherical_bessel_2nd(x:Double, j: Int): Double = {
    var bes: Double = 0.0d
    if (j == 0) {
      bes = math.sin(x)
    } else {
      var pres: Double = spherical_bessel_2nd(x, j)
      var futr: Double = spherical_bessel_2nd(x, j+1)
      var past: Double = spherical_bessel_2nd(x, j-1)
      bes = (pres/x-(futr-past))/2.0d 
    }
    bes
  }

}
