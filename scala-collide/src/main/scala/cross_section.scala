package cross_section

import scatter_math.ScatterMath
import scala.collection.mutable.ListBuffer

class CrossSection(phases: Array[Double], wavenumber: Double) {

  def test_differential_cross_section() {
    val N: Int = 100
    val dtheta: Double = math.Pi/N.toDouble
    var i: Int = 1
    for (i <- 0 to N) {
      val theta: Double = i.toDouble*dtheta
      val dcs: Double = differential_cross_section(theta)
      println(f"Theta: $theta%1.2f   DCS: $dcs%1.2f")
    }
  }

  def differential_cross_section(theta: Double): Double = {
    var reAmp: Double = 0.0d
    var imAmp: Double = 0.0d
    var i: Int = 0
    val LMax: Int = phases.length-1
    var sm = new ScatterMath()
    val leg = sm.legendre(math.cos(theta),LMax)
    for (i <- 0 to LMax) {
      var coeff: Double = (2.0d*i.toDouble+1.0d)*leg(i)*math.sin(phases(i))
      reAmp = reAmp + coeff*math.cos(phases(i))
      imAmp = imAmp + coeff*math.sin(phases(i))
    }
    (1.0d/(wavenumber*wavenumber))*(reAmp*reAmp+imAmp*imAmp).toDouble
  }

  def total_cross_section(): Double = {
    var i: Int = 0
    var tcs: Double = 0.0d

    for (i <- 0 to phases.length-1) {
      var sin_phase = math.sin(phases(i))
      tcs += (2.0d*i.toDouble+1.0d)*sin_phase*sin_phase
    } 
    tcs*(4.0d*math.Pi/(wavenumber*wavenumber))
  }



}