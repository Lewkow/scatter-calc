
package scatter

import potential.Potential
// import mongo.Mongo
import scala.collection.mutable.ListBuffer


class Scatter(energyC: Double, 
              projectileC: String, 
              targetC: String,
              potential_typeC: String) {

  import math.sqrt
  val energy: Double = energyC
  val coll = new Collision(projectileC, targetC)
  val projectile: String = projectileC
  val target: String = targetC
  val potential_type: String = potential_typeC
  val phase = new Phaseshift(energy, potential_type, coll)
  // Get some cross sections with the phases
  var cross_sections = new CrossSection(phase.phases, coll.wavenumber(energy))
  val tcs = cross_sections.total_cross_section()

  def printValues() {
    val grid_start = phase.grid.start_pos 
    val grid_end = phase.grid.end_pos 
    coll.printValues()
    println(f"Energy:       $energy%1.3f [eV]")
    println(f"TCS:          $tcs%1.3f [a0^2]")

  }

}

class Collision(projectile: String,
                target: String) {

  val masses = Map("H"->1.00782503207d,
                   "He3"->3.0160293191d,
                   "He4"->4.00260325415d,
                   "O"->15.99491461956d)
  val projectile_mass: Double = masses(projectile)
  val target_mass: Double = masses(target)
  val reduced_mass: Double = reducedMass(projectile_mass,target_mass)

  def reducedMass(proj: Double, targ: Double): Double = {
    val denom: Double = proj + targ
    val reduced_mass: Double = proj*targ/denom
    reduced_mass
  }

  def wavenumber(energy: Double): Double = {
    math.sqrt(2.0d*reduced_mass*energy)
  }

  def printLine() {
    println("--------------------------------")
  }

  def printValues() {
    printLine()
    println(s"Projectile:   $projectile") 
    println(f"Mass:         $projectile_mass%1.3f [u]")
    printLine()
    println(s"Target:       $target")
    println(f"Mass:         $target_mass%1.3f [u]")
    printLine()
    println(f"Reduced Mass: $reduced_mass%1.3f [u]")
    printLine()
  }

}

class CrossSection(phases: ListBuffer[Double], wavenumber: Double) {

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
    (1.0d/(wavenumber*wavenumber))*(reAmp*reAmp+imAmp*imAmp)
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

class Phaseshift(energy: Double, potential_type: String, coll: Collision) {
  val MAX_GRID_SIZE = 100000
  val pot = new Potential(potential_type)
  val grid = new Grid(energy, pot, coll)
  var phases = new ListBuffer[Double]
  phases = get_phases()

  def get_phases() = {
    val small_phase: Double = 5.0e-7
    var partial_wave: Int = 0
    var keep_going: Boolean = true
    var converge_counter: Int = 0
    val max_converge_counter: Int = 30
    var phase = new ListBuffer[Double]

    while (keep_going == true) {
      var current_phase = numerov(partial_wave)
      // println(s"L: $partial_wave phase: $current_phase")
      phase += current_phase
      if (math.abs(current_phase) < small_phase) {
        converge_counter = converge_counter + 1
      }
      if (converge_counter > max_converge_counter) {
        keep_going = false
      }
      partial_wave += 1
    }
    phase
  }

  def dummy_numerov(partial_wave: Int) = {
    // dummy numerov function for now
    val p: Double = (partial_wave + 1).toDouble
    // println(p, 1.0d/p)
    1.0d/(p*p)
  }

  def numerov(partial_wave: Int): Double = {
    // dummy_numerov(partial_wave)
    val big: Double = 1.0e20
    val small: Double = 1.0e-20
    val k: Double = coll.wavenumber(energy)
    val k2: Double = k*k
    val l2: Double = (partial_wave*(partial_wave+1)).toDouble
    val hh: Double = (grid.dr_grid*grid.dr_grid)/12.0d
    val effective_pot_grid: Array[Double] = grid.build_effective_pot_grid(k2, l2)
    var wave = new Array[Double](grid.N_grid+2)
    wave(0) = 0.0d
    wave(1) = 0.1d
    var i: Int = 1
    for (i <- 1 to grid.N_grid-2) {
      if (math.abs(wave(i)) > big) {
        var j: Int = 0
        for (j <- 0 to i) {
          wave(j) = wave(j)*small
        }
      }
      if (math.abs(wave(i)) < small) {
        var j: Int = 0
        for (j <- 0 to i) {
          wave(j) = wave(j)*big
        }
      }
      // println(hh, effective_pot_grid(i+1), effective_pot_grid(i), effective_pot_grid(i-1), wave(i), wave(i-1))
      wave(i+1) = ( (2.0d-10.0d*hh*effective_pot_grid(i)*wave(i) - 
        (1.0d+hh*effective_pot_grid(i-1))*wave(i-1)) ) / 
        (1.0d+hh*effective_pot_grid(i+1))
    }
    val d1: Double = (wave(grid.N_grid-3)-wave(grid.N_grid-5))/(2.0d*grid.dr_grid)
    val d2: Double = (wave(grid.N_grid-2)-wave(grid.N_grid-6))/(4.0d*grid.dr_grid)
    val d4: Double = (wave(grid.N_grid)  -wave(grid.N_grid-8))/(8.0d*grid.dr_grid)
    val d12: Double = d1 + (d1-d2)/3.0d 
    val d24: Double = d2 + (d2-d4)/3.0d 
    val grad: Double = d12 + (d12-d24)/15.0d
    val kr: Double = k*grid.pos_grid(grid.N_grid-4)
    // println(d1, d2, d4, d12, d24)
    // println(grad, kr)
    val bessel = new ScatterMath()
    var sine: Double = bessel.spherical_bessel_1st(kr, partial_wave)
    var cose: Double = bessel.spherical_bessel_2nd(kr, partial_wave)
    var dsine: Double = bessel.dx_spherical_bessel_1st(kr, partial_wave)
    var dcose: Double = bessel.dx_spherical_bessel_2nd(kr, partial_wave)
    // println(sine, cose, dsine, dcose)
    dsine = dsine*k
    dcose = dcose*k
    // println(dsine, dcose)
    val s: Double = (dsine*wave(grid.N_grid-4) - sine*grad)/k
    val c: Double = (dcose*wave(grid.N_grid-4) - cose*grad)/k
    val ph: Double = math.atan(s/c)
    // println(s, c, ph)
    ph
  }
}

class Grid(energyC: Double, pot: Potential, coll: Collision) {
  val energy: Double = energyC
  val start_pos: Double = grid_start_finder() 
  val end_pos: Double = grid_end_finder()
  val dr_grid: Double = grid_dr_finder()
  val N_grid: Int = grid_N_finder() 
  val pos_grid: Array[Double] = build_pos_grid()
  val pot_grid: Array[Double] = build_pot_grid()
  // printValues()

  def printValues() {
    println(f"Grid Start:   $start_pos%1.3f [a0]")
    println(f"Grid End:     $end_pos%1.3f [a0]")
  }

  def grid_dr_finder(): Double = {
    val k_max: Double = math.sqrt(2.0d*coll.reduced_mass*(energy+pot.depth))
    val h_max: Double = (1.0d/k_max)/20.0d
    var step: Double = 7.0e-4
    if (7.0e-4 > h_max) {
      step = h_max
    }
    step
  }  

  def grid_N_finder() = {
    val N_grid: Int = ((end_pos - start_pos)/dr_grid).toInt + 1
    N_grid
  }

  def grid_end_finder() = {
    val lim: Double = 1e-8
    val dr: Double = 0.1d
    var r: Double = start_pos
    var pot1: Double = pot.get_potential(1.0d)
    var pot2: Double = pot1
    var dpot: Double = 0.0d
    var going : Boolean = true
    var i: Int = 1
    var r_end: Double = start_pos

    while (going == true) {
      r = start_pos + i*dr
      pot2 = pot.get_potential(r)
      dpot = math.abs(pot2-pot1)
      if (dpot < lim) {
        r_end = r
        going = false
      }
      pot1 = pot2
      i = i + 1
    }
    r_end
  }

  def grid_start_finder() = {
    var i: Int = 0
    var j: Int = 0 
    val start: Double = 10.0d
    val dr: Double = 0.0001
    val lim: Double = energy*1.3d
    var going: Boolean = true
    var r: Double = start

    while (going == true) {
      r = start - j.toDouble*dr
      var pot_r: Double = pot.get_potential(r)
      if ( (pot_r >= lim) | (r <= 1.0e-5) ) {
        going = false 
      } else {
        j = j + 1
      }
    }
    r
  }

  def build_pos_grid():Array[Double] = {
    val z = new Array[Double](N_grid)    
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      z(i) = start_pos + i.toDouble*dr_grid
    }
    z 
  }

  def build_pot_grid():Array[Double] = {
    val z = new Array[Double](N_grid)
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      var x: Double = start_pos + i.toDouble*dr_grid
      z(i) = -2.0d*coll.reduced_mass*pot.get_potential(x)
    }   
    z
  }

  def build_effective_pot_grid(k2: Double, l2: Double):Array[Double] = {
    val z = new Array[Double](N_grid)
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      z(i) = pot_grid(i) + k2 - l2/(pos_grid(i)*pos_grid(i))
    }
    z
  }

}

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