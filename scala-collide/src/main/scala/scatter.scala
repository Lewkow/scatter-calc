
package scatter

import potential.Potential
import cross_section.CrossSection
import mongo.Mongo
import scatter_math.ScatterMath
import scala.collection.mutable.ListBuffer

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import org.apache.spark.sql._
import org.apache.spark.rdd._

class Scatter(energyC: Double, 
              projectileC: String, 
              targetC: String,
              potential_typeC: String) extends Serializable {

  import math.sqrt
  val energy: Double = energyC
  val coll = new Collision(projectileC, targetC)
  val projectile: String = projectileC
  val target: String = targetC
  val potential_type: String = potential_typeC
  val phase = new Phaseshift(energy, potential_type, coll)
  val collision_id: String = projectile + "_" + target + "_" +
                             potential_type
  // Get some cross sections with the phases
  var cross_sections = new CrossSection(collision_id, phase.phases, coll.wavenumber(energy), energy)
  val tcs = cross_sections.total_cross_section()

  def get_dcs(theta: Double): Double = {
    val dcs: Double = cross_sections.differential_cross_section(theta)
    dcs
  }

  def printValues() {
    val grid_start = phase.grid.start_pos 
    val grid_end = phase.grid.end_pos 
    coll.printValues()
    println(f"Energy:       $energy%1.3f [eV]")
    println(f"TCS:          $tcs%1.3f [a0^2]")

  }

}

class Collision(projectile: String,
                target: String) extends Serializable {

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

class Phaseshift(energy: Double, potential_type: String, coll: Collision) extends Serializable {
  val MAX_GRID_SIZE = 100000
  val pot = new Potential(potential_type)
  val grid = new Grid(energy, pot, coll)
  val phases = get_phases()
  // val phases = get_spark_phases()

  def get_spark_phases(): Array[Double] = {
    val conf = new SparkConf().setAppName("CIS Views").setMaster("spark://MA01322RFH0:7077")
    val sc = new SparkContext(conf)
    val small_phase: Double = 5.0e-7
    var keep_going: Boolean = true
    var converge_counter: Int = 0
    val max_converge_counter: Int = 30
    var partial_wave: Array[Int] = 0 to 1000 toArray
    val r0 = sc.parallelize(partial_wave)
    val r1 = r0.map(x => numerov(x))
    val p = r1.collect
    val phase = p.toArray
    for (i <- 0 to phase.length-1) {println(phase(i))}
    println("energy: " + energy.toString + " phaseshift: " + phase.toString)
    phase
  }

  def get_phases(): Array[Double] = {
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
    // println("energy: " + energy.toString + " phaseshift: " + phase.toString)
    phase.toArray
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
    if (ph.isNaN) {
      0.0.toDouble
    } else {
      ph
    }
  }
}

class Grid(energyC: Double, pot: Potential, coll: Collision) extends Serializable {
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
